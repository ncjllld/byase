# This file is part of BYASE.
#
# BYASE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BYASE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BYASE.  If not, see <https://www.gnu.org/licenses/>.
#
# Author: Lili Dong
#

from typing import Tuple

import numpy as np
import pandas as pd
import scipy.stats
import HTSeq
import pymc3 as pm
from pymc3.backends.tracetab import trace_to_dataframe
import theano.tensor as tt

from ..message import MessageCenter
from .task import Task
from .align import BAMParam, TaskAlign, AlignStat


class TaskInferenceError(Exception):
    """Task inference error."""
    def __init__(self, task_id, msg):
        super().__init__(msg)
        self.task_id = task_id


class _ObservedDistribution(pm.Continuous):
    """Distribution of observed data."""

    def __init__(self, _p, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.p = _p

    def logp(self, value):
        return tt.log(tt.mul(tt.dot(value, self.p)))


class TaskInference(TaskAlign):
    """Task inference."""

    def __init__(self, task: Task, bam_param: BAMParam, mcmc_samples: int, tune_samples: int, mc: MessageCenter):
        super().__init__(task, bam_param, mc)
        self.mcmc_samples = mcmc_samples
        self.tune_samples = tune_samples

    @property
    def isoform_lens(self):
        """An ndarray containing the length of each isoform."""
        return np.array([iso.length for iso in self.segment.isoforms])

    def _calc_compatible_factor(self, align_stat: AlignStat) -> np.ndarray:
        """"Calculate compatible factor."""

        # Based on compatible alleles.
        tmp = np.repeat([0], self.ploidy)
        for i in align_stat.compatible_alleles:
            tmp[i] = 1
        compatible_factor = np.repeat(tmp, self.isoforms_count)

        # Based on compatible isoforms.
        tmp = np.repeat([0], self.isoforms_count)
        for i in align_stat.compatible_isoforms:
            tmp[i] = 1
        tmp = np.tile(tmp, self.ploidy)
        compatible_factor &= tmp

        return compatible_factor

    def _transform_se_observed_data(self, aln_stat: AlignStat) -> np.ndarray:
        """Transform single-end alignment to observed data."""
        assert aln_stat.aln2 is None
        cf = self._calc_compatible_factor(aln_stat)

        mappable_reads_count = self.isoform_lens - self.bam_param.read_len + 1
        c = mappable_reads_count.astype(float)
        c[c <= 0] = np.inf
        c = 1 / c

        data = cf * np.tile(c, self.ploidy)
        return data

    def _transform_pe_observed_data(self, aln_stat: AlignStat) -> np.ndarray:
        """Transform paired-end alignments pair to observed data."""
        cf = self._calc_compatible_factor(aln_stat)

        aln1, aln2 = aln_stat.aln, aln_stat.aln2
        start = min(aln1.iv.start, aln2.iv.start)
        end = max(aln1.iv.end, aln2.iv.end)
        fragment_iv = HTSeq.GenomicInterval(self.segment.iv.chrom, start, end, '.')

        # Inferred insert sizes for each isoform.
        inferred_insert_sizes = np.repeat([0], self.isoforms_count)
        for iv, value in self.region[fragment_iv].steps():
            isoform_nums = [v for v in value if isinstance(v, int)]
            for i in isoform_nums:
                inferred_insert_sizes[i] += (iv.end - iv.start)

        iis = inferred_insert_sizes.astype(float)
        iis[iis == 0] = -np.inf

        # Mappable fragments count.
        c = self.isoform_lens.astype(float) - iis + 1
        c[c <= 0] = np.inf
        p = scipy.stats.norm.pdf(iis, self.bam_param.insert_size_mean, self.bam_param.insert_size_std)
        c = 1 / c * p

        data = cf * np.tile(c, self.ploidy)
        return data

    def _extract_observed_data(self):
        """Extract observed data."""
        data = []
        for aln_stat in self.align_stats_iterator():
            if self.bam_param.paired_end:
                d = self._transform_pe_observed_data(aln_stat)
            else:
                d = self._transform_se_observed_data(aln_stat)

            assert len(d[d < 0]) == 0
            if len(d[d > 0]) == 0:
                continue

            data.append(d)

        if len(data) == 0:
            raise TaskInferenceError(self.id, 'No observed data found.')

        return np.array(data)

    def run(self) -> Tuple[int, pd.DataFrame]:
        """Bayesian inference."""
        observed_data = self._extract_observed_data()
        self.mc.log_debug('{} reads (read pairs) are used for Bayesian inference.'.format(len(observed_data)))

        p_count = self.isoforms_count * self.ploidy

        with pm.Model():
            p = pm.Dirichlet(name='p', a=tt.stack([1 for _ in range(p_count)]), shape=p_count)
            c = np.tile(self.isoform_lens, self.ploidy)
            p_rescaled = pm.Deterministic(name='p_rescaled', var=(p * c) / tt.dot(p, c))
            _ObservedDistribution(name='observed', _p=p_rescaled, observed=observed_data)

            # Inference.
            trace = pm.sample(self.mcmc_samples, tune=self.tune_samples, chains=1, progressbar=False)

        # Convert trace to data frame.
        trace = trace_to_dataframe(trace)

        # Remove rescaled variables.
        cols = [col for col in trace.columns if not col.startswith('p_rescaled')]
        assert len(cols) == p_count
        trace = trace[cols]

        # Rename trace columns.
        rename_map = {}
        for i in range(self.isoforms_count):
            for j in range(self.ploidy):
                n = i + j * self.isoforms_count
                rename_map['p__{}'.format(n)] = 'RAW_I{}_A{}'.format(i, j)
        trace = trace.rename(columns=rename_map)

        return len(observed_data), trace
