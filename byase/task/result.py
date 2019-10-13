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

from typing import Optional

import numpy as np
import pandas as pd

from .task import Task


class TaskResultMeta:
    """Task result meta."""

    @staticmethod
    def _get_var_name(var_type: str, iso_num: Optional[int] = None, allele_num: Optional[int] = None,
                      allele_num2: Optional[int] = None):
        """Get variable name."""
        assert (iso_num is not None) or (allele_num is not None)

        if var_type == 'expression':
            var_name = 'EXP'
        else:
            assert var_type == 'raw'
            assert (iso_num is not None) and (allele_num is not None) and (allele_num2 is None)
            var_name = 'RAW'

        if iso_num is not None:
            var_name += '_I{}'.format(iso_num)
        if allele_num is not None:
            var_name += '_A{}'.format(allele_num)
        if allele_num2 is not None:
            var_name = 'DIFF_' + var_name + '_A{}'.format(allele_num2)
        return var_name

    def _get_var_raw(self, iso_num: int, allele_num: int):
        """Get raw variable name."""
        return self._get_var_name('raw', iso_num=iso_num, allele_num=allele_num)

    def get_var_expression(self, allele_num: int, iso_num: Optional[int] = None):
        """Get expression variable name."""
        return self._get_var_name('expression', iso_num=iso_num, allele_num=allele_num)

    def get_var_diff_expression(self, allele_num1: int, allele_num2: int, iso_num: Optional[int] = None):
        """Get expression difference variable name."""
        return self._get_var_name('expression', iso_num=iso_num, allele_num=allele_num1, allele_num2=allele_num2)

    @staticmethod
    def delta_iterator(ploidy: int):
        """Iterator for delta allele number pairs."""
        for m in range(ploidy - 1):
            for n in range(m + 1, ploidy):
                yield m, n


class TaskResult(Task, TaskResultMeta):
    """Task result.

    Attributes:
        trace: Inference results.
    """

    def __init__(self, task: Task, trace: pd.DataFrame):
        super().__init__(task.id, task.segment, task.phased, task.ploidy, task.snps)
        self.trace = self.construct_full_trace(trace)

    def construct_full_trace(self, trace: pd.DataFrame) -> pd.DataFrame:
        """Construct full trace."""
        cols = trace.columns.tolist()
        trace = pd.DataFrame({col: trace[col].tolist() for col in cols}, index=trace.index)[cols]
        assert len(trace.columns) == self.isoforms_count * self.ploidy

        for n in range(self.ploidy):
            cols = [self._get_var_raw(iso_num=i, allele_num=n) for i in range(self.isoforms_count)]
            var = self.get_var_expression(allele_num=n)
            trace[var] = trace[cols].sum(axis=1)

        for i in range(self.isoforms_count):
            for n in range(self.ploidy):
                col = self._get_var_raw(allele_num=n, iso_num=i)
                all_cols = [self._get_var_raw(iso_num=i, allele_num=n) for n in range(self.ploidy)]
                total = trace[all_cols].sum(axis=1)
                var = self.get_var_expression(allele_num=n, iso_num=i)
                trace[var] = trace[col] / total

        for m, n in self.delta_iterator(self.ploidy):
            col1 = self.get_var_expression(allele_num=m)
            col2 = self.get_var_expression(allele_num=n)
            var = self.get_var_diff_expression(allele_num1=m, allele_num2=n)
            trace[var] = np.abs(trace[col1] - trace[col2])
            for i in range(self.isoforms_count):
                col1 = self.get_var_expression(allele_num=m, iso_num=i)
                col2 = self.get_var_expression(allele_num=n, iso_num=i)
                var = self.get_var_diff_expression(allele_num1=m, allele_num2=n, iso_num=i)
                trace[var] = np.abs(trace[col1] - trace[col2])

        return trace

    def stats(self) -> pd.DataFrame:
        """Stats task result."""
        import pymc3 as pm

        idxes = [col for col in self.trace.columns if not col.startswith('RAW')]

        means = []
        sds = []
        hpd_2_5s = []
        hpd_97_5s = []

        cols = ['mean', 'sd', 'hpd_2.5', 'hpd_97.5']
        data = [means, sds, hpd_2_5s, hpd_97_5s]

        for var_name in idxes:
            d_i = self.trace[var_name]
            means.append(np.mean(d_i))
            sds.append(np.std(d_i))
            hpd_2_5, hpd_97_5 = pm.hpd(d_i)
            hpd_2_5s.append(hpd_2_5)
            hpd_97_5s.append(hpd_97_5)

        df = pd.DataFrame({col: data[i] for i, col in enumerate(cols)}, index=idxes, columns=cols)
        return df
