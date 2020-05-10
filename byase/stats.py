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

import os
from typing import List

import pandas as pd

from .message import MessageCenter
from .inference import InferenceTool
from .annotation import AnnotationDB
from .result import ResultDB, ResultRecord
from .task.result import TaskResult, TaskResultMeta


_STATS_DIR_NAME = 'stats'
_TRACE_DIR_NAME = 'trace'
_HEAD_MEAN = 'Mean'
_HEAD_HPD_WIDTH = '95% HPD Width'


class StatsError(Exception):
    """Stats error."""
    def __init__(self, msg):
        super().__init__(msg)


class StatsTool(TaskResultMeta):
    """Stats tool.

    Attributes:
        anno_path: The path of annotation.
        result_path: The path of result.
        out_dir: The path of output directory.

        mc: Message center.
    """

    def __init__(self, anno_path: str, result_path: str, out_dir: str,  mc: MessageCenter):
        self.anno_path = anno_path
        self.result_path = result_path
        self.out_dir = out_dir
        self.mc = mc

    @property
    def _ase_gene_level_stats_path(self):
        """Gene level ASE stats file path."""
        return os.path.join(self.out_dir, 'ASE_geneLevel.csv')

    @property
    def _ase_isoform_level_stats_path(self):
        """Isoform level ASE stats file path."""
        return os.path.join(self.out_dir, 'ASE_isoformLevel.csv')

    @staticmethod
    def _store(df: pd.DataFrame, file_path: str):
        """Store data frame."""
        df.to_csv(file_path, index=False)

    def _stats(self):
        """Generate stats."""
        gene_cols = ['Task ID', 'Gene ID', 'Gene Name', 'Location', 'Isoform Count', 'SNP Count']
        iso_cols = ['Task ID', 'Isoform ID', 'Gene ID', 'Gene Name', 'Isoform Number', 'Isoform Name', 'Location',
                    'SNP Count']

        gene_d = [[] for _ in gene_cols]
        iso_d = [[] for _ in iso_cols]

        def _add_mean(_data: list, _record: ResultRecord, _var_name: str):
            _var_stats = _record.trace_stats.loc[_var_name]
            _data.append(_var_stats['mean'])

        def _add_mean_and_hpd_width(_data: list, _record: ResultRecord, _var_name: str):
            _var_stats = _record.trace_stats.loc[_var_name]
            _data.append(_var_stats['mean'])
            _data.append(_var_stats['hpd_97.5'] - _var_stats['hpd_2.5'])

        def _add_row(_d: List[list], _cols: list, _row: list):
            assert len(_d) == len(_cols) == len(_row)
            for _i, _item in enumerate(_row):
                _d[_i].append(_item)

        n = -1
        ploidy = None
        with AnnotationDB(self.anno_path) as anno_db, ResultDB(self.result_path) as result_db:
            for task in anno_db.tasks_iterator():
                record = result_db.get_record(task.id)
                if not record.success:
                    continue

                if ploidy is None:
                    ploidy = task.ploidy
                    for i in range(ploidy):
                        for cols, d in [(gene_cols, gene_d),
                                        (iso_cols, iso_d)]:
                            cols.append('Allele {} Expression Mean'.format(i+1))
                            d.append([])
                    for i, j in self.delta_iterator(ploidy):
                        for cols, d in [(gene_cols, gene_d),
                                        (iso_cols, iso_d)]:
                            for measure in [_HEAD_MEAN, _HEAD_HPD_WIDTH]:
                                cols.append('Difference (Allele {} & {}) {}'.format(i+1, j+1, measure))
                                d.append([])
                else:
                    assert ploidy == task.ploidy

                assert ploidy == task.ploidy

                n += 1
                if n != 0 and n % 100 == 0:
                    self.mc.handle_progress('Stats: {} tasks processed...'.format(n))

                seg = task.segment
                location = '{}:{}-{}'.format(seg.iv.chrom, seg.iv.start, seg.iv.end)
                gene_row = [task.id, seg.id, seg.gene_name, location, seg.isoforms_count, len(task.snps)]
                for i in range(ploidy):
                    var = self.get_var_expression(allele_num=i)
                    _add_mean(gene_row, record, var)
                for i, j in self.delta_iterator(ploidy):
                    var = self.get_var_diff_expression(allele_num1=i, allele_num2=j)
                    _add_mean_and_hpd_width(gene_row, record, var)
                _add_row(gene_d, gene_cols, gene_row)

                for iso_num, iso in enumerate(seg.isoforms):
                    location = '{}:{}-{}'.format(iso.iv.chrom, iso.iv.start, iso.iv.end)
                    iso_row = [task.id, iso.id, seg.id, seg.gene_name, iso_num + 1, iso.name,
                               location, len(task.isoform_snps(iso_num))]
                    for i in range(ploidy):
                        var = self.get_var_expression(allele_num=i, iso_num=iso_num)
                        _add_mean(iso_row, record, var)
                    for i, j in self.delta_iterator(ploidy):
                        var = self.get_var_diff_expression(allele_num1=i, allele_num2=j, iso_num=iso_num)
                        _add_mean_and_hpd_width(iso_row, record, var)
                    _add_row(iso_d, iso_cols, iso_row)

        def _get_df(_d: List[list], _cols: list):
            assert len(_d) == len(_cols)
            return pd.DataFrame({_col: _d[_i] for _i, _col in enumerate(_cols)})[_cols]

        gene_df = _get_df(gene_d, gene_cols)
        iso_df = _get_df(iso_d, iso_cols)

        self._store(gene_df, self._ase_gene_level_stats_path)
        self._store(iso_df, self._ase_isoform_level_stats_path)

    def stats(self):
        """Generate stats."""
        already_exists = True
        for file_path in [self._ase_gene_level_stats_path, self._ase_isoform_level_stats_path]:
            if not os.path.exists(file_path):
                already_exists = False

        if not already_exists:
            self._stats()

        return {'gene-level path': self._ase_gene_level_stats_path,
                'isoform-level path': self._ase_isoform_level_stats_path}


def stats(args):
    """Generate stats."""
    result_dir = args['result_dir']

    mc = args['mc']  # type: MessageCenter

    mc.log_debug('result_dir: {}'.format(result_dir))

    inference_tool = InferenceTool(result_dir, n_process=1, param=None, mc=mc)

    stats_dir = os.path.join(result_dir, _STATS_DIR_NAME)
    if not os.path.exists(stats_dir):
        os.mkdir(stats_dir)

    stats_tool = StatsTool(anno_path=inference_tool.anno_path, result_path=inference_tool.result_path,
                           out_dir=stats_dir, mc=mc)
    return stats_tool.stats()


def trace(args):
    """Get MCMC traces."""
    result_dir = args['result_dir']
    task_id = args['task_id']

    mc = args['mc']  # type: MessageCenter

    mc.log_debug('result_dir: {}'.format(result_dir))
    mc.log_debug('task_id: {}'.format(task_id))

    trace_dir = os.path.join(result_dir, _TRACE_DIR_NAME)
    if not os.path.exists(trace_dir):
        os.mkdir(trace_dir)

    out_path = os.path.join(trace_dir, '{}.csv'.format(task_id))

    inference_tool = InferenceTool(result_dir, n_process=1, param=None, mc=mc)
    anno_path = inference_tool.anno_path
    result_path = inference_tool.result_path

    with AnnotationDB(anno_path) as anno_db, ResultDB(result_path) as result_db:
        task = anno_db.get_task(task_id)
        record = result_db.get_record(task_id)
        result = TaskResult(task, record.trace)
        result.trace.to_csv(out_path, index=False)
