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

import numpy as np
import scipy.stats
import HTSeq

from .message import MessageCenter


_INSERT_SIZE_NAME = 'insert_size.txt'
_INSERT_SIZE_STAT_NAME = 'insert_size_stat.txt'
_INSERT_SIZE_STAT_PLOT_NAME = 'insert_size_stat.png'


class PEUtilError(Exception):
    """Paired-end utils error."""
    def __init__(self, msg):
        super().__init__(msg)


class PEUtilPathError(PEUtilError):
    """Path error."""
    def __init__(self, path, msg):
        super().__init__('[PATH: {}] {}'.format(path, msg))


class PEUtilParseError(PEUtilError):
    """Parsing error."""
    def __init__(self, path, msg):
        super().__init__('[PARSE: {}] {}'.format(path, msg))


def _extract_long_continuous_regions(gff_path: str, min_region_len: int, out_path: str, mc: MessageCenter):
    """Extract long continuous regions with length of at least min_region_len."""
    mc.log_debug('gff_path: {}'.format(gff_path))
    mc.log_debug('min_region_len: {}'.format(min_region_len))
    mc.log_debug('out_path: {}'.format(out_path))

    mc.handle_progress('Calculating long continuous regions...')

    region = HTSeq.GenomicArray('auto', stranded=False, typecode='i')

    if not os.path.exists(gff_path):
        raise PEUtilPathError(gff_path, 'File not exists.')

    gff = HTSeq.GFF_Reader(gff_path)
    n = -1
    for ft in gff:
        n += 1
        if n != 0 and n % 100000 == 0:
            mc.handle_progress('{} lines read from GFF file...'.format(n))
        if ft.type == 'exon':
            region[ft.iv] += 1

    with open(out_path, 'w') as o:
        for iv, v in region.steps():
            if v != 0:
                region_len = iv.end - iv.start
                if region_len >= min_region_len:
                    o.write('{0}\t{1}\t{2}\t{3}\n'.format(iv.chrom, iv.start, iv.end, region_len))


def _calc_insert_size(regions_path: str, bam_paths: List[str], out_path: str, mc: MessageCenter):
    """Calculate insert sizes of each read pairs mapped to long continuous regions."""
    mc.log_debug('regions_path: {}'.format(regions_path))
    mc.log_debug('bam_paths: {}'.format(', '.join(bam_paths)))
    mc.log_debug('out_path: {}'.format(out_path))

    mc.handle_progress('Collecting insert sizes...')

    if not os.path.exists(regions_path):
        raise PEUtilPathError(regions_path, 'File not exists.')

    bam_readers = []
    for bam_path in bam_paths:
        if not os.path.exists(bam_path):
            raise PEUtilPathError(bam_path, 'File not exists.')
        bam_reader = HTSeq.BAM_Reader(bam_path)
        bam_readers.append(bam_reader)

    with open(regions_path) as f, open(out_path, 'w') as o:
        n = -1
        for row in f:
            row = row.strip('\n')
            n += 1
            if n != 0 and n % 1000 == 0:
                mc.handle_progress('{} regions processed...'.format(n))
            cells = row.split('\t')
            try:
                chrom, start, end = cells[:3]
                start, end = int(start), int(end)
            except ValueError:
                raise PEUtilParseError(regions_path, 'Incorrect file format.')
            insert_sizes = []
            for bam_reader in bam_readers:
                alns = bam_reader.fetch(chrom, start, end)
                for aln1, aln2 in HTSeq.pair_SAM_alignments_with_buffer(alns):
                    aln1 = aln1  # type: HTSeq.SAM_Alignment
                    aln2 = aln2  # type: HTSeq.SAM_Alignment
                    if (aln1 is None) or (aln2 is None):
                        continue

                    should_skip = False
                    for aln in [aln1, aln2]:
                        assert aln.aligned
                        assert aln.iv.start < aln.iv.end
                        if aln.not_primary_alignment or aln.pcr_or_optical_duplicate or aln.failed_platform_qc:
                            should_skip = True
                            break
                        if (aln.iv.start < start) or (end < aln.iv.end):
                            should_skip = True
                            break
                    if should_skip:
                        continue

                    insert_sizes.append(str(np.abs(aln1.inferred_insert_size)))
            o.write('{0}:{1}-{2}({3})\t{4}\n'.format(chrom, start, end, end - start, ';'.join(insert_sizes)))


def _stat_insert_size(insert_size_path: str, out_stat_path: str, out_plot_path: str, mc: MessageCenter):
    """Stat insert sizes."""
    mc.log_debug('insert_size_path: {}'.format(insert_size_path))
    mc.log_debug('out_stat_path: {}'.format(out_stat_path))
    mc.log_debug('out_plot_path: {}'.format(out_plot_path))

    if not os.path.exists(insert_size_path):
        raise PEUtilPathError(insert_size_path, 'File not exists.')

    mc.handle_progress('Calculating insert size stat...')

    insert_sizes = []
    with open(insert_size_path) as f:
        for row in f:
            lens = row.strip('\n').split('\t')[-1]
            if lens == '':
                continue
            lens = lens.split(';')
            for i in lens:
                insert_sizes.append(int(i))
    insert_sizes = np.array(insert_sizes)

    # Remove outliers.
    m = np.mean(insert_sizes)
    insert_sizes = insert_sizes[insert_sizes < 2 * m]

    m = np.mean(insert_sizes)
    s = np.std(insert_sizes)
    d = s / np.sqrt(m)

    with open(out_stat_path, 'w') as o:
        o.write('mean={}\nstd={}\ndispersion={}\npairs_count={}\n'.format(m, s, d, len(insert_sizes)))

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)

    ax.hist(insert_sizes, bins=np.arange(min(insert_sizes), max(insert_sizes), 1), density=True,
            histtype='stepfilled', alpha=0.5)

    x = np.arange(min(insert_sizes), max(insert_sizes), 0.1)
    y = scipy.stats.norm.pdf(x, loc=m, scale=s)
    ax.plot(x, y, linestyle='--')

    ax.set_title(r'Distribution of insert size: $\mu={:.2f}$, $\sigma={:.2f}$'.format(m, s), fontsize=10)
    ax.set_xlabel('Insert size')
    ax.set_ylabel('Probability density')

    fig.savefig(out_plot_path)


def extract_region(args):
    """Extract long continuous regions."""
    gff_path = args['gff']
    min_region_len = args['min_region_len']
    out = args['out']
    mc = args['mc']
    _extract_long_continuous_regions(gff_path, min_region_len, out, mc)


def stat_insert_size(args):
    regions_path = args['region']
    bam_paths = args['bams']
    out_dir = args['out_dir']

    mc = args['mc']

    if not os.path.isdir(out_dir):
        raise PEUtilPathError(out_dir, 'Invalid directory.')

    insert_size_path = '{}/{}'.format(out_dir, _INSERT_SIZE_NAME)
    insert_size_stat_path = '{}/{}'.format(out_dir, _INSERT_SIZE_STAT_NAME)
    insert_size_stat_plot_path = '{}/{}'.format(out_dir, _INSERT_SIZE_STAT_PLOT_NAME)

    _calc_insert_size(regions_path, bam_paths, insert_size_path, mc)
    _stat_insert_size(insert_size_path, insert_size_stat_path, insert_size_stat_plot_path, mc)
