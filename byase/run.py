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

import argparse
import warnings
import traceback

from .message import MessageCenter


def _add_gen_annotation_parser(subparsers):
    parser = subparsers.add_parser('gen-task', help='Generate tasks.')

    parser.add_argument('-g', '--gff', required=True, help='The path of the GFF annotation file.')
    parser.add_argument('-G', '--gene-feature', default='gene',
                        help='The feature type name of genes, '
                             'in the 3rd column of the GFF file. (DEFAULT: gene)')
    parser.add_argument('-I', '--isoform-feature', default='transcript',
                        help='The feature type name of isoforms, '
                             'in the 3rd column of the GFF file. (DEFAULT: transcript)')
    parser.add_argument('-A', '--gene-name-attr', default='ID',
                        help='The attribute that describes gene names, '
                             'extracted from the 9th column of the GFF file. (DEFAULT: ID)')
    parser.add_argument('-B', '--isoform-name-attr', default='ID',
                        help='The attribute that describes isoform names, '
                             'extracted from the 9th column of the GFF file. (DEFAULT: ID)')

    parser.add_argument('-v', '--vcf', required=True, help='The path of the VCF file.')
    parser.add_argument('-s', '--sample', required=True, help='The name of the target sample in the VCF file.')
    parser.add_argument('-p', '--ploidy', type=int, default=2, help='The ploidy of the target sample. (DEFAULT: 2)')
    parser.add_argument('-C', '--add_chrom_prefix', action='store_true',
                        help='Add "chr" prefix to chromosome names when extracting SNPs from the VCF file.')

    parser.add_argument('-o', '--out-dir', required=True, help='The path of the output directory.')


def _add_extract_region_parser(subparsers):
    parser = subparsers.add_parser('extract-region', help='Extract long continuous regions.')

    parser.add_argument('-g', '--gff', required=True, help='The path of the GFF3 annotation file.')
    parser.add_argument('-m', '--min-region-len', type=int, default=1000,
                        help='Minimum length of regions. (DEFAULT: 1000)')

    parser.add_argument('-o', '--out', required=True, help='The path of the output file.')


def _add_stat_insert_size_parser(subparsers):
    parser = subparsers.add_parser('stat-insert-size', help='Stat insert size of paired-end BAM files.')

    parser.add_argument('-r', '--region', required=True,
                        help='The path of the file containing long continuous regions.')
    parser.add_argument('-b', '--bams', nargs='+', required=True, help='The paths of the BAMs.')

    parser.add_argument('-o', '--out-dir', required=True, help='The path of output directory.')


def _add_inference_parser(subparsers):
    parser = subparsers.add_parser('inference', help='Perform Bayesian inference on tasks.')

    parser.add_argument('-t', '--task', required=True, help='The path of task directory.')

    parser.add_argument('-b', '--bam', nargs='+', required=True, help='The paths of BAM files.')
    parser.add_argument('-L', '--read-len', type=int, required=True, help='Read length.')

    parser.add_argument('-P', '--pe', action='store_true', help='Specify the BAM files are paired-end.')
    parser.add_argument('-M', '--insert-size-mean', type=float, help='The mean of insert size.')
    parser.add_argument('-S', '--insert-size-std', type=float, help='The standard deviation of insert size.')

    parser.add_argument('-o', '--out-dir', required=True, help='The path of output directory.')

    parser.add_argument('-n', '--process', type=int, default=1, help='The number of processes used for '
                                                                     'parallel computing. (DEFAULT: 1)')
    parser.add_argument('-c', '--count', type=int, help='The count of tasks to be inferred. '
                                                        'If not specified, inference will be done on all tasks.')


def _add_inference_resume_parser(subparsers):
    parser = subparsers.add_parser('resume', help='Resume Bayesian inference.')

    parser.add_argument('-o', '--out-dir', required=True, help='The path of previous output directory.')

    parser.add_argument('-n', '--process', type=int, default=1, help='The number of processes used for '
                                                                     'parallel computing. (DEFAULT: 1)')
    parser.add_argument('-c', '--count', type=int, help='The count of tasks to be inferred. '
                                                        'If not specified, inference will be done on all tasks.')


def _add_stats_parser(subparsers):
    parser = subparsers.add_parser('stats', help='Stats results.')

    parser.add_argument('-d', '--result-dir', required=True, help='The paths of result directory.')


def _add_plot_parser(subparsers):
    parser = subparsers.add_parser('plot', help='Plot.')

    parser.add_argument('-d', '--result-dir', required=True, help='The path of result directory.')
    parser.add_argument('-i', '--task-id', required=True, help='The ID of the task to plot.')


def _parse_args():
    parser = argparse.ArgumentParser(prog='byase')
    subparsers = parser.add_subparsers(dest='sub_command')

    _add_gen_annotation_parser(subparsers)
    _add_extract_region_parser(subparsers)
    _add_stat_insert_size_parser(subparsers)
    _add_inference_parser(subparsers)
    _add_inference_resume_parser(subparsers)
    _add_stats_parser(subparsers)
    _add_plot_parser(subparsers)

    args = parser.parse_args()

    return args


def main():
    # Hide HTSeq warnings.
    warnings.filterwarnings('ignore')

    args = vars(_parse_args())
    mc = MessageCenter()
    args['mc'] = mc

    sub_command = args['sub_command']

    try:
        if sub_command == 'gen-task':
            from .annotation import generate_annotation
            generate_annotation(args)

        elif sub_command == 'extract-region':
            from .pe_utils import extract_region
            extract_region(args)

        elif sub_command == 'stat-insert-size':
            from .pe_utils import stat_insert_size
            stat_insert_size(args)

        elif sub_command == 'inference':
            from .inference import inference
            inference(args)

        elif sub_command == 'resume':
            from .inference import inference_resume
            inference_resume(args)

        elif sub_command == 'stats':
            from .stats import stats
            stats(args)

        elif sub_command == 'plot':
            from .plot import plot_task
            plot_task(args)

    except Exception as e:
        mc.log_error(e)
        traceback.print_exc()
