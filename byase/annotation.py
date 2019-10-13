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
from typing import List, Optional

import HTSeq

from .message import MessageCenter
from .component import Segment, Isoform, SNP
from .task import Task
from .db import DB, Field


class AnnotationError(Exception):
    """Annotation error."""
    def __init__(self, msg):
        super().__init__(msg)


class AnnotationPathError(AnnotationError):
    """Path error."""
    def __init__(self, path, msg):
        super().__init__('[PATH: {}] {}'.format(path, msg))


class AnnotationParseError(AnnotationError):
    """Parsing error."""
    def __init__(self, path, msg):
        super().__init__('[PARSE: {}] {}'.format(path, msg))


class AnnotationDBError(AnnotationError):
    """Database error."""
    def __init__(self, path, msg):
        super().__init__('[DATABASE: {}] {}'.format(path, msg))


def _read_gff(gff_path: str, segment_feature_name: str, isoform_feature_name: str,
              gene_name_attr: str, isoform_name_attr: str, mc: MessageCenter) -> List[Segment]:
    """Read all segments from GFF file.

    Args:
        gff_path: The path of GFF file.
        segment_feature_name: Segment feature name in GFF hierarchy.
        isoform_feature_name: Isoform feature name in GFF hierarchy.
        gene_name_attr: Name of attribute which describes gene name.
        isoform_name_attr: Name of attribute which describes isoform name.
        mc: Message center.
    """

    mc.log_debug('gff_path: {}'.format(gff_path))
    mc.log_debug('segment_feature_name: {}'.format(segment_feature_name))
    mc.log_debug('isoform_feature_name: {}'.format(isoform_feature_name))
    mc.log_debug('gene_name_attr: {}'.format(gene_name_attr))
    mc.log_debug('isoform_name_attr: {}'.format(isoform_name_attr))

    mc.handle_progress('Reading GFF file...')

    if not os.path.exists(gff_path):
        raise AnnotationPathError(gff_path, 'File not exists.')

    segments = []
    segments_dict = {}
    isoforms_dict = {}

    gff = HTSeq.GFF_Reader(gff_path)
    n = -1
    for ft in gff:
        n += 1
        if n != 0 and n % 100000 == 0:
            mc.handle_progress('{} lines read from GFF file...'.format(n))

        if ft.type == segment_feature_name:
            if gene_name_attr not in ft.attr:
                raise AnnotationParseError(gff_path, '"{}" not in attributes.'.format(gene_name_attr))
            segment = Segment(ft.attr['ID'], ft.attr[gene_name_attr], ft.iv)
            if segment.id in segments_dict:
                raise AnnotationParseError(gff_path, 'Segment ID is not unique: {}'.format(segment.id))
            segments_dict[segment.id] = segment
            segments.append(segment)
        elif ft.type == isoform_feature_name:
            if isoform_name_attr not in ft.attr:
                raise AnnotationParseError(gff_path, '"{}" not in attributes.'.format(isoform_name_attr))
            isoform = Isoform(ft.attr['ID'], ft.attr[isoform_name_attr], ft.iv)
            if isoform.id in isoforms_dict:
                raise AnnotationParseError(gff_path, 'Isoform ID is not unique: {}'.format(isoform.id))
            isoforms_dict[isoform.id] = isoform
            parent_id = ft.attr['Parent']
            if parent_id not in segments_dict:
                raise AnnotationParseError(gff_path, 'Cannot find the parent of isoform "{}", the segment feature '
                                                     'type name may be incorrect.'.format(isoform.id))
            parent = segments_dict[parent_id]
            assert parent.iv.chrom == isoform.iv.chrom
            assert parent.iv.strand == isoform.iv.strand
            parent.isoforms.append(isoform)
        elif ft.type == 'exon':
            parent_id = ft.attr['Parent']
            if parent_id not in isoforms_dict:
                raise AnnotationParseError(gff_path, 'Cannot find the parent of exon "{}", the isoform feature '
                                                     'type name may be incorrect.'.format(ft.attr['ID']))
            parent = isoforms_dict[parent_id]
            assert parent.iv.chrom == ft.iv.chrom
            assert parent.iv.strand == ft.iv.strand
            parent.exons.append(ft.iv)

    return segments


class HeteroSNPReader:
    """Heterozygous SNP reader.

    Attributes:
        vcf_path: The path of VCF file.
        sample: The name of sample to be extracted.
        ploidy: The ploidy of the sample
        add_chrom_prefix: If add "chr" prefix to chromosome.
        mc: Message center.
    """

    def __init__(self, vcf_path: str, sample: str, ploidy: int, add_chrom_prefix: bool, mc: MessageCenter):
        self.vcf_path = vcf_path
        self.sample = sample
        self.ploidy = ploidy
        self.add_chrom_prefix = add_chrom_prefix
        self.mc = mc

        if not os.path.exists(self.vcf_path):
            raise AnnotationPathError(self.vcf_path, 'File not exists.')

    def __iter__(self):
        self.mc.log_debug('vcf_path: {}'.format(self.vcf_path))
        self.mc.log_debug('sample: {}'.format(self.sample))
        self.mc.log_debug('ploidy: {}'.format(self.ploidy))
        self.mc.log_debug('add_chrom_prefix: {}'.format(self.add_chrom_prefix))

        vcf = HTSeq.VCF_Reader(self.vcf_path)
        vcf.parse_meta()

        self.mc.handle_progress('Reading VCF file...')

        n = -1
        for vc in vcf:
            n += 1
            if n != 0 and n % 500000 == 0:
                self.mc.handle_progress('{} lines read from VCF file...'.format(n))

            if self.sample not in vc.samples:
                raise AnnotationParseError(self.vcf_path, 'Sample "{}" not in VCF file.'.format(self.sample))

            gt = vc.samples[self.sample]['GT']
            if '.' in gt:
                continue
            if '/' in gt:
                phased = False
                if '|' in gt:
                    gt = gt.replace('|', '/')
                sep = '/'
            else:
                assert '|' in gt
                phased = True
                sep = '|'

            gt = gt.split(sep)
            if len(gt) != self.ploidy:
                raise AnnotationParseError(self.vcf_path,
                                           'The ploidy({}) may be inconsistent with the '
                                           'sample "{}"({}).'.format(self.ploidy, self.sample, len(gt)))

            ref_alt = [vc.ref] + vc.alt
            alleles = [ref_alt[int(g)] for g in gt]

            for allele in alleles:
                if len(allele) != 1:
                    continue
            if len(set(alleles)) < 2:
                continue

            chrom = vc.pos.chrom
            if self.add_chrom_prefix:
                chrom = 'chr{}'.format(chrom)
            pos = vc.pos.pos - 1

            snp = SNP(chrom, pos, alleles, phased)
            assert self.ploidy == snp.ploidy
            yield snp


_SEGMENT_DB_FILENAME = 'segment.db'
_ISOFORM_DB_FILENAME = 'isoform.db'
_SNP_DB_FILENAME = 'snp.db'
_TASK_DB_FILENAME = 'task.db'


class AnnotationDB:
    """Annotation Database."""

    _segment_scheme = [Field('id', 'S'), Field('gene_name', 'S'), Field('chrom', 'S'),
                       Field('start', 'I'), Field('end', 'I'), Field('strand', 'S'),
                       Field('isoforms_count', 'I'), Field('isoforms', 'S'), Field('isoforms_db_offsets', 'S')]
    _isoform_schema = [Field('id', 'S'), Field('name', 'S'), Field('start', 'I'), Field('end', 'I'),
                       Field('exons_count', 'I'), Field('exons_starts', 'S'), Field('exons_ends', 'S')]
    _snp_schema = [Field('id', 'S'), Field('chrom', 'S'), Field('pos', 'I'), Field('phased', 'I'),
                   Field('ploidy', 'I'), Field('alleles', 'S')]
    _task_schema = [Field('id', 'S'), Field('segment', 'S'), Field('phased', 'I'), Field('ploidy', 'I'),
                    Field('snps_count', 'I'), Field('snps', 'S')]

    def __init__(self, dir_path: str, initialize: bool = False, mc: Optional[MessageCenter] = None):
        self._mc = mc

        if not os.path.isdir(dir_path):
            raise AnnotationPathError(dir_path, 'Invalid directory.')

        read_only = not initialize

        self._segment_db = DB(os.path.join(dir_path, _SEGMENT_DB_FILENAME),
                              AnnotationDB._segment_scheme if initialize else None, read_only=read_only)
        self._isoform_db = DB(os.path.join(dir_path, _ISOFORM_DB_FILENAME),
                              AnnotationDB._isoform_schema if initialize else None, read_only=read_only)
        self._snp_db = DB(os.path.join(dir_path, _SNP_DB_FILENAME),
                          AnnotationDB._snp_schema if initialize else None, read_only=read_only)
        self._task_db = DB(os.path.join(dir_path, _TASK_DB_FILENAME),
                           AnnotationDB._task_schema if initialize else None, read_only=read_only)

        if not initialize:
            assert self._segment_db.read_only
            assert self._isoform_db.read_only
            assert self._snp_db.read_only
            assert self._task_db.read_only

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._segment_db.close()
        self._isoform_db.close()
        self._snp_db.close()
        self._task_db.close()

    def get_all_segment_ids(self):
        """Get all segment IDs."""
        return self._segment_db.get_all_item_ids()

    def get_segment_basic_info(self, segment_id):
        """Get segment basic information."""
        item = self._segment_db.get_item(segment_id)
        _, gene_name, chrom, start, end, _, isoforms_count, _, _ = item
        return {'gene_name': gene_name, 'location': '{}:{}-{}'.format(chrom, start, end),
                'isoforms_count': isoforms_count}

    def get_segment(self, segment_id):
        """Get segment by ID."""
        item = self._segment_db.get_item(segment_id)
        segment_id, gene_name, chrom, start, end, strand, isoforms_count, isoform_ids, isoform_offsets = item
        iv = HTSeq.GenomicInterval(chrom, start, end, strand)
        segment = Segment(segment_id, gene_name, iv)

        isoform_ids = isoform_ids.split(';')
        isoform_offsets = [int(offset) for offset in isoform_offsets.split(';')]
        assert isoforms_count == len(isoform_ids)
        assert isoforms_count == len(isoform_offsets)
        for n, isoform_offset in enumerate(isoform_offsets):
            item = self._isoform_db.get_item_by_offset(isoform_offset)
            isoform_id, isoform_name, start, end, exons_count, exons_starts, exons_ends = item
            assert isoform_id == isoform_ids[n]
            iv = HTSeq.GenomicInterval(chrom, start, end, strand)
            isoform = Isoform(isoform_id, isoform_name, iv)
            segment.isoforms.append(isoform)

            exons_starts = [int(s) for s in exons_starts.split(';')]
            exons_ends = [int(e) for e in exons_ends.split(';')]
            assert exons_count == len(exons_starts) == len(exons_ends)
            for i in range(exons_count):
                exon_iv = HTSeq.GenomicInterval(chrom, exons_starts[i], exons_ends[i], strand)
                isoform.exons.append(exon_iv)
            assert len(isoform.exons) == exons_count
        assert len(segment.isoforms) == isoforms_count
        return segment

    def get_snp(self, snp_id):
        """Get SNP by ID."""
        snp_id, chrom, pos, phased_code, ploidy, alleles = self._snp_db.get_item(snp_id)
        alleles = alleles.split(';')
        assert len(alleles) == ploidy
        if phased_code == 1:
            phased = True
        else:
            assert phased_code == 0
            phased = False
        snp = SNP(chrom, pos, alleles, phased)
        return snp

    def get_all_task_ids(self):
        """Get all task IDs."""
        return self._task_db.get_all_item_ids()

    def _construct_task(self, item):
        """Construct task from database item."""
        task_id, segment_id, phased_code, ploidy, snps_count, snp_ids = item
        if phased_code == 1:
            phased = True
        else:
            assert phased_code == 0
            phased = False
        snp_ids = snp_ids.split(';')
        assert len(snp_ids) == snps_count
        segment = self.get_segment(segment_id)
        snps = [self.get_snp(snp_id) for snp_id in snp_ids]
        task = Task(task_id, segment, phased, ploidy, snps)
        return task

    def tasks_iterator(self):
        """Tasks iterator."""
        for item in self._task_db.item_iterator():
            yield self._construct_task(item)

    def tasks_db_fields_iterator(self):
        """Tasks database fields iterator."""
        for item in self._task_db.item_iterator():
            yield {self._task_schema[n].name: item[n] for n in range(len(self._task_schema))}

    def get_task(self, task_id):
        """Get task by ID."""
        item = self._task_db.get_item(task_id)
        return self._construct_task(item)

    def store_segments(self, segments: List[Segment]):
        """Store segments into the database."""
        for segment in segments:
            for isoform in segment.isoforms:
                self._isoform_db.store_item(
                    (isoform.id, isoform.name, isoform.iv.start, isoform.iv.end, len(isoform.exons),
                     ';'.join([str(e.start) for e in isoform.exons]),
                     ';'.join([str(e.end) for e in isoform.exons])))
        for segment in segments:
            isoform_offsets = []
            for isoform in segment.isoforms:
                offset = self._isoform_db.find_item_offset(isoform.id)
                assert offset is not None
                isoform_offsets.append(str(offset))
            self._segment_db.store_item(
                (segment.id, segment.gene_name, segment.iv.chrom, segment.iv.start, segment.iv.end, segment.iv.strand,
                 len(segment.isoforms), ';'.join([i.id for i in segment.isoforms]), ';'.join(isoform_offsets)))

    def store_snps(self, snps):
        """Store SNPs into the database."""
        for snp in snps:
            snp = snp  # type: SNP
            self._snp_db.store_item(
                (snp.id, snp.iv.chrom, snp.iv.start, 1 if snp.phased else 0, len(snp.alleles), ';'.join(snp.alleles)))

    def create_tasks(self, ploidy: int, snps):
        """Create tasks."""
        vc_sites = HTSeq.GenomicArray('auto', stranded=False, typecode='O')

        for snp in snps:
            snp = snp  # type: SNP
            vc_sites[snp.iv] = snp

        self._mc.handle_progress('Creating tasks...')
        selected_snp_ids = set()
        selected_snps = []
        n = -1
        segment_ids = self.get_all_segment_ids()
        for segment_id in segment_ids:
            segment = self.get_segment(segment_id)
            phased_snp_ids, unphased_snp_ids = [], []
            for isoform in segment.isoforms:
                for exon in isoform.exons:
                    for vc_iv, snp in vc_sites[exon].steps():
                        if snp is not None:
                            if snp.id not in selected_snp_ids:
                                selected_snp_ids.add(snp.id)
                                selected_snps.append(snp)
                            if snp.phased and snp.id not in phased_snp_ids:
                                phased_snp_ids.append(snp.id)
                            if not snp.phased and snp.id not in unphased_snp_ids:
                                unphased_snp_ids.append(snp.id)
            if len(phased_snp_ids) > 0:
                n = self._create_task(n, segment.id, True, ploidy, phased_snp_ids)
            for unphased_snp_id in unphased_snp_ids:
                n = self._create_task(n, segment.id, False, ploidy, [unphased_snp_id])

        self.store_snps(selected_snps)

    def _create_task(self, n: int, segment_id: str, phased: bool, ploidy: int, snp_ids: List[str]) -> int:
        """Create task."""
        assert self._mc is not None

        n += 1
        if n != 0 and n % 10000 == 0:
            self._mc.handle_progress('{} tasks created...'.format(n))
        if phased:
            task_id = '{}-PHASED'.format(segment_id)
        else:
            assert len(snp_ids) == 1
            task_id = '{}-{}'.format(segment_id, snp_ids[0])
        self._task_db.store_item((task_id, segment_id, 1 if phased else 0, ploidy, len(snp_ids), ';'.join(snp_ids)))
        return n


def generate_annotation(args):
    """Generate annotation."""
    gff_path = args['gff']
    segment_feature_name = args['gene_feature']
    isoform_feature_name = args['isoform_feature']
    gene_name_attr = args['gene_name_attr']
    isoform_name_attr = args['isoform_name_attr']

    vcf_path = args['vcf']
    vcf_sample = args['sample']
    ploidy = args['ploidy']
    add_chrom_prefix = args['add_chrom_prefix']

    out_dir = args['out_dir']
    mc = args['mc']  # type: MessageCenter

    if not os.path.isdir(out_dir):
        raise AnnotationPathError(out_dir, 'Invalid directory.')

    with AnnotationDB(out_dir, initialize=True, mc=mc) as anno_db:
        segments = _read_gff(gff_path, segment_feature_name, isoform_feature_name,
                             gene_name_attr, isoform_name_attr, mc)
        anno_db.store_segments(segments)
        del segments
        anno_db.create_tasks(ploidy, HeteroSNPReader(vcf_path, vcf_sample, ploidy, add_chrom_prefix, mc))
