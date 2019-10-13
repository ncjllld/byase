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

from typing import List, Tuple, Set, Optional
from typing import NamedTuple

import numpy as np
import HTSeq

from ..message import MessageCenter
from .task import SNP, Task


_PHRED_QUAL_THRESHOLD = 20


class BAMParam:
    """BAM parameters.

    Attributes:
        paths: The paths of BAM files.
        read_len: Read length.
        paired_end: If the BAMs are paired-end.
        insert_size_mean: Mean of insert-size in paired-end data.
        insert_size_std: Std of insert-size in paired-end data.
    """

    def __init__(self, bam_paths: List[str], read_len: int, paired_end: bool,
                 insert_size_mean: Optional[float], insert_size_std: Optional[float]):
        self.paths = bam_paths
        self.read_len = read_len
        self.paired_end = paired_end
        self.insert_size_mean = insert_size_mean
        self.insert_size_std = insert_size_std

        if self.paired_end:
            assert self.insert_size_mean is not None
            assert self.insert_size_std is not None
        else:
            assert self.insert_size_mean is None
            assert self.insert_size_std is None


class _MatchedRegion(NamedTuple):
    iv: HTSeq.GenomicInterval
    seq: Optional[str]
    qual: Optional[np.ndarray]


class _SkippedRegion(NamedTuple):
    iv: HTSeq.GenomicInterval


class AlignStat:
    """Alignment stats.

    Attributes:
        aln: The original SAM alignment.
        aln2: The other SAM alignment if the BAM file is paired-end.
        matched_ivs: Matched genomic intervals.
        compatible_isoforms: Compatible isoforms.
        compatible_alleles: Compatible alleles.
    """

    def __init__(self, aln: HTSeq.SAM_Alignment, aln2: Optional[HTSeq.SAM_Alignment],
                 matched_ivs: List[HTSeq.GenomicInterval],
                 compatible_isoforms: Set[int], compatible_alleles: Set[int]):
        self.aln = aln
        self.aln2 = aln2
        self.matched_ivs = matched_ivs
        self.compatible_isoforms = compatible_isoforms
        self.compatible_alleles = compatible_alleles

    @property
    def invalid(self):
        """If the alignment is invalid."""
        return (len(self.compatible_isoforms) == 0) or (len(self.compatible_alleles) == 0)


class TaskAlign(Task):
    """Task alignment.

        Attributes:
            region: A GenomicArrayOfSets object containing isoforms and SNPs info.
            bam_param: The parameters of BAM files.
            mc: Message center.
        """

    def __init__(self, task: Task, bam_param: BAMParam, mc: MessageCenter):
        super().__init__(task.id, task.segment, task.phased, task.ploidy, task.snps)
        self.bam_param = bam_param
        self.mc = mc

        # Construct the region.
        self.region = HTSeq.GenomicArrayOfSets('auto', stranded=False)
        for i, isoform in enumerate(self.segment.isoforms):
            for exon in isoform.exons:
                self.region[exon] += i
        for snp in self.snps:
            self.region[snp.iv] += snp

    def _validate_aln(self, aln: HTSeq.SAM_Alignment) -> bool:
        """Validate alignment"""
        if not aln.aligned:
            return False
        if aln.not_primary_alignment or aln.pcr_or_optical_duplicate or aln.failed_platform_qc:
            return False
        if (aln.iv.start < self.segment.iv.start) or (self.segment.iv.end < aln.iv.end):
            return False
        return True

    def _extract_matched_and_skipped_regions(
            self, aln: HTSeq.SAM_Alignment) -> Tuple[List[_MatchedRegion], List[_SkippedRegion]]:
        """Get matched and skipped regions."""
        matched_regions = []  # type: List[_MatchedRegion]
        skipped_regions = []  # type: List[_SkippedRegion]

        last_matched_end = None

        for co in aln.cigar:
            co = co  # type: HTSeq.CigarOperation
            # CIGAR string 'D'(deletion) should be considered as matched to keep consistency.
            if co.type in ['M', '=', 'X', 'I', 'D']:
                co_iv = co.ref_iv  # type: HTSeq.GenomicInterval

                if co.type in ['M', '=', 'X']:
                    matched_seq = aln.read_as_aligned.seq[co.query_from: co.query_to].decode()
                    matched_qual = aln.read_as_aligned.qual[co.query_from: co.query_to]
                    matched_region = _MatchedRegion(iv=co_iv, seq=matched_seq, qual=matched_qual)
                    matched_regions.append(matched_region)

                if (last_matched_end is not None) and (last_matched_end < co_iv.start):
                    skipped_iv = HTSeq.GenomicInterval(self.segment.iv.chrom, last_matched_end, co_iv.start, '.')
                    skipped_region = _SkippedRegion(iv=skipped_iv)
                    skipped_regions.append(skipped_region)

                last_matched_end = co_iv.end
        return matched_regions, skipped_regions

    def _calc_compatible_isoforms(self, matched_regions: List[_MatchedRegion],
                                  skipped_regions: List[_SkippedRegion]) -> Set[int]:
        """Calculate compatible isoforms."""

        # Calculate compatible isoforms from matched regions.
        compatible_set = None

        for matched_region in matched_regions:
            for iv, value in self.region[matched_region.iv].steps():
                isoforms_set = set([v for v in value if isinstance(v, int)])
                if compatible_set is None:
                    compatible_set = isoforms_set
                else:
                    compatible_set &= isoforms_set

        if compatible_set is None:
            return set()

        # Calculate incompatible isoforms from skipped regions.
        incompatible_set = set()

        for skipped_region in skipped_regions:
            for iv, value in self.region[skipped_region.iv].steps():
                isoforms_set = set([v for v in value if isinstance(v, int)])
                incompatible_set |= isoforms_set

        # Incorporate compatible isoforms with incompatible isoforms.
        compatible_set -= incompatible_set
        return compatible_set

    def _calc_compatible_alleles(self, matched_regions: List[_MatchedRegion]) -> Set[int]:
        """Calculate compatible alleles."""
        compatible_set = set([i for i in range(self.ploidy)])

        for matched_region in matched_regions:
            for iv, value in self.region[matched_region.iv].steps():
                snps = [v for v in value if isinstance(v, SNP)]
                for snp in snps:
                    snp = snp  # type: SNP
                    nt_idx = snp.pos - matched_region.iv.start
                    nt = matched_region.seq[nt_idx]

                    # Filter base quality.
                    nt_qual = matched_region.qual[nt_idx]
                    if nt_qual < _PHRED_QUAL_THRESHOLD:
                        continue

                    alleles_set = set([i for (i, snp_a) in enumerate(snp.alleles) if snp_a == nt])

                    # If nt is not any of the possible alleles in SNP, it should be considered as invalid.
                    if len(alleles_set) == 0:
                        continue

                    compatible_set &= alleles_set
        return compatible_set

    def _calc_se_align_stat(self, aln: HTSeq.SAM_Alignment) -> AlignStat:
        """Calculate stats for single-end SAM alignment."""
        matched_regions, skipped_regions = self._extract_matched_and_skipped_regions(aln)
        compatible_isoforms = self._calc_compatible_isoforms(matched_regions, skipped_regions)
        compatible_alleles = self._calc_compatible_alleles(matched_regions)
        align_stat = AlignStat(aln=aln, aln2=None,
                               matched_ivs=[region.iv for region in matched_regions],
                               compatible_isoforms=compatible_isoforms,
                               compatible_alleles=compatible_alleles)
        return align_stat

    def _calc_pe_align_stat(self, aln1: HTSeq.SAM_Alignment, aln2: HTSeq.SAM_Alignment) -> AlignStat:
        """Calculate stats for paired-end-end SAM alignments pair."""
        aln1_stat = self._calc_se_align_stat(aln1)
        aln2_stat = self._calc_se_align_stat(aln2)
        assert aln1_stat.aln2 is None
        assert aln2_stat.aln2 is None

        aln_stat = AlignStat(aln=aln1_stat.aln, aln2=aln2_stat.aln,
                             matched_ivs=aln1_stat.matched_ivs + aln2_stat.matched_ivs,
                             compatible_isoforms=aln1_stat.compatible_isoforms & aln2_stat.compatible_isoforms,
                             compatible_alleles=aln1_stat.compatible_alleles & aln2_stat.compatible_alleles)
        return aln_stat

    def align_stats_iterator(self):
        """Iterator for (valid) alignment stats."""
        for bam_path in self.bam_param.paths:
            bam_reader = HTSeq.BAM_Reader(bam_path)
            alns = bam_reader.fetch(self.segment.iv.chrom, self.segment.iv.start, self.segment.iv.end)
            if self.bam_param.paired_end:
                for aln1, aln2 in HTSeq.pair_SAM_alignments_with_buffer(alns):
                    if (aln1 is None) or (aln2 is None):
                        continue
                    if self._validate_aln(aln1) and self._validate_aln(aln2):
                        aln_stat = self._calc_pe_align_stat(aln1, aln2)
                        if aln_stat.invalid:
                            continue
                        yield aln_stat
            else:
                for aln in alns:
                    if aln is None:
                        continue
                    if self._validate_aln(aln):
                        aln_stat = self._calc_se_align_stat(aln)
                        if aln_stat.invalid:
                            continue
                        yield aln_stat
