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

from typing import List

from ..component import Segment, SNP


class Task:
    """Task.

    Attributes:
        id: Task ID.
        segment: The target segment.
        snps: The SNPs located on the segment.
        ploidy: The ploidy of the SNP.
        phased: If the task is phased.
    """

    def __init__(self, task_id: str, segment: Segment, phased: bool, ploidy: int, snps: List[SNP]):
        self.id = task_id
        self.segment = segment
        self.phased = phased
        self.ploidy = ploidy
        self.snps = snps

        if not self.phased:
            assert len(self.snps) == 1

        for snp in self.snps:
            assert self.phased == snp.phased
            assert self.ploidy == snp.ploidy

    @property
    def isoforms_count(self):
        """Isoforms count."""
        return len(self.segment.isoforms)

    def isoform_snps(self, isoform_num: int) -> List[SNP]:
        """SNPs located on isoform."""
        isoform = self.segment.isoforms[isoform_num]
        snps = []
        for snp in self.snps:
            for exon in isoform.exons:
                if exon.start <= snp.pos <= exon.end:
                    snps.append(snp)
                    break
        return snps
