import pysam
from typing import Generator


def filter_by_alignment_length(alignments: Generator[pysam.AlignedSegment, None, None], min_aln_length: int) -> Generator[pysam.AlignedSegment, None, None]:
    """ Takes a generator of alignments (as AlignedSegments) and removes
        alignments that have an alignment length shorter then min_aln_length

        Requires the 'al' tag set in the AlignedSegment. Will fail otherwise.

    :param alignments: A collection of alignments. Preferably a generator
    :param min_aln_length: The minimum length of an alignment to pass this filter
    :return: A generator of AlignedSegments
    """

    assert min_aln_length >= 0
    return filter(lambda aln: aln.get_tag('al') >= min_aln_length, alignments)


def filter_by_query_coverage(alignments: Generator[pysam.AlignedSegment, None, None], min_query_coverage: float) -> Generator[pysam.AlignedSegment, None, None]:
    """ Takes a generator of alignments (as AlignedSegments) and removes
        alignments that have a query coverage lower then min_query_coverage

        Requires the 'qc' tag set in the AlignedSegment. Will fail otherwise.

    :param alignments: A collection of alignments. Preferably a generator
    :param min_query_coverage: The minimum query coverage to pass the filter
    :return: A generator of AlignedSegments
    """

    assert min_query_coverage >= 0.0
    assert min_query_coverage <= 1.0
    return filter(lambda aln: aln.get_tag('qc') >= min_query_coverage, alignments)


def filter_unmapped(alignments: Generator[pysam.AlignedSegment, None, None]) -> Generator[pysam.AlignedSegment, None, None]:
    """ Takes a generator of alignments (as AlignedSegments) and removes
        alignments that are unaligned using the SAM flag 4.

    :param alignments: A collection of alignments. Preferably a generator
    :return: A generator of AlignedSegments
    """

    return filter(lambda aln: aln.is_unmapped is not True, alignments)


def filter_by_percent_identity(alignments: Generator[pysam.AlignedSegment, None, None], min_percid: float) -> Generator[pysam.AlignedSegment, None, None]:
    """ Takes a generator of alignments (as AlignedSegments) and removes
        alignments that have a percent identity lower then min_percid

        Requires the 'id' tag set in the AlignedSegment. Will fail otherwise.

    :param alignments: A collection of alignments. Preferably a generator
    :param min_percid: The minimum percent identity that an alignment has to have pass this filter. Range between 0.0 and 1.0
    :return: A generator of AlignedSegments
    """
    assert min_percid <= 1.0
    assert min_percid >= 0.0

    return filter(lambda aln: aln.get_tag('id') >= min_percid, alignments)



def add_features(alignments: Generator[pysam.AlignedSegment, None, None]) -> Generator[pysam.AlignedSegment, None, None]:
    """ Adds Percent ID and Alignment length to a collection of AlignedSegments

    :param alignments: A collection of alignments. Preferably a generator
    :return: alignments: A collection of alignments.
    """

    return map(lambda aln : add_tags(aln), alignments)



def add_tags(alignedSegment: pysam.AlignedSegment) -> pysam.AlignedSegment:
    """ Takes an AlignedSegment and add percent identity and alignment length as tags
    alignment length = MID
    mismatches = NM
    percent identity = (MID - NM) / MID

    The percent identity is a value between 0.0 and 1.0

    If the segment is unmapped then it is returned as with a percent identity of 0
    and an alignment length of 0.

    :param alignedSegment: The pysam AlignedSegment object
    :return: alignedSegment: The updated pysam AlignedSegment object
    """
    if alignedSegment.is_unmapped:
        alignedSegment.set_tag('id', 0.0, 'f')
        alignedSegment.set_tag('al', 0, 'i')
        alignedSegment.set_tag('qc', 0.0, 'f')
        return alignedSegment
    alnlength = sum(map(lambda i : alignedSegment.get_cigar_stats()[0][i], range(0, 3)))
    query_covered_bases = sum(map(lambda i: alignedSegment.get_cigar_stats()[0][i], range(0, 2)))
    query_length = alignedSegment.infer_read_length()
    mismatches = alignedSegment.get_tag('NM')
    percid = (alnlength - mismatches) / float(alnlength)
    qcov = query_covered_bases / float(query_length)
    alignedSegment.set_tag('id', percid, 'f')
    alignedSegment.set_tag('qc', qcov, 'f')
    alignedSegment.set_tag('al', alnlength, 'i')
    return alignedSegment


