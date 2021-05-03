import gzip
import Bio.SeqIO.QualityIO
from Bio.Seq import Seq
import Bio.SeqIO.FastaIO as FastaIO
from typing import Generator


class FastA(object):
    '''
    Standard data container for fasta sequences
    '''
    __slots__ = ['header', 'sequence']
    def __init__(self, header: str, sequence: str) -> None:
        self.header = header
        self.sequence = sequence



def revcomp(sequence: str) -> str:
    '''
    Reverse complement a standard nucleotide sequence.
    :param sequence:
    :return:
    '''
    return str(Seq(sequence).reverse_complement())


def stream_fa(sequence_file: str) -> Generator[FastA, None, None]:
    '''
    Read a fastq file either gzipped or not and return it as a stream of tuples
    (Header, Sequence, Quality)
    :param infile:
    :return: Generator[FastA, None, None]
    '''

    if sequence_file.endswith('fq.gz') or sequence_file.endswith('fastq.gz'):
        with gzip.open(sequence_file, 'rt') as handle:
            for header, sequence, qual in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle):
                yield FastA(header, sequence)
    elif sequence_file.endswith('fq') or sequence_file.endswith('fastq'):
        with open(sequence_file) as handle:
            for header, sequence, qual in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle):
                yield FastA(header, sequence)
    elif sequence_file.endswith('fasta.gz') or sequence_file.endswith('fa.gz'):
        with gzip.open(sequence_file, 'rt') as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield FastA(header, sequence)
    elif sequence_file.endswith('fasta') or sequence_file.endswith('fa'):
        with open(sequence_file) as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield FastA(header, sequence)
    else:
        raise Exception(f'{sequence_file} not a sequence file.')
