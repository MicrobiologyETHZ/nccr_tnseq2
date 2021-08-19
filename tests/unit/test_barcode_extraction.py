import subprocess
import shlex
from tnseq2.src import sequence, extract_barcodes

TESTDATA = "./tests/test_files"


def get_test_data():
    r1 = f'{TESTDATA}/library_13_1_1.fq'
    r2 = f'{TESTDATA}/library_13_1_2.fq'
    return r1, r2


def test_extract_barcode_host():
    r1, r2 = get_test_data()
    tp2 = 'GTGTATAAGAGACAG'
    r = list(sequence.stream_fa(r1))[1000]
    expectedBarcode = 'GACGGCTATACTCAAAG'
    expectedHost = 'GTCCTGAACTCGCGACGCAAATAAACGGTATCTTGAGGTATCTAATAACGAGAAATGCTTTAAAATATTAATTTCCGGGTAAGCATCATCATCATAGAAAAATAC'
    outBarcode, outHost = extract_barcodes.extract_barcode_host(r, tp2)
    assert outBarcode == expectedBarcode
    assert outHost == expectedHost


# todo add tests for WISH barcode structure
