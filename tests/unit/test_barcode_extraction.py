import subprocess
import shlex
from tnseq2.src import mapping, sequence

import logging

logging.basicConfig(filename='tests/test_barcode_extraction.log', filemode='w', level=logging.INFO)

import pickle

TESTDATA = "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/test_data/tnseq/"


def get_test_data():
    r1 = f'{TESTDATA}/mapping/library_13_1_1.fq'
    r2 = f'{TESTDATA}/mapping/library_13_1_2.fq'
    return r1, r2


def test_extract_barcode_host1():
    r1, r2 = get_test_data()
    tnSeq = 'GTGTATAAGAGACAG'
    r = list(sequence.stream_fa(r1))[1000]
    expectedBarcode = 'GACGGCTATACTCAAAG'
    expectedHost = 'GTCCTGAACTCGCGACGCAAATAAACGGTATCTTGAGGTATCTAATAACGAGAAATGCTTTAAAATATTAATTTCCGGGTAAGCATCATCATCATAGAAAAATAC'
    outBarcode, outHost = mapping.extract_barcode_host(r, tnSeq)
    assert outBarcode == expectedBarcode
    assert outHost == expectedHost


def test_extract_barcodes():
        r1, r2 = get_test_data()
        reads1 = sequence.stream_fa(r1)
        reads2 = sequence.stream_fa(r2)
        tp2= 'GTGTATAAGAGACAG'
        bc2tp2=13
        bcLen=17
        before=True
        mode='tn5'
        min_host_bases = 20
        barcodes = mapping.extract_barcodes(zip(reads1, reads2), tp2)
        with open(f'{TESTDATA}/mapping/test_extract_barcodes.out', 'rb') as po:
            expectedBarcode = pickle.load(po)
        assert barcodes == expectedBarcode

def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,)
    out, err = proc.communicate()
    return out, err, proc

def test_prepare_extract_barcodes():
    r1, r2 = get_test_data()
    outFasta = f'{TESTDATA}/mapping/tmp/pepare_extract_barcodes.out'
    expectedFasta = f'{TESTDATA}/mapping/pepare_extract_barcodes.out'
    mapping.prepare_extract_barcodes(r1, r2, outFasta)
    cmd_str2 = f'cmp --silent  {outFasta} {expectedFasta}'
    out, err, proc = capture(cmd_str2)
    assert proc.returncode == 0
