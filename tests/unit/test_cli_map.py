import subprocess
import shlex
from pathlib import Path


TESTDATA = "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/test_data/tnseq/"
def get_test_data_map():

    r1 = Path(TESTDATA)/"mapping/library_13_1_1.fq"
    r2 = Path(TESTDATA)/"mapping/library_13_1_2.fq"
    Sdb = Path(TESTDATA)/"ref/Salmonella_genome_FQ312003.1_SL1344.fasta"
    expectedMap = Path(TESTDATA)/"mapping/library_13_1_barcode_ed.txt"
    outMap = Path(TESTDATA)/"mapping/tmp/library_13_1_barcode.txt"
    return r1, r2, Sdb, expectedMap, outMap


# def get_test_data_quant():
#     r = "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_quantify/1315-107-library11_1-TV3371B-inoculum.fasta.gz"
#     outMap = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testMap.tsv"
#     expectedQuant = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/outQuant.tsv"
#     outQuant = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testQuant.tsv"
#     return r, outMap, expectedQuant, outQuant

def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,)
    out, err = proc.communicate()
    return out, err, proc


def test_cli_map():
    r1, r2, Sdb, expectedMap, outMap = get_test_data_map()
    cmd_str = f'python tnseq2/tnseq.py map -r1 {r1} -r2 {r2} -o {outMap} -db {Sdb} -f '
    _, err, proc = capture(cmd_str)
    if proc.returncode !=0:
        print(to_str(err))
    assert proc.returncode == 0
    cmd_str2 = f'cmp --silent  {outMap} {expectedMap}'
    out, err, proc = capture(cmd_str2)
    assert proc.returncode == 0



def to_str(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value

#
# def test_cli_quantify():
#
#     r, outMap, expectedQuant, outQuant = get_test_data_quant()
#     cmd_str = f'tnseq quantify -r {r}  -o {outQuant} -b {outMap} -f'
#     _, err, proc = capture(cmd_str)
#     assert proc.returncode == 0
#     out, err, proc = capture(f'sort -k1 {outQuant}')
#     with open(expectedQuant, 'r') as o:
#        expectedLines = o.read()
#     assert to_str(out) == expectedLines


if __name__ == "__main__":
    test_cli_map()