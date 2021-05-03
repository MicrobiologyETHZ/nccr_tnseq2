#from tnseq import sequence
#from tnseq import tnseq
#from tnseq import extract_barcodes

from tnseq2.src.mapping import *
import subprocess
import shlex

TESTDATA = "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/test_data/tnseq/"


def get_test_data():
    r1 = f'{TESTDATA}/mapping/library_13_1_1.fq'
    r2 = f'{TESTDATA}/mapping/library_13_1_2.fq'
    genome = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    expectedMap = f'{TESTDATA}/mapping/library_13_1_barcode_ed.txt'
    outMap = f'{TESTDATA}/mapping/tmp/library_13_1_barcode.txt'
    return r1, r2, genome, expectedMap, outMap

def get_test_data_map():
    fasta_file = f'{TESTDATA}/mapping/tnseq_host_reads.fasta'
    blastn_file = f'{TESTDATA}/mapping/tnseq_host_reads.blastn'
    Sdb = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    return fasta_file, blastn_file, Sdb


def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,)
    out, err = proc.communicate()
    return out, err, proc

# def test_map_host():
#     fasta_file, blastn_file, Sdb = get_test_data_map()
#     expectedOut = f'{TESTDATA}/mapping/map_host.out'
#     actualOut = f'{TESTDATA}/mapping/tmp/map_host.out'
#     mapping.map_host(fasta_file, blastn_file, Sdb, actualOut)
#     cmd_str2 = f'cmp --silent  {actualOut} {expectedOut}'
#     out, err, proc = capture(cmd_str2)
#     assert proc.returncode == 0


def test_map_host_map():
    #r1, r2, genome, expectedMap, outMap = get_test_data()
    fasta_file, blastn_file, genome = get_test_data_map()
    actualBlastn = f'{TESTDATA}/mapping/tmp/blastn.out'
    expectedBlastn = blastn_file
    map_host_map(fasta_file, actualBlastn, genome)
    #cmd_str2 = f'cmp --silent  {actualBlastn} {expectedBlastn}'
    #out, err, proc = capture(cmd_str2)
    #assert proc.returncode == 0

test_map_host_map()


def test_map_host_new():
    fasta_file, blastn_file, Sdb = get_test_data_map()
    barcode_file = f'{TESTDATA}/new_map.csv'
    map_host_new(fasta_file, blastn_file, Sdb, 1, barcode_file)



#
# def test_barcode_to_host_sequences():
#     r1, r2, Sdb, outMap, expectedMap = get_test_data()
#     logging.basicConfig(filename ="/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/mapBarcodeToHost.log", format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)
#     inserts = zip(sequence.stream_fa(r1), sequence.stream_fa(r2))
#     alignments = map_barcodes.barcode_to_host_sequences(inserts)
#     with open("/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/mapBarcodeToHost.txt", 'rb') as fp:
#         expectedAlignments = pickle.load(fp)
#     assert alignments == expectedAlignments
#
#
# def test_barcodes_host_to_fasta():
#     r1, r2, Sdb, outMap, expectedMap = get_test_data_map()
#     fasta_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/temp.fasta"
#     bar_2_seq = map_barcodes.barcodes_host_to_fasta(r1, r2, fasta_file)
#     with open("/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/barcodes_host_to_fasta.txt", 'rb') as fp:
#         expectedDict = pickle.load(fp)
#     assert expectedDict == bar_2_seq # this does not check the fasta file




#



#
# def test_extract_best_hits():
#     expected_dict = {'1_bc_AAGAGTGCTGGATCCTG_cnt_75': [['FQ312003.1', '2393467', 'minus']],
#                      '2_bc_AAGAGTGCTGGATCCTG_cnt_52': [['FQ312003.1', '2393467', 'minus']],
#                      '3_bc_AAGAGTGCTGGATCCTG_cnt_9': [['FQ312003.1', '2393467', 'minus']],
#                      '4_bc_AAGAGTGCTGGATCCTG_cnt_9': [['FQ312003.1', '2393467', 'minus']],
#                      '5_bc_AAGAGTGCTGGATCCTG_cnt_8': [['FQ312003.1', '2393467', 'minus']],
#                      '144264_bc_ATAGGGCATCACATATC_cnt_1': [['FQ312003.1', '68415', 'minus']],
#                      '144265_bc_TCTCCAAGTTGACGGTT_cnt_1': [['FQ312003.1', '3109727', 'minus']],
#                      '144266_bc_CCGCGTTGACGCTGGAC_cnt_1': [['FQ312003.1', '1174483', 'minus']],
#                      '144267_bc_AGCCTTGCTGCTGGTTT_cnt_1': [['FQ312003.1', '358092', 'plus']],
#                      '144268_bc_CGCGGGTAAAGGACGAG_cnt_1': [['NC_017720.1', '27484', 'plus']],
#                      '246_bc_CGATACCCGTAACGCGT_cnt_3': [['NC_017720.1', '33037', 'plus']],
#                      '247_bc_CGATACCCGTAACGCGT_cnt_3': [['NC_017720.1', '33037', 'plus']],
#                      '248_bc_CGATACCCGTAACGCGT_cnt_3': [['NC_017720.1', '33037', 'plus']],
#                      '249_bc_CGATACCCGTAACGCGT_cnt_3': [['NC_017720.1', '33037', 'plus']],
#                      '250_bc_CGATACCCGTAACGCGT_cnt_3': [['NC_017720.1', '33037', 'plus']]}
#     blast_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/expected.blastn"
#     q2besthits, _ = map_barcodes.extract_best_hits(blast_file)
#     assert expected_dict == q2besthits

#test_extract_best_hits()

#
# def test_is_simiar_barcode():
#     barcode1 = ''
#     barcode2 = ''
#     return None
#
#
# def test_map_annotate():
#     fasta_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/test.fasta"
#     blast_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/expected.blastn"
#     header_2_sequence, query_2_besthits, _ = map_barcodes.map_host_load(fasta_file, blast_file)
#     barcode_2_potential_position, barcode_2_count = map_barcodes.map_annotate(header_2_sequence, query_2_besthits)
#     expected_b2position, expected_bc2_count = ({'AAGAGTGCTGGATCCTG': ((2393467, 'FQ312003.1', 'minus'), 153.0),
#                            'CGATACCCGTAACGCGT': ((33037, 'NC_017720.1', 'plus'), 15.0),
#                            'ATAGGGCATCACATATC': ((68415, 'FQ312003.1', 'minus'), 1.0),
#                            'TCTCCAAGTTGACGGTT': ((3109727, 'FQ312003.1', 'minus'), 1.0),
#                            'CCGCGTTGACGCTGGAC': ((1174483, 'FQ312003.1', 'minus'), 1.0),
#                            'AGCCTTGCTGCTGGTTT': ((358092, 'FQ312003.1', 'plus'), 1.0),
#                            'CGCGGGTAAAGGACGAG': ((27484, 'NC_017720.1', 'plus'), 1.0)}, Counter({'AAGAGTGCTGGATCCTG': 153, 'CGATACCCGTAACGCGT': 15,
#                                   'ATAGGGCATCACATATC': 1, 'TCTCCAAGTTGACGGTT': 1,
#                                   'CCGCGTTGACGCTGGAC': 1, 'AGCCTTGCTGCTGGTTT': 1,
#                                   'CGCGGGTAAAGGACGAG': 1}))
#
#     assert barcode_2_potential_position == expected_b2position
#     assert barcode_2_count == expected_bc2_count
#
# def test_map_host_barcode_merge():
#     fasta_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/test.fasta"
#     blast_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/expected.blastn"
#     header_2_sequence, query_2_besthits, _ = map_barcodes.map_host_load(fasta_file, blast_file)
#     barcode_2_potential_position, barcode_2_count = map_barcodes.map_annotate(header_2_sequence, query_2_besthits)
#     out_barcode_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/expected.barcode.out"
#     barcode_max_edit_distance = 0
#     map_barcodes.map_host_barcode_merge(barcode_2_count, barcode_2_potential_position, out_barcode_file, barcode_max_edit_distance)
#


# def test_multiple_best_hits():
#     r1, r2, Sdb, outMap, expectedMap = get_test_data_map()
#     fasta_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/multiple_blast_hits.fasta"
#     blast_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/multiple_blast_hits.blastn"
#     #map_barcodes.blast_host(fasta_file, blast_file, Sdb) # This was created with blastn-short
#     header_2_sequence, query_2_besthits = map_barcodes.map_host_load(fasta_file, blast_file)
#     barcode_2_potential_position, barcode_2_count = map_barcodes.map_annotate(header_2_sequence, query_2_besthits)
#     return barcode_2_potential_position, barcode_2_count



# def test_one_bc_multiple_locations():
#     r1, r2, Sdb, outMap, expectedMap = get_test_data_map()
#     fasta_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/test_one_bc_two_locations.fasta"
#     blast_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/test_one_bc_two_locations.blastn"
#     map_barcodes.blast_host(fasta_file, blast_file, Sdb) # This was created with blastn-short
#     header_2_sequence, query_2_besthits = map_barcodes.map_host_load(fasta_file, blast_file)
#     barcode_2_potential_position, barcode_2_count = map_barcodes.map_annotate(header_2_sequence, query_2_besthits)
#     return header_2_sequence, query_2_besthits




#print(test_multiple_best_hits())

# def test_map_annotate():
#     fasta_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/test.fasta"
#     blast_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/expected.blastn"
#
#     header_2_sequence, query_2_besthits = map_barcodes.map_host_load(fasta_file, blast_file)
#     breakpoint()
#     a, b = map_barcodes.map_annotate(header_2_sequence, query_2_besthits)
#     return a, b

#print(test_map_annotate())

