import tnseq2.src.sequence
import tnseq2.src.commandline
import tnseq2.src.mapping
import logging
import collections
import sys
import pathlib
import argparse
from argparse import RawTextHelpFormatter
from typing import Tuple, List, Dict, Counter, Set


def mutate(sequence_to_mutate: str) -> List[str]:
    '''
    Receives a nucleotide sequence and returns a list with all possible 1 base mutations
    of the initial nucleotide sequence
    '''
    mutated_sequences: List[str] = []

    length: int = len(sequence_to_mutate)
    for pos in range(0, length):
        prefix: str = sequence_to_mutate[:pos]
        suffix: str = sequence_to_mutate[pos + 1:]
        mutated_sequences.append(prefix + 'A' + suffix)
        mutated_sequences.append(prefix + 'C' + suffix)
        mutated_sequences.append(prefix + 'T' + suffix)
        mutated_sequences.append(prefix + 'G' + suffix)
    return mutated_sequences

#
# def editdistance(seq1: str, seq2: str) -> int:
#     '''
#     Calculate the edit distance between 2 sequences with identical length.
#     Will throw an error if the length of both sequences differs
#     :param seq1:
#     :param seq2:
#     :return:
#     '''
#
#     if seq1 == seq2:
#         return 0
#     if len(seq1) != len(seq2):
#         raise Exception(f'{seq1} and {seq2} have different length. Edit distance can be computed on same length sequences only.')
#     dist = 0
#     for letter1, letter2 in zip(seq1, seq2):
#         if letter1 != letter2:
#             dist += 1
#     return dist



# def extract_barcode_host(r1: tnseq2.sequence.FastA, tp2: str) -> Tuple[str, str]:
#     '''
#    Extract barcode and host sequence from read with tp2.
#    Return (None, None) if the barcode sequence is not complete (17bp)
#     :param r1:
#     :param tp2:
#     :return:
#     '''
#
#     splits: List[str] = r1.sequence.split(tp2)
#     if len(splits[0]) >= 30:
#         barcode: str = splits[0][-30:-13]
#         host_sequence: str = splits[1]
#         return barcode, host_sequence
#     else:
#         return None, None


#
#
# def extract_barcodes(inserts: Generator[Tuple[tnseq2.sequence.FastA, tnseq2.sequence.FastA], None, None], mode: str='tn5', min_host_bases: int = 20) -> List[Tuple[str, str]]:
#     '''
#     |r1=(0-k)bp|BARCODE=17bp|tp1=13bp|tp2=15bp|r2=8bp|HOST=(15-j)bp||ADAPTER=(0-n)bp|
#     |----------|============|--------|========|------|=============||---------------|
#     :param inserts:
#     :param mode:
#     :return:
#     '''
#     if mode == 'tn5':
#         tp2 = 'GTGTATAAGAGACAG'
#         tp2_mut = mutate(tp2)
#     else:
#         logging.error('Unkown transposon')
#         shutdown(1)
#
#
#     total_inserts: int = 0
#     inserts_with_tp2: int = 0
#     inserts_without_tp2: int = 0
#     inserts_with_tp2_but_short_bc: int = 0
#     inserts_with_tp2_with_good_bc: int = 0
#     inserts_with_tp2_with_good_bc_short_host: int = 0
#     tp2mut = 0
#
#     sequences: List[Tuple[str, str]] = []
#     r1: tnseq2.sequence.FastA = None
#     r2: tnseq2.sequence.FastA = None
#
#     for total_inserts, (r1, r2) in enumerate(inserts):
#         if total_inserts % 100000 == 0:
#             logging.info(f'Processed {total_inserts} reads')
#         if tp2 in r1.sequence:
#             inserts_with_tp2 += 1
#             barcode, host_sequence = extract_barcode_host(r1, tp2)
#             if barcode:
#                 inserts_with_tp2_with_good_bc += 1
#                 if len(host_sequence) < min_host_bases:
#                     inserts_with_tp2_with_good_bc_short_host += 1
#                 else:
#                     sequences.append((barcode, host_sequence))
#             else:
#                 inserts_with_tp2_but_short_bc += 1
#
#         else:
#             inserts_without_tp2 += 1
#
#
#     logging.info(f'Processed {total_inserts} reads')
#     logging.info('Extraction statistics:')
#     logging.info(f'\tTotal inserts:\t{total_inserts}\t100.0%')
#     logging.info(f'\tInserts w transposon:\t{inserts_with_tp2}\t{int(100 * ((inserts_with_tp2 * 100.0)/total_inserts))/100.0}%')
#     logging.info(f'\t\tand w barcode:\t{inserts_with_tp2_with_good_bc}\t{int(100 * ((inserts_with_tp2_with_good_bc * 100.0)/total_inserts))/100.0}%')
#     logging.info(f'\t\t\tand good host sequence:\t{len(sequences)}\t{int(100 * ((len(sequences) * 100.0)/total_inserts))/100.0}% --> Used for downstream analysis')
#     logging.info(f'\t\t\tand short host sequence:\t{inserts_with_tp2_with_good_bc_short_host}\t{int(100 * ((inserts_with_tp2_with_good_bc_short_host * 100.0)/total_inserts))/100.0}%')
#     logging.info(f'\t\tand w/o barcode:\t{inserts_with_tp2_but_short_bc}\t{int(100 * ((inserts_with_tp2_but_short_bc * 100.0)/total_inserts))/100.0}%')
#     logging.info(f'\tInserts w/o transposon:\t{inserts_without_tp2}\t{int(100 * ((inserts_without_tp2 * 100.0)/total_inserts))/100.0}%')
#
#     return sequences
#
#
#
#



#
#
# def prepare_extract_barcodes(in_r1_file: str, in_r2_file: str, temp_fasta_file: str, mode: str= 'tn5') -> None:
#     '''
#     Wrapper function to extract barcodes and host sequences from
#     sequence file. Sequences will then be written to fasta file
#     for downstream analysis
#     :param in_r1_file:
#     :param in_r2_file:
#     :param temp_fasta_file:
#     :param mode:
#     :return:
#     '''
#
#
#     fq1_stream: Generator[tnseq2.sequence.FastA, None, None] = tnseq2.sequence.stream_fa(in_r1_file)
#     fq2_stream: Generator[tnseq2.sequence.FastA, None, None] = tnseq2.sequence.stream_fa(in_r2_file)
#     logging.info('----------------')
#     logging.info('Step 1.1: Start extraction')
#
#     inserts: Iterator[Tuple[tnseq2.sequence.FastA, tnseq2.sequence.FastA]]= zip(fq1_stream, fq2_stream)
#     alignments = extract_barcodes(inserts, mode)
#
#     logging.info('Step 1.1: Finished extraction')
#     logging.info('----------------')
#     logging.info('Step 1.2: Start writing dereplicated barcode/host pairs.')
#     with open(temp_fasta_file, 'w') as handle:
#         barcode_2_sequences = collections.defaultdict(list)
#         for alignment in alignments:
#             barcode_2_sequences[alignment[0]].append(alignment[1])
#         tot = 1
#         for barcode, sequences in barcode_2_sequences.items():
#             for sequence, cnt in  collections.Counter(sequences).most_common():
#                 handle.write(f'>{tot}_bc_{barcode}_cnt_{cnt}\n{sequence}\n')
#                 tot += 1
#     logging.info(f'Step 1.2: Finished writing {tot-1} dereplicated barcode/host pairs.')
#     logging.info('----------------')
#
#     # with open(out_barcode, 'w') as handle:
#     #     cnter = collections.Counter(map(lambda x: x[0], alignments))
#     #     for bc, cnt in cnter.most_common():
#     #         handle.write(f'{bc}\t{cnt}\n')
#
#
#



def quantify_read_barcode_map_files(in_bc_file: str) -> Tuple[Dict[str, Tuple[int, str, str]], Dict[str, float]]:
    '''
    Reading barcode to position file into memory

    This only works with a file that looks the following:
    BARCODE             COUNT   STARTPOS    CHROMOSOME  ORIENTATION NORM_COUNT

    ACGCAGACCCTCACTTT       24579   5151    NC_017719.1     minus   24485.0
    CAGGTACTCAGACAACG       14807   7883    NC_017719.1     minus   14753.0
    CCCGAACCCCTGGCAGT       11731   4773    NC_017719.1     minus   11691.0
    AAAACCTCCCTGCCCAT       11650   7626    NC_017719.1     minus   11598.0
    GCAAAAGGCCATAAATG       8894    4095    NC_017719.1     minus   8828.166666666668
    GCATTCTCACCGGTCGA       6198    4776    NC_017719.1     plus    6164.0
    ACCGAAGGTACCCGTAT       4943    4188    NC_017719.1     minus   4911.333333333333
    AAAAGACTCTTTAGCCC       4705    4532    NC_017719.1     minus   4685.0
    GTAGGAGGTCTCAAAAA       4530    3113    NC_017719.1     minus   4508.0
    GGACCAGCGTATACAAC       2478    7071    NC_017719.1     plus    2467.0
    :param in_bc_file:
    :return:
    '''
    barcode_2_pos: Dict[str, Tuple[int, str, str]] = {}
    barcode_2_abundance:  Dict[str, float] = collections.Counter()
    with open(in_bc_file) as handle:
        for line in handle:
            line: str = line.strip()
            splits: List[str] = line.split()
            barcode_2_pos[splits[0]] = (int(splits[2]), splits[3], splits[4])
            barcode_2_abundance[splits[0]] = float(splits[5])
    return barcode_2_pos, barcode_2_abundance



def prefix(sequence, minlength=15):
    prefixes = []
    for j in range(minlength,len(sequence)+1):
            prefixes.append(sequence[:j])
    return prefixes


def quantify_load_fq_barcodes(in_fq_file: str, tp2: str='GTGTATAAGAGACAG', pad: int=0) -> Counter[int]:
    '''
    Load the barcode counts from a fasta/fastq file.
    This currently only works with tn5.
    |r1=(0-k)bp|BARCODE=17bp|tp1=13bp|tp2=15bp|rest=(0-j)bp|
    |----------|============|--------|========|------------|


    :param in_fq_file:
    :param tp2:
    :return:
    '''

    fq1_stream = tnseq2.src.sequence.stream_fa(in_fq_file)
    with_tp2: int = 0
    without_tp2: int = 0
    with_tp_but_short: int = 0
    cnter: Counter[str] = collections.Counter()
    for total_inserts, r1 in enumerate(fq1_stream, 1):
        if total_inserts % 1000000 == 0:
            logging.info(f'\tReads processed:\t{total_inserts}')
        if tp2 in r1.sequence:
            startpos_tp2 = r1.sequence.index(tp2)
            if startpos_tp2 != 30 + pad:
                with_tp_but_short += 1
                continue
            with_tp2 += 1
            barcode = r1.sequence[:17]
            cnter[barcode] += 1
        else:
            without_tp2 += 1
    logging.info(f'\tReads processed:\t{total_inserts}')
    logging.info(f'\tFastA/Q Stats:')
    logging.info(f'\t\tTotal Reads:\t{total_inserts}')
    logging.info(f'\t\tReads w transposon and good barcode:\t{with_tp2}')
    logging.info(f'\t\tReads w transposon but short barcode:\t{with_tp_but_short}')
    logging.info(f'\t\tReads w/o transposon:\t{without_tp2}')

    return cnter


def quantify_extract_annotated_correct(barcode_2_position, barcode_2_abundance, in_fq_file, out_barcode_count_file, max_edit_distance=2, mode='tn5', pad: int=0):
    '''
    Method to extract barcodes from sequence file and match them with known
    barcodes from the map step.

    1. Take a sequence file and iterate. Find complete barcodes and count them
    2. Compare these barcodes against known barcodes. Start with edit distance = 0, then increase
        till max_edit_distance is reached.
    3. Print matched and unmatched barcodes together with abundances to the output file

    :param barcode_2_position:
    :param barcode_2_abundance:
    :param in_fq_file:
    :param out_barcode_count_file:
    :param max_edit_distance:
    :return:
    '''

    logging.info('----------------')
    logging.info('Step 2.1. Load barcodes from sequence file')
    cnter: Counter[str] = quantify_load_fq_barcodes(in_fq_file, pad=pad)
    logging.info(f'\tFound {len(cnter)} distinct barcodes in sequence file.')
    logging.info('Step 2.1. Finished')
    logging.info('----------------')

    logging.info('Step 2.2. Match barcodes against known barcodes')

    current_edit_distance: int = 0
    found_barcode_2_abundance: Counter[str] = collections.Counter()
    most_common: List[Tuple[str, int]] = barcode_2_abundance.most_common()
    unmatched_barcodes: Set[str] = set(cnter.keys())
    unmatched_barcodes_tmp: Set[str] = set()
    editdistance_count: Counter[int] = collections.Counter()

    while current_edit_distance <= max_edit_distance:
        logging.info(f'\tStart matching {len(unmatched_barcodes)} barcodes with edit distance {current_edit_distance} against known barcodes')

        for total_count, barcode in enumerate(unmatched_barcodes, 1):
            if total_count % 10000 == 0:
                logging.info(f'\t{total_count}/{len(unmatched_barcodes)} tested')
            if current_edit_distance == 0:
                    if barcode in barcode_2_position:
                        found_barcode_2_abundance[barcode] += cnter[barcode]
                        editdistance_count[current_edit_distance] += 1
                    else:
                        unmatched_barcodes_tmp.add(barcode)
            else:
                found: bool = False
                for known_barcode, known_barcode_count in most_common:
                    if editdistance(known_barcode, barcode) <= current_edit_distance:
                        found_barcode_2_abundance[known_barcode] += cnter[barcode]
                        editdistance_count[current_edit_distance] += 1
                        found = True
                        break
                if not found:
                    unmatched_barcodes_tmp.add(barcode)
        logging.info(f'\t{total_count}/{len(unmatched_barcodes)} tested')
        unmatched_barcodes = unmatched_barcodes_tmp
        unmatched_barcodes_tmp: Set[str] = set()
        current_edit_distance += 1
        logging.info(
            f'\tFinished matching barcodes with edit distance {current_edit_distance - 1}')
    editdistance_count[-1] = len(unmatched_barcodes)
    for edit_distance, count in editdistance_count.most_common():
        if edit_distance == -1:
            logging.info(f'\t{count} barcodes not matched with maximum edit distance {max_edit_distance}')
        else:
            logging.info(f'\t{count} barcodes matched with edit distance {edit_distance}')
    for unmatched_barcode in unmatched_barcodes:
        found_barcode_2_abundance[unmatched_barcode] += cnter[unmatched_barcode]

    logging.info(f'\tBarcode Stats:')
    logging.info(f'\t\tTotal annotated barcodes:\t{len(found_barcode_2_abundance) - len(unmatched_barcodes)}')
    logging.info(f'\t\tTotal unknown barcodes:\t{len(unmatched_barcodes)}')
    logging.info(f'\t\tTotal barcodes:\t{len(found_barcode_2_abundance)}')

    logging.info(f'\tFinished matching detected counting barcodes to annotated barcodes (max {max_edit_distance} edit distance)')
    logging.info('Step 2.2. Finished')
    logging.info('----------------')



    logging.info(f'Step 2.3. Start writing output file.')
    with open(out_barcode_count_file, 'w') as handle:
        for barcode, abundance in found_barcode_2_abundance.items():
            if barcode in barcode_2_position:
                tmp = '\t'.join(map(lambda x: str(x), barcode_2_position[barcode]))
                handle.write(f'{barcode}\t{abundance}\t{tmp}\n')
            else:
                handle.write(f'{barcode}\t{abundance}\tNone\n')
    logging.info('Step 2.3. Finished')
    logging.info('----------------')
    logging.info('Step 2. Finished')





# def parse_mapping(blastfile: str) -> Dict[str, List[List[str]]]:
#     '''
#     Load blastn file and extract for each query the best hit
#     by bitscore. Will return a list with one or multiple best
#     hits.
#     :param blastfile:
#     :return:
#     '''
#     query_2_hits: DefaultDict[str, List[List[str]]] = collections.defaultdict(list)
#     with open(blastfile) as handle:
#         for line in handle:
#             splits = line.strip().split('\t')
#             query_2_hits[splits[0]].append(splits)
#     query_2_besthits: Dict[str, List[List[str]]] = {}
#     for query, hits in query_2_hits.items():
#         bestscore: float = max(map(lambda x: float(x[9]), hits))
#         best_score_hits: List[List[str]] = list(map(lambda y: [y[1], y[6], y[11]], filter(lambda x : float(x[9]) >= bestscore, hits)))
#
#         query_2_besthits[query] = best_score_hits
#     return query_2_besthits
#

#
# def is_similar_barcode(barcode1: str, pos1: int, barcode2: str, pos2: int, barcode_max_edit_distance: int = 0) -> bool:
#     '''
#     Takes 2 barcodes and their positions in the host genomes as input
#     and checks if the barcodes have a edit distance smaller then
#     {barcode_max_edit_distance}. If so then they should map to
#     roughly the same position (+/- 5 bases) in the reference
#     genome. Then they are considered similar. Otherwise the
#     method returns false.
#     :param barcode1:
#     :param pos1:
#     :param barcode2:
#     :param pos2:
#     :param barcode_max_edit_distance:
#     :return:
#     '''
#     if editdistance(barcode1, barcode2) <= barcode_max_edit_distance:
#         if pos1[1] == pos2[1] and pos1[2] == pos2[2]:
#             if (pos1[0] - 5) <= pos2[0] <= (pos1[0] + 5):
#                 return True
#     return False

# def map_host_map(temp_fasta_file: str, temp_blastn_file: str, blastdb: str, blast_threads: int = 1) -> None:
#     '''
#     Map reads against the  blast database
#     :param temp_fasta_file:
#     :param temp_blastn_file:
#     :param blastdb:
#     :param out_barcode_file:
#     :param blast_threads:
#     :return:
#     '''
#     command = f'blastn -task blastn -db {blastdb} -out {temp_blastn_file} -query {temp_fasta_file} -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand" -num_threads {blast_threads}'
#     logging.info(f'Blastn command:\t{command}')
#     tnseq2.commandline.check_call(command)

# def map_host_load(temp_fasta_file: str, temp_blastn_file: str) -> Tuple[Dict[str, str], Dict[str, List[List[str]]]]:
#     '''
#     Helper method to load blast and fasta file
#     :param temp_fasta_file:
#     :param temp_blastn_file:
#     :return:
#     '''
#     fq1_stream = tnseq2.sequence.stream_fa(temp_fasta_file)
#     header_2_sequence: Dict[str, str] = {}
#     for read in fq1_stream:
#         header_2_sequence[read.header] = read.sequence
#     query_2_besthits: Dict[str, List[List[str]]] = parse_mapping(temp_blastn_file)
#     return header_2_sequence, query_2_besthits
#
#
#


# def map_annotate(header_2_sequence : Dict[str, str], query_2_besthits: Dict[str, List[List[str]]]) -> Tuple[Dict[str, Tuple[Tuple[str, str, str], int]], Counter[str]]:
#     '''
#     Find the best position/orientation for each barcode.
#     Each barcode can consist of multiple host sequences.
#     The correct host sequences can be identified the following:
#     - They appear in large quantity
#     - Even if they dont appear in large quantity they map against
#         the sample position in the host genome
#
#     For each barcode we look into the alignments for each associated
#     host sequence and pick the most likely position. The most likely
#     position is the one that appears most often.
#
#     :param header_2_sequence:
#     :param query_2_besthits:
#     :return:
#     '''
#     barcode_2_alignments: Dict[str, List[Tuple[int, List[str]]]] = {}
#     barcode_2_count: Counter[str] = collections.Counter()
#     barcode_2_potential_position: Dict[str, Tuple[Tuple[str, str, str], int]] = {}
#     for header in header_2_sequence:
#         splits: List[str] = header.split('_')
#         barcode: str = splits[2]
#         cnt: int = int(splits[4])
#         barcode_2_count[barcode] += cnt
#         barcode_2_alignments[barcode] = []
#
#     if len(query_2_besthits) == 0:
#         logging.error('There were 0 alignments found in the blast output found. Most likely a problem with the input file, the blast installation and/or blast database.')
#         shutdown(1)
#     for query, besthits in query_2_besthits.items():
#         splits: List[str] = query.split('_')
#         barcode: str = splits[2]
#         cnt: int = int(splits[4])
#         barcode_2_alignments[barcode].append((cnt, besthits))
#
#     # for each barcode i get all alignments
#     # Then group the alignments by start position (beg plus, end minus)
#     # (start, ref, orientation) = count
#     for barcode, total_count in barcode_2_count.most_common():
#         alignments: List[Tuple[int, List[str]]] = barcode_2_alignments[barcode]
#         if len(alignments) == 0:
#             barcode_2_potential_position[barcode] = (('-1', 'None', 'plus'), total_count)
#             continue
#         potential_positions: Counter[Tuple[str, str, str]] = collections.Counter()
#         for alignment in alignments:
#             alignment_count: str = alignment[0]
#             seq_alignments: List[List[str]] = alignment[1]
#             weight = float(alignment_count) / len(seq_alignments)
#             for seqaln in seq_alignments:
#                 orientation: str = seqaln[2]
#                 reference: str = seqaln[0]
#                 startpos: int = int(seqaln[1])
#                 potential_positions[(startpos, reference, orientation)] += weight
#         mc: Tuple[Tuple[str, str, str], int] = potential_positions.most_common()[0]
#         barcode_2_potential_position[barcode] = mc
#
#     return barcode_2_potential_position, barcode_2_count
#


# def map_host_barcode_merge(barcode_2_count, barcode_2_potential_position, out_barcode_file, barcode_max_edit_distance: int=0):
#     '''
#     Gets a list of barcodes/abundance and their respective postion in the
#     host genome. Barcodes are then sorted by abundance as more abundant
#     barcode are potentially more credible. The algorithm then checks all
#     less abundant barcode if they have only (0,1) edit distance and if they
#     map to the same position on the genome. If so, then these barcodes are merged
#     If not then we continue
#
#     Finally, barcode abundances are written to the output file.
#     :param barcode_2_count:
#     :param barcode_2_potential_position:
#     :param out_barcode_file:
#     :return:
#     '''
#
#     barcodes: List[str] = list(map(lambda x: x[0], barcode_2_count.most_common()))
#     initial_barcode_count = len(barcodes)
#     final_barcodes: List[str] = []
#     if barcode_2_potential_position == 0:
#         final_barcodes = barcodes
#     else:
#         while True:
#             top_barcode: str = barcodes.pop(0)
#             final_barcodes.append(top_barcode)
#             if len(barcodes) == 0:
#                 break
#             barcodes_to_group = set()
#             for barcode in barcodes:
#                 if is_similar_barcode(top_barcode, barcode_2_potential_position[top_barcode][0], barcode, barcode_2_potential_position[barcode][0], barcode_max_edit_distance):
#                     barcodes_to_group.add(barcode)
#             barcodes_2_keep = []
#             for barcode in barcodes:
#                 if barcode not in barcodes_to_group:
#                     barcodes_2_keep.append(barcode)
#             barcodes = barcodes_2_keep
#         final_barcodes = set(final_barcodes)
#     final_barcode_count = len(final_barcodes)
#     logging.info(f'Initial number of barcodes:\t{initial_barcode_count}')
#     logging.info(f'Corrected number of barcodes:\t{final_barcode_count}')
#
#     with open(out_barcode_file, 'w') as handle:
#         for barcode, potpos in barcode_2_potential_position.items():
#             if barcode in final_barcodes:
#                 tmp = '\t'.join(map(lambda x: str(x), potpos[0]))
#                 handle.write(f'{barcode}\t{barcode_2_count[barcode]}\t{tmp}\t{potpos[1]}\n')


# def map_host(temp_fasta_file: str, temp_blastn_file: str, blastdb: str, out_barcode_file: str, blast_threads: int = 1, barcode_max_edit_distance: int=1):
#     logging.info('----------------')
#     logging.info('Step 2.1: Start blastn alignment of dereplicated barcode/host pairs')
#
#     map_host_map(temp_fasta_file, temp_blastn_file, blastdb, blast_threads)
#
#     logging.info('Step 2.1: Finished blastn alignment of dereplicated barcode/host pairs')
#     logging.info('----------------')
#     logging.info('Step 2.2: Start annotating barcodes')
#
#     header_2_sequence, query_2_besthits = map_host_load(temp_fasta_file, temp_blastn_file)
#     logging.info(
#         f'Host sequences with alignments:\t{len(query_2_besthits)}/{len(header_2_sequence)}\t{str(100.0 * len(query_2_besthits) / len(header_2_sequence))[:5]}%')
#
#     barcode_2_potential_position, barcode_2_count = map_annotate(header_2_sequence, query_2_besthits)
#
#     logging.info('Step 2.2: Finished annotating barcodes')
#     logging.info('----------------')
#     logging.info('Step 2.2: Start merging barcodes by edit distance and host position')
#
#     map_host_barcode_merge(barcode_2_count, barcode_2_potential_position, out_barcode_file, barcode_max_edit_distance)
#
#     logging.info('Step 2.2: Finished merging barcodes by edit distance and host position')
#

def demultiplex_tnseq():
    parser = argparse.ArgumentParser(description = 'Demultiplex samples according to sequencing barcodes.\n\n' \
                                                   'A file mapping barcodes to sample names should be provided.\n' \
                                                   'It should consist of two tab-separated columns with no header.\n' \
                                                   'E.g.:\n' \
                                                   'AGTT    AK387\n' \
                                                   'CCTT    AK388', formatter_class = RawTextHelpFormatter)
    parser.add_argument('-bc', action='store', required=True, help='Barcode to sample name mapping file (required)')
    parser.add_argument('-r1', action='store', required=True, help='R1 FastA/Q file. Gzip input allowed (required)')
    parser.add_argument('-o', action='store', required=False, help='Output directory')

    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        parser.print_help()
        shutdown(1)

    if args.o is None:
        args.o = ""

    index_2_sample = {}
    with open(args.bc) as handle:
        for line in handle:
            line = line.strip()
            index,samplename = line.split('\t')
            index_2_sample[index] = samplename

    k2 = 'GTGTATAAGAGACAG'

    index_2_writer = {}
    for index, samplename in index_2_sample.items():
        index_2_writer[index] = open(args.o + '/' + samplename + '.fasta', 'w')

    k2_found = 0
    good_index = 0

    for total_inserts, r1 in enumerate(tnseq2.src.sequence.stream_fa(args.r1), 1):
        if k2 in r1.sequence:
            k2_found += 1
            startpos_k2 = r1.sequence.index(k2)
            if startpos_k2 >= 30:
                index = r1.sequence[startpos_k2+15:startpos_k2+23]
                if index in index_2_sample:
                    good_index += 1
                    index_2_writer[index].write(f'>{r1.header}\n{r1.sequence}\n')

        if total_inserts % 100000 == 0:
            print(f'Found/Good/Total\t{k2_found}/{good_index}/{total_inserts}')

    for index, writer in index_2_writer.items():
        writer.close()
 


def map_tnseq():
    parser = argparse.ArgumentParser(description = 'Extract host sequences, align to host genome and map to barcode sequences.\n' \
                                                  'The current tool allows for tn5 transposons and identifies the following structure\n\n' \
                                                  '|r1=(0-k)bp|BARCODE=17bp|tp1=13bp|tp2=15bp|r2=8bp|HOST=(15-j)bp||ADAPTER=(0-n)bp|\n' \
                                                  '|----------|============|--------|========|------|=============||---------------|\n\n' \
                                                  'The tool searches for the (tp2=GTGTATAAGAGACAG) and finds the barcode\n' \
                                                  'and the host sequence. Then we align the host sequence against the\n' \
                                                  'host genome and map the position to the barcode. Finally, we perform\n'
                                                   'barcode correction based on abundance and edit distance.' \
                                                  'Write hansr@ethz.ch for any requests', formatter_class = RawTextHelpFormatter)


    parser.add_argument('-r1', action='store', required=True, help='R1 FastA/Q file. Gzip input allowed. (required)')
    parser.add_argument('-r2', action='store', required=True, help='R2 FastA/Q file. Gzip input allowed. (required)')
    parser.add_argument('-o',  action='store', required=True, help='Output barcode file. (Input for quantify step). (required)')
    #parser.add_argument('-s',  action='store', required=True, help='Output file with general summarizing statistics. (required)')
    parser.add_argument('-db', action='store', required=True, help='Blast database name of the host genome. (required)')
    parser.add_argument('-t', action='store', help='Temp folder. Use folder of output file as default.')
    parser.add_argument('-c', action='store', help='Max edit distance for barcode correction. (default=0)', type=int, default=0, choices=[0, 1])
    parser.add_argument('-n', action='store', help='Number of threads for blast. (default=1)', type=int,
                        default=1)
    parser.add_argument('-f', action='store_true', default=False, help='Force overwrite existing files.')

    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        parser.print_help()
        shutdown(1)
    startup()
    # Input
    in_r1_file: str = args.r1
    in_r2_file: str = args.r2
    blastdb: str = args.db
    # Output
    out_barcode_file: str = args.o
    temp_folder: str = args.t
    # Params
    force: bool = args.f
    barcode_correction: int = args.c
    mode: str = 'tn5'
    blastdb_threads: int = args.n


    if pathlib.Path(out_barcode_file).exists() and not force:
        logging.error(f'File {out_barcode_file} exists. Aborting!. Use -f flag to overwrite existing files.')
        shutdown(1)

    if not pathlib.Path(in_r1_file).exists():
        logging.error(f'File {in_r1_file} does not exist. Aborting!')
        shutdown(1)
    if not pathlib.Path(in_r2_file).exists():
        logging.error(f'File {in_r2_file} does not exist. Aborting!')
        shutdown(1)
    if temp_folder:
        if not pathlib.Path(temp_folder).exists():
            pathlib.Path(temp_folder).mkdir(parents=False, exist_ok=True)
    else:
        temp_folder = str(pathlib.Path(out_barcode_file).resolve().parent)
    temp_fasta_file = pathlib.Path(temp_folder).joinpath('tnseq_host_reads.fasta')
    if temp_fasta_file.exists() and not force:
        logging.error(f'File {str(temp_fasta_file)} exists. Aborting!. Use -f flag to overwrite existing files.')
        shutdown(1)

    temp_blastn_file = pathlib.Path(temp_folder).joinpath('tnseq_host_reads.blastn')
    if temp_blastn_file.exists() and not force:
        logging.error(f'File {str(temp_blastn_file)} exists. Aborting!. Use -f flag to overwrite existing files.')
        shutdown(1)

    logging.info('Input:')
    logging.info(f'\tR1 file:\t{in_r1_file}')
    logging.info(f'\tR2 file:\t{in_r2_file}')
    logging.info(f'\tBlastDB:\t{blastdb}')
    logging.info('Output:')
    logging.info(f'\tBarcode file:\t{out_barcode_file}')
    logging.info(f'\tTemp folder:\t{temp_folder}')
    logging.info(f'\tTemp fasta file:\t{temp_fasta_file}')
    logging.info(f'\tTemp blastn file:\t{temp_blastn_file}')
    logging.info('Parameters:')
    logging.info(f'\tMax barcode distance:\t{barcode_correction}')
    logging.info(f'\tBlast threads:\t{blastdb_threads}')

    logging.info('----------------------------------------------------------')
    logging.info('Step 1: Start extraction of barcodes/host sequences')
    tnseq2.src.mapping.prepare_extract_barcodes(in_r1_file, in_r2_file, temp_fasta_file)
    #prepare_extract_barcodes(in_r1_file, in_r2_file, temp_fasta_file, mode)
    logging.info('Step 1: Finished extraction of barcodes/host sequences')
    logging.info('----------------------------------------------------------')
    logging.info('Step 2: Start alignment/assignment of host sequences')
    tnseq2.src.mapping.map_host(str(temp_fasta_file), str(temp_blastn_file), blastdb, out_barcode_file, blast_threads=blastdb_threads, barcode_max_edit_distance=barcode_correction)
    logging.info('Step 2: Finished alignment/assignment of host sequences')
    logging.info('----------------------------------------------------------')
    shutdown(0)




def quantify_tnseq():


    parser = argparse.ArgumentParser(description='Extract barcodes from sequence file, map them against the provided barcodes\n' \
                                                 'and quantify them using barcode counts.' \
                                                 'The current tool allows for tn5 transposons and identifies the following structure\n\n' \
                                                 '|r1=(0-k)bp|BARCODE=17bp|tp1=13bp|tp2=15bp|rest=(0-j)bp|\n' \
                                                 '|----------|============|--------|========|------------|\n\n' \
                                                 'The tool searches for the (tp2=GTGTATAAGAGACAG) and finds the barcode\n' \
                                                 'based on the distance to the tp2 sequence. This barcode sequence is matched\n' \
                                                 'against known barcodes allowing for n (0-2) edit distances and quantified.\n' \
                                                 'Write hansr@ethz.ch for any requests', formatter_class=RawTextHelpFormatter)

    parser.add_argument('-r', action='store', required=True, help='FastA/Q file. Gzip input allowed. (required)')
    parser.add_argument('-o', action='store', required=True, help='Output barcode counts file. (required)')
    parser.add_argument('-b', action='store', required=True,
                        help='Known barcode file (Output from map step). (required)')
    parser.add_argument('-c', action='store', help='Max edit distance for barcode correction. (default=2)', type=int,
                        default=2, choices=[0, 1, 2])
    parser.add_argument('-f', action='store_true', default=False, help='Force overwrite existing files.')
    parser.add_argument('-s', action='store', required=False, help='Number of bases stuffing added to the 13bp spacer', type=int, default=0)

    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        parser.print_help()
        shutdown(1)
    startup()
    # Input
    in_fq_file: str = args.r
    in_bc_file: str = args.b
    # Output
    out_barcode_count_file: str = args.o

    # Params
    force: bool = args.f
    max_edit_distance: int = args.c
    mode: str = 'tn5'

    if pathlib.Path(out_barcode_count_file).exists() and not force:
        logging.error(f'File {out_barcode_count_file} exists. Aborting!. Use -f flag to overwrite existing files.')
        shutdown(1)
    if not pathlib.Path(in_bc_file).exists():
        logging.error(f'File {in_bc_file} does not exist. Aborting!')
        shutdown(1)

    if not pathlib.Path(in_fq_file).exists():
        logging.error(f'File {in_fq_file} does not exist. Aborting!')
        shutdown(1)

    logging.info('Input:')
    logging.info(f'\tReads file:\t{in_fq_file}')
    logging.info(f'\tBarcode map file:\t{in_bc_file}')
    logging.info('Output:')
    logging.info(f'\tBarcode count file:\t{out_barcode_count_file}')
    logging.info('Parameters:')
    logging.info(f'\tMax barcode distance:\t{max_edit_distance}')

    logging.info('----------------------------------------------------------')
    logging.info('Step 1: Reading barcode map file')
    barcode_2_position, barcode_2_abundance = quantify_read_barcode_map_files(in_bc_file)
    logging.info(f'\tFound {len(barcode_2_position)} annotated barcodes')
    logging.info('Step 1: Finished Reading barcode map file.')

    logging.info('----------------------------------------------------------')
    logging.info('Step 2: extracting/annotating/correcting sample specific barcodes')
    quantify_extract_annotated_correct(barcode_2_position, barcode_2_abundance, in_fq_file, out_barcode_count_file, max_edit_distance=max_edit_distance, mode=mode, pad=args.s)
    logging.info('Step 2: Finished extracting/annotating/correcting sample specific barcodes')
    logging.info('----------------------------------------------------------')
    shutdown(0)


def startup():
    '''
    Generic Startup commands

    '''
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)
    logging.info('Starting tnseq pipeline')


def shutdown(status=0):
    '''
    Generic shutdown commands
    '''
    logging.info('Finishing tnseq pipeline with status:\t{}'.format(status))
    sys.exit(status)

def main():

    parser = argparse.ArgumentParser(description='A toolkit to map and quantify Tn-Seq data',
                                     usage='''tnseq <command> [<args>]

    Command options
        map         Map barcodes to positions in the host genome. (First step)
        quantify    Quantify barcode abundances. (Second Step)
        demultiplex Demultiplex samples according to sequencing barcodes. (Third step)
    ''')

    parser.add_argument('command', help='Subcommand to run: map|quantify|demultiplex')

    args = parser.parse_args(sys.argv[1:2])
    if args.command == 'map':
        #map -r1 /Volumes/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_map/RB-TnSeq_mapping_library_1.1.fq.gz -r2 /Volumes/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_map/RB-TnSeq_mapping_library_1.2.fq.gz -o /Volumes/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_quantify/barcodes.txt -db /Volumes/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/blastdb/Salmonella_genome_FQ312003.1_SL1344.fasta -c 1 -n 4 -f
        map_tnseq()
    elif args.command == 'quantify':
        #quantify -r /Volumes/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_quantify/1315-19-library11_1-TV3379-u443-d1-feces.fasta.gz -o /Volumes/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_quantify/1315-19-library11_1-TV3379-u443-d1-feces.count -b /Volumes/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_quantify/barcodes.txt -c 1 -f
        quantify_tnseq()
    elif args.command == 'demultiplex':
        demultiplex_tnseq()
    else:
        parser.print_help()
        shutdown(0)
    shutdown(0)

if __name__ == '__main__':
    main()





