
def map_host(temp_fasta_file: str, temp_blastn_file: str, blastdb: str, out_barcode_file: str, blast_threads: int = 1, barcode_max_edit_distance: int=1):
    logging.info('----------------')
    logging.info('Step 2.1: Start blastn alignment of dereplicated barcode/host pairs')

    map_host_map(temp_fasta_file, temp_blastn_file, blastdb, blast_threads)

    logging.info('Step 2.1: Finished blastn alignment of dereplicated barcode/host pairs')
    logging.info('----------------')
    logging.info('Step 2.2: Start annotating barcodes')

    header_2_sequence, query_2_besthits = map_host_load(temp_fasta_file, temp_blastn_file)
    logging.info(
        f'Host sequences with alignments:\t{len(query_2_besthits)}/{len(header_2_sequence)}\t{str(100.0 * len(query_2_besthits) / len(header_2_sequence))[:5]}%')

    barcode_2_potential_position, barcode_2_count = map_annotate(header_2_sequence, query_2_besthits)

    logging.info('Step 2.2: Finished annotating barcodes')
    logging.info('----------------')
    logging.info('Step 2.2: Start merging barcodes by edit distance and host position')

    map_host_barcode_merge(barcode_2_count, barcode_2_potential_position, out_barcode_file, barcode_max_edit_distance)

    logging.info('Step 2.2: Finished merging barcodes by edit distance and host position')





def parse_mapping(blastfile: str) -> Dict[str, List[List[str]]]:
    '''
    Load blastn file and extract for each query the best hit
    by bitscore. Will return a list with one or multiple best
    hits.
    :param blastfile:
    :return:
    '''
    query_2_hits: DefaultDict[str, List[List[str]]] = collections.defaultdict(list)
    with open(blastfile) as handle:
        for line in handle:
            splits = line.strip().split('\t')
            query_2_hits[splits[0]].append(splits)
    query_2_besthits: Dict[str, List[List[str]]] = {}
    for query, hits in query_2_hits.items():
        bestscore: float = max(map(lambda x: float(x[9]), hits))
        best_score_hits: List[List[str]] = list(map(lambda y: [y[1], y[6], y[11]], filter(lambda x: float(x[9]) >= bestscore, hits)))

        query_2_besthits[query] = best_score_hits
    return query_2_besthits



def map_host_load(temp_fasta_file: str, temp_blastn_file: str) -> Tuple[Dict[str, str], Dict[str, List[List[str]]]]:
    '''
    Helper method to load blast and fasta file
    :param temp_fasta_file:
    :param temp_blastn_file:
    :return:
    '''
    fq1_stream = tnseq2.sequence.stream_fa(temp_fasta_file)
    header_2_sequence: Dict[str, str] = {}
    for read in fq1_stream:
        header_2_sequence[read.header] = read.sequence
    query_2_besthits: Dict[str, List[List[str]]] = parse_mapping(temp_blastn_file)
    return header_2_sequence, query_2_besthits



def map_annotate(header_2_sequence : Dict[str, str], query_2_besthits: Dict[str, List[List[str]]]) -> Tuple[Dict[str, Tuple[Tuple[str, str, str], int]], Counter[str]]:
    '''
    Find the best position/orientation for each barcode.
    Each barcode can consist of multiple host sequences.
    The correct host sequences can be identified the following:
    - They appear in large quantity
    - Even if they dont appear in large quantity they map against
        the sample position in the host genome

    For each barcode we look into the alignments for each associated
    host sequence and pick the most likely position. The most likely
    position is the one that appears most often.

    :param header_2_sequence:
    :param query_2_besthits:
    :return:
    '''
    barcode_2_alignments: Dict[str, List[Tuple[int, List[str]]]] = {}
    barcode_2_count: Counter[str] = collections.Counter()
    barcode_2_potential_position: Dict[str, Tuple[Tuple[str, str, str], int]] = {}
    for header in header_2_sequence:
        splits: List[str] = header.split('_')
        barcode: str = splits[2]
        cnt: int = int(splits[4])
        barcode_2_count[barcode] += cnt
        barcode_2_alignments[barcode] = []

    if len(query_2_besthits) == 0:
        logging.error('There were 0 alignments found in the blast output found. Most likely a problem with the input file, the blast installation and/or blast database.')
        shutdown(1)
    for query, besthits in query_2_besthits.items():
        splits: List[str] = query.split('_')
        barcode: str = splits[2]
        cnt: int = int(splits[4])
        barcode_2_alignments[barcode].append((cnt, besthits))

    # for each barcode i get all alignments
    # Then group the alignments by start position (beg plus, end minus)
    # (start, ref, orientation) = count
    for barcode, total_count in barcode_2_count.most_common():
        alignments: List[Tuple[int, List[str]]] = barcode_2_alignments[barcode]
        if len(alignments) == 0:
            barcode_2_potential_position[barcode] = (('-1', 'None', 'plus'), total_count)
            continue
        potential_positions: Counter[Tuple[str, str, str]] = collections.Counter()
        for alignment in alignments:
            alignment_count: str = alignment[0]
            seq_alignments: List[List[str]] = alignment[1]
            weight = float(alignment_count) / len(seq_alignments)
            for seqaln in seq_alignments:
                orientation: str = seqaln[2]
                reference: str = seqaln[0]
                startpos: int = int(seqaln[1])
                potential_positions[(startpos, reference, orientation)] += weight
        mc: Tuple[Tuple[str, str, str], int] = potential_positions.most_common()[0]
        barcode_2_potential_position[barcode] = mc

    return barcode_2_potential_position, barcode_2_count

def map_host_barcode_merge(barcode_2_count, barcode_2_potential_position, out_barcode_file, barcode_max_edit_distance: int=0):
    '''
    Gets a list of barcodes/abundance and their respective postion in the
    host genome. Barcodes are then sorted by abundance as more abundant
    barcode are potentially more credible. The algorithm then checks all
    less abundant barcode if they have only (0,1) edit distance and if they
    map to the same position on the genome. If so, then these barcodes are merged
    If not then we continue

    Finally, barcode abundances are written to the output file.
    :param barcode_2_count:
    :param barcode_2_potential_position:
    :param out_barcode_file:
    :return:
    '''

    barcodes: List[str] = list(map(lambda x: x[0], barcode_2_count.most_common()))
    initial_barcode_count = len(barcodes)
    final_barcodes: List[str] = []
    if barcode_2_potential_position == 0:
        final_barcodes = barcodes
    else:
        while True:
            top_barcode: str = barcodes.pop(0)
            final_barcodes.append(top_barcode)
            if len(barcodes) == 0:
                break
            barcodes_to_group = set()
            for barcode in barcodes:
                if is_similar_barcode(top_barcode, barcode_2_potential_position[top_barcode][0], barcode, barcode_2_potential_position[barcode][0], barcode_max_edit_distance):
                    barcodes_to_group.add(barcode)
            barcodes_2_keep = []
            for barcode in barcodes:
                if barcode not in barcodes_to_group:
                    barcodes_2_keep.append(barcode)
            barcodes = barcodes_2_keep
        final_barcodes = set(final_barcodes)
    final_barcode_count = len(final_barcodes)
    logging.info(f'Initial number of barcodes:\t{initial_barcode_count}')
    logging.info(f'Corrected number of barcodes:\t{final_barcode_count}')

    with open(out_barcode_file, 'w') as handle:
        for barcode, potpos in barcode_2_potential_position.items():
            if barcode in final_barcodes:
                tmp = '\t'.join(map(lambda x: str(x), potpos[0]))
                handle.write(f'{barcode}\t{barcode_2_count[barcode]}\t{tmp}\t{potpos[1]}\n')


def is_similar_barcode(barcode1: str, pos1: int, barcode2: str, pos2: int, barcode_max_edit_distance: int = 0) -> bool:
    '''
    Takes 2 barcodes and their positions in the host genomes as input
    and checks if the barcodes have a edit distance smaller then
    {barcode_max_edit_distance}. If so then they should map to
    roughly the same position (+/- 5 bases) in the reference
    genome. Then they are considered similar. Otherwise the
    method returns false.
    :param barcode1:
    :param pos1:
    :param barcode2:
    :param pos2:
    :param barcode_max_edit_distance:
    :return:
    '''
    if editdistance(barcode1, barcode2) <= barcode_max_edit_distance:
        if pos1[1] == pos2[1] and pos1[2] == pos2[2]:
            if (pos1[0] - 5) <= pos2[0] <= (pos1[0] + 5):
                return True
    return False
