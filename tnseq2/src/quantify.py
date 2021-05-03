from tnseq2.src.sequence import stream_fa, FastA
from tnseq2.src.extract_barcodes import extract_barcode_host
from tnseq2.src.commandline import *
from typing import List, Tuple, Generator, Dict, Iterator, DefaultDict, Counter, Set
import logging
import collections
import sys
import pandas as pd
from pathlib import Path

def quantify_load_fq_barcodes(in_fq_file: str, tp2: str='GTGTATAAGAGACAG',
                              bc2tp2:int=13, bcLen:int=17, before:bool=True) -> Counter[str]:
    '''
    Count barcodes in a sequencing file (fasta/fastq)
    This currently the default:

    -----|BARCODE|----------|TN end sequence (tp2)|---Host------
    -----|-bcLen-|--bc2tp2--|---------tp2---------|-------------
    -----|-17bp--|---13bp---|---------15bp--------|----?--------
    ---(-30)---(-13)-------(0)---------------------------------


    :param in_fq_file:
    :param tp2:
    :return:
    '''
    fq1_stream = stream_fa(in_fq_file)
    # Counting reads with and without transposon sequence
    with_tp2: int = 0
    without_tp2: int = 0
    with_tp_but_short: int = 0
    cnter: Counter[str] = collections.Counter()
    for total_inserts, r1 in enumerate(fq1_stream, 1):

        if total_inserts % 1000000 == 0:
            logging.info(f'\tReads processed:\t{total_inserts}')
        if tp2 in r1.sequence:
            barcode, _ = extract_barcode_host(r1, tp2, bc2tp2, bcLen, before)
            if not barcode:
                with_tp_but_short += 1
            else:
                with_tp2 += 1
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


def editdistance(seq1: str, seq2: str) -> int:
    '''
    Calculate the edit distance between 2 sequences with identical length.
    Will throw an error if the length of both sequences differs
    :param seq1:
    :param seq2:
    :return:
    '''

    if seq1 == seq2:
        return 0
    if len(seq1) != len(seq2):
        raise Exception(
            f'{seq1} and {seq2} have different length. Edit distance can be computed on same length sequences only.')
    dist = 0
    for letter1, letter2 in zip(seq1, seq2):
        if letter1 != letter2:
            dist += 1
    return dist


def get_similar(bc_to_id, bc_mapped):
    # slow slow slow
    distances = {}
    for total, bc in enumerate(bc_to_id):
        min_dist = 1000
        match = ''
        if total % 1000 == 0:
            print(f'\tBarcodes processed:\t{total}')
        for mbc in bc_mapped:
            actual_dist = editdistance(bc, mbc)
            if actual_dist < min_dist:
                min_dist = actual_dist
                match = mbc
        distances[bc] = [min_dist, match]

    dist = pd.DataFrame(distances).T.reset_index()
    print(dist.head())
    dist.columns = ['barcode', 'editdistance', 'match']

    return dist


def annotate_barcodes(cnter, barcode_map_file, edit_cutoff=3):
    cnts_df = pd.DataFrame.from_dict(cnter, orient='index').reset_index()
    cnts_df.columns = ['barcode', 'cnt']
    cnts_df = cnts_df[cnts_df['cnt'] > 1]
    if cnts_df.empty:
        logging.error('No barcodes with counts > 1 found')
        sys.exit(1)

    bc_df = pd.read_csv(barcode_map_file)
    bc_df.columns = 'barcode,libcnt,sstart,send,sseqid,sstrand,multimap,ShortName,locus_tag'.split(',')

    annotated_cnts = cnts_df.merge(bc_df, how='left', on='barcode')

    with_ids = annotated_cnts[annotated_cnts.sstart.notnull()]
    print(with_ids.shape)
    no_ids = annotated_cnts[annotated_cnts.sstart.isna()]
    print(no_ids.shape)
    distances = get_similar(no_ids.barcode.values, bc_df.barcode.values)
    no_matches = distances[distances.editdistance >= edit_cutoff]
    with_matches = distances[distances.editdistance < edit_cutoff]
    with_matches = (with_matches.merge(cnts_df, how='left', on='barcode')
                    .drop(['barcode', 'editdistance'], axis=1)
                    .rename({'match': 'barcode'}, axis=1))

    never_ided = no_ids[no_ids.barcode.isin(no_matches.barcode.values)]
    all_ids = pd.concat([with_ids[['barcode', 'cnt']].drop_duplicates(), with_matches])
    all_ids = all_ids.groupby('barcode').cnt.sum().reset_index()
    final_ids = all_ids.merge(
        annotated_cnts[['barcode', 'libcnt', 'sstart', 'sseqid', 'sstrand', 'multimap', 'ShortName', 'locus_tag']],
        on='barcode', how='left')
    final_ids.ShortName.fillna(final_ids.barcode, inplace=True)
    return final_ids, never_ided.dropna(axis=1)


def quantify(fasta_file, map_file, outdir, prefix):
    cnter = quantify_load_fq_barcodes(fasta_file)
    cnts, no_ids = annotate_barcodes(cnter, map_file, edit_cutoff=3)
    cnts.to_csv(Path(outdir) / f'{prefix}_counts_mapped.csv')
    no_ids.to_csv(Path(outdir) / f'{prefix}_counts_unmapped.csv')
    return cnts, no_ids


if __name__ == "__main__":
    fasta_file = sys.argv[1]
    map_file = sys.argv[2]
    outdir = sys.argv[3]
    prefix = sys.argv[4]
    quantify(fasta_file, map_file, outdir, prefix)



