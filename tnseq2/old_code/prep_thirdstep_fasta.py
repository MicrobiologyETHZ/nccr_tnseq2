import sys
import tnseq.sequence
from typing import Generator, Tuple, List

# Multiplexing barcode file format has two tab-separated columns, no headers:
# <BARCODE>   <SAMPLE>
# e.g.:
# AGTT  AK387
# CCTT  AK388

index_2_sample = {}
with open(sys.argv[2]) as handle:
    for line in handle:
        line = line.strip()
        index,samplename = line.split('\t')
        index = index + index
        index_2_sample[index] = samplename

k2 = 'GTGTATAAGAGACAG'

index_2_writer = {}
for index, samplename in index_2_sample.items():
    index_2_writer[index] = open(sys.argv[3] + '/' + samplename + '.fasta', 'w')

k2_in_r1 = 0
k2_in_r1_but_barcode_short = 0

good_index = 0
bad_index = 0

for total_inserts, r1 in enumerate(tnseq.sequence.stream_fa(sys.argv[1]), 1):
    if k2 in r1.sequence:
        k2_in_r1 += 1
        startpos_k2 = r1.sequence.index(k2)
        if startpos_k2 != 30:
            k2_in_r1_but_barcode_short += 1
            continue
        barcode = r1.sequence[:17]
        index = r1.sequence[startpos_k2+15: startpos_k2+23]
        if index in index_2_sample:
            good_index += 1
            index_2_writer[index].write(f'>{r1.header}\n{r1.sequence}\n')
        else:
            bad_index += 1

    if (good_index + bad_index) % 100000 == 0:
        print(f'Total\t{total_inserts}')
        print(f'Good Index \t{good_index}')
        print(f'Bad Index:\t{bad_index}')
        print('---------')

for index, writer in index_2_writer.items():
    writer.close()

