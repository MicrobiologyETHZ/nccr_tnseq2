from tnseq2.src.sequence import stream_fa, FastA
from Bio.Seq import Seq
# def demultiplex_tnseq():
#     parser = argparse.ArgumentParser(description='Demultiplex samples according to sequencing barcodes.\n\n' \
#                                                  'A file mapping barcodes to sample names should be provided.\n' \
#                                                  'It should consist of two tab-separated columns with no header.\n' \
#                                                  'E.g.:\n' \
#                                                  'AGTT    AK387\n' \
#                                                  'CCTT    AK388', formatter_class=RawTextHelpFormatter)
#     parser.add_argument('-bc', action='store', required=True, help='Barcode to sample name mapping file (required)')
#     parser.add_argument('-r1', action='store', required=True, help='R1 FastA/Q file. Gzip input allowed (required)')
#     parser.add_argument('-o', action='store', required=False, help='Output directory')
#
#     try:
#         args = parser.parse_args(sys.argv[2:])
#     except:
#         parser.print_help()
#         shutdown(1)
#
#     if args.o is None:
#         args.o = ""


def demux_tnseq(r1_fastq, demux_bc_file, out_dir='.', name='',  tn='GTGTATAAGAGACAG:17:13:before', rc=True):
    index_2_sample = {}
    with open(demux_bc_file) as handle:
        for line in handle:
            line = line.strip()
            index, samplename = line.split('\t')
            if name:
                samplename = name + "_" + samplename
            if rc:
                index = str(Seq(index).reverse_complement())
            index_2_sample["{}{}".format(index, index)] = samplename
    k2 = tn.split(':')[0]
    bcLen = int(tn.split(':')[1])
    bc2tp2 = int(tn.split(':')[2])
    index_2_writer = {}
    for index, samplename in index_2_sample.items():
        index_2_writer[index] = open(out_dir + '/' + samplename + '.fasta', 'w')
    k2_found = 0
    good_index = 0
    for total_inserts, r1 in enumerate(stream_fa(r1_fastq), 1):
        if k2 in r1.sequence:
            k2_found += 1
            startpos_k2 = r1.sequence.index(k2)
            if startpos_k2 >= bcLen + bc2tp2:
                index = r1.sequence[startpos_k2 + len(k2):startpos_k2 + len(k2)+8]
                if index in index_2_sample:
                    good_index += 1
                    index_2_writer[index].write(f'>{r1.header}\n{r1.sequence}\n')

        if total_inserts % 100000 == 0:
            print(f'Found/Good/Total\t{k2_found}/{good_index}/{total_inserts}')

    for index, writer in index_2_writer.items():
        writer.close()

