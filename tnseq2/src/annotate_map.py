from pathlib import Path
import sys
import pandas as pd


def add_gene_annotations(bedfile_path, barcode_file):
    bedfile = pd.read_table(bedfile_path, header=None)
    bedfile.columns = ['seq', 'position', 'gene_info', 'barcode']
    bedfile['ShortName'] = bedfile['gene_info'].str.extract(r'(Name=.+?;)', expand=False).str.strip(
        'Name=').str.strip(';')
    bedfile['locus_tag'] = bedfile['gene_info'].str.extract(r'(locus_tag=.+?$)', expand=False).str.strip(
        'locus_tag=').str.strip(';')
    potential_positions = pd.read_csv(barcode_file)
    final_map = Path(barcode_file).with_suffix('.annotated.csv')
    fdf = potential_positions.merge(bedfile[['barcode', 'ShortName', 'locus_tag']], how='left', on='barcode')
    fdf.to_csv(final_map, index=False)


if __name__ == "__main__":
    bedfile = sys.argv[1]
    barcode_file = sys.argv[2]
    add_gene_annotations(bedfile, barcode_file)

