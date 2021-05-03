import pandas as pd
from pathlib import Path
import sys


def read_tnseq_count_file(f):
    df = pd.read_csv(f, index_col=0)
    df['sampleName'] = f.stem.split("_counts")[0]
    return df

# def get_metafile(meta_dir):
#     meta_file = (pd.read_table(Path(meta_dir) / f"{dnaid}_metadata.txt", header=None,
#                                names=['sampleID', 'library', 'experiment', 'DN1', 'mouse', 'day', 'organ'],
#                                dtype={'sampleID': str})
#                  .drop('DN1', axis=1))
#     return meta_file


def merge_counts(count_dir, meta_file='', dnaid='', k='_mapped'):
    # concat all the count files
    counts = [read_tnseq_count_file(f) for f in (Path(count_dir)).iterdir() if k in f.name]
    df = pd.concat(counts)
    if meta_file:
        # map columns using the meta file
        meta = pd.read_table(meta_file)
        if 'sampleName' not in meta.columns:
            print('No sampleName in metadata found')
            sys.exit(1)
        df = df.merge(meta, on='sampleName')
        df['sampleID'] = df['mouse'] + "_" + df['day']
    if dnaid and dnaid not in df.columns:
        df['dnaid'] = dnaid
    return df


def merge_controls(count_dir, control_file, meta_file,  dnaid, k='unmapped'):
    controls = pd.read_table(control_file, header=None, names=['DN', 'barcode', 'phenotype', 'conc'])
    counts = merge_counts(count_dir, meta_file, dnaid, k=k)
    return controls.merge(counts, how='left', on='barcode')


def final_merge(count_dir,  meta_file='', control_file='', dnaid='', out_file='',  merged_control='', k='_mapped', ck='unmapped'):
    if not out_file:
        out_file = Path(count_dir)/"merged_counts.csv"
    df = merge_counts(count_dir, meta_file=meta_file, dnaid=dnaid, k=k)
    df.to_csv(out_file)
    if control_file:
        if not merged_control:
            merged_control = Path(count_dir) / 'merged_controls.csv'
        control_df = merge_controls(count_dir, control_file, meta_file,  dnaid, k=ck)
        control_df.to_csv(merged_control)


if __name__ == "__main__":
    count_dir = sys.argv[1]
    meta_file = sys.argv[2]
    dnaid = sys.argv[3]
    #out_file = sys.argv[4]
    control_file = sys.argv[4]
    merge_counts(count_dir, meta_file, dnaid)
    merge_controls(count_dir, control_file, meta_file, dnaid)



