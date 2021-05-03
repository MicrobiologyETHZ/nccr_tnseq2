import pandas as pd
import plotnine as p9
from pathlib import Path
from scipy import stats
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math
import plotnine as p9

import skmisc
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA
from diffexpr.py_deseq import py_DESeq2


from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')

from scipy.stats import norm
import statsmodels
import scipy

from datetime import date


today = date.today().strftime("%Y_%m_%d")


def load_files(samples, results_dir):
    df_list = []
    for sample in samples:
        df_list.append(pd.read_csv(Path(results_dir)/sample/f"{sample}_merged_counts.csv", index_col = 0))
    return pd.concat(df_list)


def subset_experiment(df,  dnaid,  exp, col1='experiment', col2='dnaid'):
    '''
    example query string : '(exp=="TV5490A") & (dnaid == "dnaid2023")'
    '''
    query_string = f'({col1} == "{exp}") & ({col2} == "{dnaid}")'
    return df.copy().query(query_string)


def calculate_correlation(exp_df, controls, for_each='sampleID', how='log', cutoff=0.9, phenotype='wt'):
    """
    Subset counts for control barcodes
    Calculate correlation on log counts (log), log counts, but keep 0 (log_w_0), or raw data (raw)

    """
    control_cnts = (controls.merge(exp_df, left_on='barcode', right_on='barcode')
                    .drop(['DN', 'position', 'seq', 'strand', 'locus', 'gene'], axis=1))

    if how == 'raw':
        col1 = 'conc'
        col2 = 'cnt'
    else:
        control_cnts['logConc'] = np.log10(control_cnts['conc'])
        col1 = 'logConc'
        if how == 'log':
            control_cnts['logCnts'] = np.log10(control_cnts['cnt'])
        elif how == 'log_w_0':
            control_cnts['logCnts'] = np.log10(control_cnts['cnt'].replace({0: 1}))
        col2 = 'logCnts'

    corr_df = control_cnts.groupby(['phenotype', for_each])[[col1, col2]].corr()
    corr_df = corr_df.reset_index()
    corr_df = corr_df[corr_df['level_2'] == col1].drop(['level_2', col1], axis=1)
    corr_df.columns = ['phenotype', 'sampleID', 'R']

    good_samples = corr_df[(corr_df.R > cutoff) & (corr_df.phenotype == phenotype)].sampleID.values
    return corr_df, good_samples


def filter_inoculum(exp_df, filter_below=0):
    filt_df = exp_df.copy().pivot(index='barcode', columns='sampleID', values='cnt')
    filt_df = filt_df.fillna(0)
    filt_df = filt_df[(filt_df['inoculum_d0'] >=filter_below)& (filt_df['unenriched_inoculum_d0'] >= filter_below)]
    return filt_df


def filter_samples(exp_df, good_samples):
    return exp_df.copy()[exp_df.sampleID.isin(good_samples)]


def generate_DE_dataset(exp_df, good_samples, filter_below = 1):
    sample_data = exp_df[['sampleID', 'mouse', 'day', 'organ', 'dnaid']].set_index('sampleID').drop_duplicates()
    sample_data = sample_data.loc[sample_data.index.intersection(good_samples)]
    expr_data = filter_samples(exp_df, good_samples)
    expr_data = filter_inoculum(expr_data, filter_below=filter_below)
    expr_data = expr_data[list(sample_data.index)].reset_index()
    return sample_data, expr_data


def calculate_fitness(edf, sdf):
    dds = py_DESeq2(count_matrix=edf,
                    design_matrix=sdf,
                    design_formula='~ day',
                    gene_column='barcode')  # <- telling DESeq2 this should be the gene ID column

    dds.run_deseq()
    days = list(sdf['day'].unique())
    days.remove('d0')
    all_results = []
    for d in days:
        dds.get_deseq_result(contrast=['day', d, 'd0'])
        res = dds.deseq_result
        res['day'] = d
        all_results.append(res)
    fdf = pd.concat(all_results)
    n_samples = sdf.groupby('day').mouse.nunique().to_dict()
    fdf['n_samples'] = fdf.day.map(n_samples)
    return fdf


def calculate_2dist_zscore(u1, s1, u2, s2):
    return (u1 - u2) / np.sqrt((s1 ** 2) + (s2 ** 2))


def calculte_comparisons(fitness, df, controls, cntrl_type='wt'):
    """

    fitness: DESeq2 output, log2FoldChange value for each barcode comparing each time point with inoculum
    df: df for 1 experiment and 1 dnaid

    """
    days = list(df['day'].unique())
    days.remove('d0')
    # Get all entries that were mapped to a gene
    gene_bc = set(df[df.gene != '-'].barcode.values)
    gene_df = fitness.loc[fitness.index.intersection(gene_bc)]
    # Add gene annotation to the fitness table
    gene_df = gene_df.merge(df[['barcode', 'gene']], how='left', on='barcode').drop_duplicates()
    # Calculate mean log2FoldChange and sigma for each gene (for all )
    gene_mean = gene_df.groupby(['gene', 'day']).agg({'log2FoldChange': ['mean'], 'lfcSE': [sigma]}).reset_index()
    gene_mean.columns = ['gene', 'day', 'gene_FC', 'sigma']

    # Get all the WITS barcodes
    controls_bc = set(controls[controls.phenotype == 'wt'].barcode.values)
    cntrl_df = fitness.loc[fitness.index.intersection(controls_bc)]
    # Calculate mean log2FoldChange and sigma for the control barcodes (for all barcodes)
    cntrl_mean = cntrl_df.groupby(['day']).agg({'log2FoldChange': ['mean'], 'lfcSE': [sigma]})
    cntrl_mean.columns = ['cntrl_FC', 'cntrl_sigma']
    cntrl_mean = cntrl_mean.reset_index()
    # Calculate zscore and competitive index (CI) for each gene
    gene_mean = gene_mean.merge(cntrl_mean, how='left', on='day')
    gene_mean['zscore'] = gene_mean.apply(
        lambda x: calculate_2dist_zscore(x['gene_FC'], x['sigma'], x['cntrl_FC'], x['cntrl_sigma']), axis=1)
    gene_mean['ci'] = gene_mean.apply(lambda x: 2 ** x['gene_FC'] / 2 ** x['cntrl_FC'], axis=1)
    gene_mean = gene_mean[['gene', 'day', 'zscore', 'ci']]

    # Get all barcodes that were not mapped to a gene
    others_bc = set(df[(df.gene == '-')].barcode.values)
    other_df = fitness.loc[fitness.index.intersection(others_bc)]
    other_df = other_df.merge(cntrl_mean, how='left', on='day')
    # Calculate zscore and CI for each barcode
    other_df['zscore'] = other_df.apply(
        lambda x: calculate_2dist_zscore(x['log2FoldChange'], x['lfcSE'], x['cntrl_FC'], x['cntrl_sigma']), axis=1)
    other_df['ci'] = other_df.apply(lambda x: 2 ** x['log2FoldChange'] / 2 ** x['cntrl_FC'], axis=1)
    other_df = other_df[['barcode', 'day', 'zscore', 'ci']].rename({'barcode': 'gene'}, axis=1)

    # Concatenate the gene and barcode results
    results = pd.concat([gene_mean, other_df])
    # Calculate p-values for the genes/barcodes
    results['pval'] = results.zscore.apply(lambda x: scipy.stats.norm.sf(abs(x)) * 2)

    # Spread results over days
    results = (results.pivot(index=['gene'], columns=['day'], values=['zscore', 'ci', 'pval'])
               .reset_index())
    # Rename columns
    results.columns = ['gene'] + [f'{day}_zscore' for day in days] + [f'{day}_ci' for day in days] + [f'{day}_pval' for
                                                                                                      day in days]
    # Adjust p-values for multiple testing
    for day in days:
        results[f'{day}_padj'] = \
        statsmodels.stats.multitest.multipletests(results[f'{day}_pval'], alpha=0.05, method='fdr_bh')[1]

    return results.set_index('gene')


def to_list(x):
    bc_list = list(x)
    if len(bc_list) == 1:
        return bc_list[0]
    return ", ".join(list(x))


def final_fitness_table(fitness, exp_df):
    barcode_info = exp_df[['barcode', 'locus', 'gene', 'library']].drop_duplicates().set_index('barcode')
    fit2 = fitness.merge(barcode_info, how='left', left_index=True, right_index=True)

    fit2_gene = fit2[(fit2.locus.notnull()) & (fit2.locus != '-')]
    fit3 = fit2_gene.groupby(['library', 'gene', 'locus', 'day', 'n_samples', ]).agg(
        {'barcode': ['count', to_list], 'log2FoldChange': ['mean', 'std']}).reset_index()
    fit3.columns = ['library', 'gene', 'locus', 'day', 'num_samples', 'num_barcodes', 'barcode', 'mean_fitness',
                    'std_fitness']

    fit4 = (fit3.pivot(index=['gene', 'locus', 'num_barcodes', 'barcode', 'library'], columns=['day'],
                       values=['mean_fitness', 'std_fitness', 'num_samples'])
            .reset_index())

    #     # if not mapped to gene
    fit2_no_gene = fit2[(fit2.locus.isna()) | (fit2.locus == '-')]
    fit_no_gene = fit2_no_gene.groupby(['barcode', 'library', 'day', 'n_samples', ]).agg(
        {'barcode': ['count', to_list], 'log2FoldChange': ['mean', 'std']}).reset_index()
    fit_no_gene.columns = ['gene', 'library', 'day', 'num_samples', 'num_barcodes', 'barcode', 'mean_fitness',
                           'std_fitness']

    fit_no_gene = (fit_no_gene.pivot(index=['gene', 'num_barcodes', 'barcode', 'library'], columns=['day'],
                                     values=['mean_fitness', 'std_fitness', 'num_samples'])
                   .reset_index())

    days = list(exp_df.day.unique())
    days.remove('d0')
    day_columns = [f'{day}_fitness' for day in days] + [f'{day}_std_fitness' for day in days] + [f'{day}_num_samples'
                                                                                                 for day in days]
    fit4.columns = ['gene', 'locus', 'num_barcodes', 'barcode', 'library'] + day_columns
    fit_no_gene.columns = ['gene', 'num_barcodes', 'barcode', 'library'] + day_columns

    fit = pd.concat([fit4, fit_no_gene])
    fit = fit.merge(exp_df[['barcode', 'position', 'seq']].drop_duplicates(), how='left', on='barcode')
    cnrls = get_control_fitness(fitness, controls)
    cnrls['gene'] = cnrls.index

    return pd.concat([fit, cnrls]).set_index('gene')


def get_control_fitness(fitness, controls):
    cnrls = controls.merge(fitness, how='left', on='barcode').assign(control='yes').drop(['DN'], axis=1)
    cnrls = cnrls[['barcode', 'phenotype', 'conc', 'log2FoldChange', 'n_samples', 'day', 'control']]
    cnrls = cnrls[cnrls.day.notnull()]
    days = list(cnrls.day.unique())
    cnrls = cnrls.pivot(index=['barcode', 'phenotype', 'conc', 'control'], columns=['day'],
                        values=['log2FoldChange', 'n_samples']).reset_index()
    day_columns = [f'{day}_fitness' for day in days] + [f'{day}_num_samples' for day in days]
    cnrls.columns = ['barcode', 'phenotype', 'conc', 'control'] + day_columns
    return cnrls.set_index('barcode')


def process_results(exp_df, controls, good_samples, filter_below=1000, cntrl_type='wt'):
    print('Generating dataset for DESeq')
    sdf, edf = generate_DE_dataset(exp_df, good_samples, filter_below=filter_below)
    print('Calculating Fitness')
    fitness = calculate_fitness(edf, sdf)
    print('Calculating Z-Scores')
    results = calculte_comparisons(fitness, exp_df, controls, cntrl_type)
    print('Merging Results')
    fit_tab = final_fitness_table(fitness, exp_df)
    final = fit_tab.merge(results, how='outer', left_index=True, right_index=True)
    days = list(sdf['day'].unique())
    days.remove('d0')
    col_order = ['num_samples', 'fitness', 'std_fitness', 'ci', 'zscore', 'pval', 'padj']
    final_cols = ['locus', 'num_barcodes', 'library', 'barcode', 'position', 'seq'] + [f'{d}_{s}' for d in days for s in
                                                                                       col_order] + ['control',
                                                                                                     'phenotype',
                                                                                                     'conc']
    print("Done")
    return fit_tab, final[final_cols]


def run_dnaid(sample_df, controls, dnaid, filter_below=1000, cntrl_type='wt', cutoff=0.9 ):
    fresults = []
    experiments = list(sample_df[sample_df.dnaid==dnaid].experiment.unique())
    for exp in experiments:
        exp_df = subset_experiment(sample_df, dnaid, exp)
        wits, good_samples = calculate_correlation(exp_df, controls, cutoff=cutoff)
        print(exp)
        print(good_samples)
        fit, final = process_results(exp_df, controls, good_samples, filter_below=filter_below, cntrl_type = cntrl_type)
        fresults.append(final.assign(experiment=exp))
    return pd.concat(fresults)
