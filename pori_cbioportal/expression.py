import os
import tempfile
from typing import Dict

import pandas
import seaborn
from ipr.connection import IprConnection
from matplotlib import pyplot as plt

from .util import logger, read_csv

GENE_NAME = 'Hugo_Symbol'
GENE_ID = 'Entrez_Gene_Id'


def load_zscore_data(filename: str) -> pandas.DataFrame:
    df = read_csv(filename, dtype={GENE_NAME: 'string', GENE_ID: 'string'})
    df = df.rename(columns={GENE_NAME: 'gene', GENE_ID: 'gene_id'})
    df['type'] = 'zscore'
    # calculate the percentiles
    ranks = df.copy().set_index('gene_id')
    ranks = ranks.rank(0, pct=True, numeric_only=True).apply(lambda x: round(x * 100))
    ranks['type'] = 'percentile'
    ranks['gene_id'] = ranks.index
    ranks = ranks.reset_index(drop=True)
    ranks = ranks.merge(df[['gene_id', 'gene']], on='gene_id')

    df = pandas.concat([df, ranks])
    return df


def plot_expression_density(
    expression_df: pandas.DataFrame, sample_id: str, gene: str, plot_name: str
) -> None:
    """
    Draw the expression density plot for a given gene in a given sample
    """
    logger.info(f'generating expression density plot for {gene}')
    plt.figure()
    df = expression_df[expression_df.gene == gene]
    df = df[df.type == 'zscore']
    sample_columns = [c for c in expression_df.columns if c not in ['gene', 'gene_id', 'type']]
    df = df[sample_columns]
    current = df[sample_id].values[0]
    values = df.values.tolist()[0]
    seaborn.set_style("whitegrid")
    ax = seaborn.distplot(values, axlabel="RNA z-score")

    # add the sample marker
    for p in ax.patches:
        if current >= p.get_x() and current <= p.get_x() + p.get_width():
            p.set_color('black')
            plt.text(
                p.get_x() + (p.get_width() / 2),
                p.get_height(),
                '*',
                horizontalalignment='center',
                verticalalignment='bottom',
            )
    ax.figure.savefig(plot_name, bbox_inches="tight", dpi=600)
    plt.close()


def upload_expression_density_plots(
    ipr_conn: IprConnection, expression_df: pandas.DataFrame, sample_id: str, content: Dict
) -> None:
    """
    Given a report that has been created, generate expression density reports
    for expression outliers which have matches in the kb
    """
    kb_expresssion_matches = set()
    for match in content.get('kbMatches', []):
        if match['variantType'] == 'exp':
            kb_expresssion_matches.add(match['variant'])

    expression_genes_to_plot = set()
    for exp in content.get('expressionVariants', []):
        if exp['key'] in kb_expresssion_matches:
            expression_genes_to_plot.add(exp['gene'])

    logger.info(f'making {len(expression_genes_to_plot)} expression plots')
    if not expression_genes_to_plot:
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        files = {}
        data = {}

        for gene in sorted(list(expression_genes_to_plot)):
            plot_name = os.path.join(tmpdir, f'{gene}.png')
            plot_expression_density(expression_df, sample_id, gene, plot_name)
            key = f'expDensity.{gene}'
            files[key] = plot_name
            data[f'{key}_title'] = f'{gene} RNA Expression Z-Scores'
            data[
                f'{key}_caption'
            ] = f'Cohort RNA Expression values of {gene}. The asterisk indicates the bin containing the expression value for this patient.'
        report_id = content['ident']
        ipr_conn.post_images(report_id, files, data)
