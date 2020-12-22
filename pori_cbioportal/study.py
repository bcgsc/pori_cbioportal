"""
Process study download data files
"""
import re
from typing import Iterable

import pandas
from graphkb.match import INPUT_COPY_CATEGORIES, INPUT_EXPRESSION_CATEGORIES
from ipr import main
from ipr.connection import IprConnection

from .expression import load_zscore_data, upload_expression_density_plots
from .util import logger

GENE_NAME = 'Hugo_Symbol'
GENE_ID = 'Entrez_Gene_Id'


def load_copy_variants(filename_discrete, filename_log2cna=None) -> pandas.DataFrame:
    logger.info(f'reading: (discrete CNA) {filename_discrete}')
    discrete_df = pandas.read_csv(
        filename_discrete, delimiter='\t', dtype={GENE_NAME: 'string', GENE_ID: 'string'}
    )
    discrete_df['type'] = 'discrete'

    if filename_log2cna:
        logger.info(f'reading: (continuous CNA) {filename_log2cna}')
        log2cna_df = pandas.read_csv(
            filename_log2cna, delimiter='\t', dtype={GENE_NAME: 'string', GENE_ID: 'string'}
        )
        log2cna_df['type'] = 'log2cna'

        assert log2cna_df.columns.tolist() == discrete_df.columns.tolist()
        copy_varints_df = pandas.concat([log2cna_df, discrete_df])
    else:
        copy_varints_df = discrete_df
    copy_varints_df = copy_varints_df.rename(columns={GENE_NAME: 'gene', GENE_ID: 'gene_id'})
    return copy_varints_df


def load_small_mutations(filename, ploidy=2) -> pandas.DataFrame:
    logger.info(f'reading: {filename}')
    mutations_df = pandas.read_csv(
        filename,
        delimiter='\t',
        dtype={
            't_alt_count': float,
            't_ref_count': float,
            'n_ref_count': float,
            'n_alt_count': float,
            'n_depth': float,
            't_depth': float,
            'Chromosome': 'string',
            'HGVSg': 'string',
            'HGVSc': 'string',
            'SYMBOL': 'string',
            'Transcript_ID': 'string',
            'HGVSp_Short': 'string',
            'Tumor_Sample_Barcode': 'string',
            'Matched_Norm_Sample_Barcode': 'string',
            'Reference_Allele': 'string',
            'Tumor_Seq_Allele1': 'string',
            'Start_Position': float,
            'End_Position': float,
        },
        na_values=[
            '#N/A',
            '#N/A N/A',
            '#NA',
            '-1.#IND',
            '-1.#QNAN',
            '-NaN',
            '-nan',
            '1.#IND',
            '1.#QNAN',
            '<NA>',
            'N/A',
            'NA',
            'NULL',
            'NaN',
            'n/a',
            'nan',
            'null',
            '.',
        ],
    )
    mutations_df = mutations_df.rename(
        columns={
            'Chromosome': 'chromosome',
            'HGVSg': 'hgvsGenomic',
            'HGVSc': 'hgvsCds',
            'SYMBOL': 'gene',
            'Transcript_ID': 'transcript',
            'HGVSp_Short': 'proteinChange',
            'Tumor_Sample_Barcode': 'sample_id',
            'Matched_Norm_Sample_Barcode': 'normalLibrary',
            't_alt_count': 'tumourAltCount',
            't_ref_count': 'tumourRefCount',
            'n_ref_count': 'normalRefCount',
            'n_alt_count': 'normalAltCount',
            'n_depth': 'normalDepth',
            't_depth': 'tumourDepth',
            'Reference_Allele': 'refSeq',
            'Tumor_Seq_Allele1': 'altSeq',
            'Start_Position': 'startPosition',
            'End_Position': 'endPosition',
            'NCBI_Build': 'ncbiBuild',
        }
    )

    def normalize_protein_change(proteinChange):
        if pandas.isnull(proteinChange):
            return ''
        if proteinChange.endswith('_splice'):
            proteinChange = proteinChange.replace('_splice', 'spl')
        if not proteinChange.startswith('p.'):
            proteinChange = 'p.' + proteinChange
        return proteinChange

    mutations_df['proteinChange'] = mutations_df['proteinChange'].apply(normalize_protein_change)

    def hgvs_protein(row):
        if pandas.isnull(row['proteinChange']) or not row['proteinChange'].strip():
            return ''
        return row['gene'] + ':' + row['proteinChange']

    def reject_bad_cds(cds):
        if pandas.isnull(cds):
            return ''
        return cds

    mutations_df['hgvsProtein'] = mutations_df.apply(hgvs_protein, axis=1)
    mutations_df['hgvsCds'] = mutations_df.hgvsCds.apply(reject_bad_cds)
    mutations_df['refSeq'] = mutations_df.refSeq.replace('-', '')
    mutations_df['altSeq'] = mutations_df.altSeq.replace('-', '')
    mutations_df['hgvsGenomic'] = ''  # TODO: Compose genomic hgvs notation where possible
    mutations_df = mutations_df[
        [
            'gene',
            'chromosome',
            'transcript',
            'refSeq',
            'altSeq',
            'tumourRefCount',
            'tumourAltCount',
            'tumourDepth',
            'normalRefCount',
            'normalAltCount',
            'normalDepth',
            'startPosition',
            'endPosition',
            'hgvsGenomic',
            'hgvsCds',
            'hgvsProtein',
            'proteinChange',
            'sample_id',
            'ncbiBuild',
        ]
    ]

    def choose_main_variant(row):
        if row.proteinChange:
            return row.proteinChange
        if row.hgvsCds:
            return row.hgvsCds.split(':')[1]
        if row.hgvsGenomic:
            return row.hgvsGenomic.split(':')[1]
        return row.proteinChange

    mutations_df['proteinChange'] = mutations_df.apply(choose_main_variant, axis=1)

    return mutations_df.drop_duplicates()


def load_fusions(filename) -> pandas.DataFrame:
    logger.info(f'reading: {filename}')
    mutations_df = pandas.read_csv(
        filename,
        delimiter='\t',
        dtype={
            'frame': 'string',
            'Tumor_Sample_Barcode': 'string',
            'Fusion': 'string',
            'DNA_support': 'string',
            'RNA_support': 'string',
        },
    )
    mutations_df = mutations_df.rename(
        columns={'Frame': 'frame', 'Tumor_Sample_Barcode': 'sample_id', 'Fusion': 'fusion'}
    )

    def gene1(fusion):
        return fusion.split('-')[0]

    def gene2(fusion):
        if '-' not in fusion:
            return ''
        return '-'.join(fusion.split('-')[1:])

    mutations_df['gene1'] = mutations_df.fusion.apply(gene1)
    mutations_df['gene2'] = mutations_df.fusion.apply(gene2)

    def omic_support(detected_in):
        return bool('DNA' in detected_in and 'RNA' in detected_in)

    def detected_in(row):
        if pandas.isnull(row['DNA_support']) or pandas.isnull(row['RNA_support']):
            return ''
        raise NotImplementedError(
            'these values have never been not blank yet. need to implement handling'
        )

    # TODO: determine if these are ever not null and what format they take
    mutations_df['detectedIn'] = mutations_df.apply(detected_in, axis=1)
    mutations_df['omicSupport'] = mutations_df.detectedIn.apply(omic_support)

    mutations_df = mutations_df[['gene1', 'gene2', 'frame', 'sample_id']]
    mutations_df['exon1'] = ''
    mutations_df['exon2'] = ''
    mutations_df['breakpoint'] = ''
    mutations_df['eventType'] = 'fusion'
    return mutations_df.drop_duplicates()


def load_clinical_data(patients_filename, samples_filename) -> pandas.DataFrame:
    logger.info(f'reading: {patients_filename}')
    patients_df = pandas.read_csv(patients_filename, delimiter='\t', comment='#')
    patients_df = patients_df.rename(
        columns={
            'SEX': 'gender',
            'PATIENT_ID': 'patientId',
            'SUBTYPE': 'subType',
            'OTHER_PATIENT_ID': 'alternateIdentifier',
        }
    )

    def days_to_age(days):
        if pandas.isnull(days):
            return None
        else:
            return abs(int(round(days / 365, 0)))

    patients_df['age'] = patients_df.DAYS_TO_BIRTH.apply(days_to_age)
    patients_df = patients_df.reset_index()
    logger.info(f'reading: {samples_filename}')
    samples_df = pandas.read_csv(samples_filename, delimiter='\t', comment='#')
    samples_df = samples_df.rename(
        columns={
            'PATIENT_ID': 'patientId',
            'SAMPLE_ID': 'sample_id',
            'CANCER_TYPE_DETAILED': 'diagnosis',
            'Tissue Source Site': 'biopsySite',
        }
    )
    samples_df = samples_df.reset_index()
    clinical_df = patients_df[
        ['gender', 'patientId', 'subType', 'age', 'alternateIdentifier']
    ].merge(
        samples_df[['patientId', 'sample_id', 'diagnosis']],
        left_on='patientId',
        right_on='patientId',
        how='left',
        suffixes=(False, False),
    )
    return clinical_df


def create_report(
    study_id: str,
    patient_id: str,
    sample_id: str,
    clinical_df: pandas.DataFrame,
    expression_df: pandas.DataFrame,
    small_mutations_df: pandas.DataFrame,
    copy_variants_df: pandas.DataFrame,
    fusions_df: pandas.DataFrame,
    gene_conflicts: Iterable[str],
    username: str,
    password: str,
    ipr_url: str,
    zscore_threshold=2,
    percentile_threshold=90,
    copy_amplification_threshold=2,
    copy_homd_threshold=-2,
    **kwargs,
):
    def categorize_expression(row):
        if row.diseaseZScore >= zscore_threshold:
            if row.diseasePercentile >= percentile_threshold:
                return INPUT_EXPRESSION_CATEGORIES.UP
        elif row.diseaseZScore <= (zscore_threshold * -1):
            if row.diseasePercentile <= 100 - percentile_threshold:
                return INPUT_EXPRESSION_CATEGORIES.DOWN
        return ''

    if sample_id in expression_df.columns:
        expression = expression_df[['gene', 'gene_id', sample_id, 'type']].copy()
        expression = expression[~expression.gene.isin(gene_conflicts)]
        expression = expression[~pandas.isnull(expression.gene)]
        expression = pandas.pivot_table(
            expression, columns=['type'], values=[sample_id], index=['gene']
        ).reset_index()
        # flatten the multi-index
        expression.columns = [tup[-1] if tup[-1] else tup[-2] for tup in expression.columns.values]
        expression = expression.rename(
            columns={'percentile': 'diseasePercentile', 'zscore': 'diseaseZScore'}
        )
        expression['kbCategory'] = expression.apply(categorize_expression, axis=1)
        expression = expression.fillna('')
        expression = expression.drop_duplicates(['gene', 'kbCategory']).to_dict('records')
    else:
        logger.warning(f'no expression data found for sample ({sample_id})')
        expression = []

    small_mutations = small_mutations_df[small_mutations_df['sample_id'] == sample_id].copy()
    small_mutations = small_mutations.drop(columns=['sample_id'])
    small_mutations = small_mutations[~small_mutations.gene.isin(gene_conflicts)]
    small_mutations = small_mutations.to_dict('records')

    def categorize_copy_change(copy_change):
        if copy_change >= copy_amplification_threshold:
            return INPUT_COPY_CATEGORIES.AMP
        elif copy_change <= copy_homd_threshold:
            return INPUT_COPY_CATEGORIES.DEEP
        elif copy_change > 0:
            return INPUT_COPY_CATEGORIES.GAIN
        elif copy_change < 0:
            return INPUT_COPY_CATEGORIES.LOSS
        else:
            return ''

    if sample_id in copy_variants_df.columns:
        copy_variants = copy_variants_df[['gene', sample_id, 'type']].copy()
        copy_variants = copy_variants.rename(columns={sample_id: 'copyChange'})
        discrete = copy_variants[copy_variants.type == 'discrete'].drop(columns=['type'])
        continuous = copy_variants[copy_variants.type == 'log2cna'].drop(columns=['type'])
        continuous = continuous.rename(columns={'copyChange': 'log2Cna'})
        copy_variants = discrete.merge(continuous, on=['gene'], how='left')
        copy_variants['kbCategory'] = copy_variants['copyChange'].apply(categorize_copy_change)
        copy_variants = copy_variants[~copy_variants.gene.isin(gene_conflicts)]
        copy_variants = copy_variants.drop_duplicates(
            ['gene', 'kbCategory', 'copyChange', 'log2Cna']
        )
        copy_variants = copy_variants.fillna('')
        copy_variants = copy_variants.to_dict('records')
    else:
        copy_variants = []

    fusions = fusions_df[fusions_df.sample_id == sample_id].copy().drop(columns=['sample_id'])
    fusions = fusions[~fusions.gene1.isin(gene_conflicts)]
    fusions = fusions[~fusions.gene2.isin(gene_conflicts)]
    fusions = fusions.fillna('')
    fusions = fusions.to_dict('records')

    clinical = clinical_df[clinical_df.sample_id == sample_id].to_dict('records')
    if len(clinical) != 1:
        raise NotImplementedError(
            f'unable to handle cases ({patient_id}) without or with multiple clinical records per sample ({sample_id}_'
        )
    clinical = clinical[0]
    content = main.create_report(
        patient_id=patient_id,
        kb_disease_match=clinical['diagnosis'],
        project='TEST',
        expression_variant_rows=expression,
        copy_variant_rows=copy_variants,
        structural_variant_rows=fusions,
        small_mutation_rows=small_mutations,
        optional_content={
            'comparators': [{'analysisRole': 'expression (disease)', 'name': study_id}],
            'patientInformation': {
                'gender': clinical.get('gender'),
                'diagnosis': re.sub(r' \(NOS\)$', ', NOS', clinical.get('diagnosis')),
                'physician': '',
                'caseType': 'not specified',
                'biopsySite': clinical.get('biopsySite'),
                'age': clinical.get('age'),
            },
            'alternateIdentifier': clinical.get('alternateIdentifier'),
            'biopsyName': sample_id,
            'subtyping': clinical.get('subtyping'),
            'ploidy': clinical.get('ploidy'),
        },
        username=username,
        password=password,
        ipr_url=ipr_url,
        **kwargs,
    )
    ipr_conn = IprConnection(username, password, ipr_url)

    upload_expression_density_plots(ipr_conn, expression_df, sample_id, content)


def find_conflicting_gene_names(
    continuous_copy_variants_filename,
    discrete_copy_variants_filename,
    small_mutations_filename,
    expression_filename,
    fusions_filename,
):
    """
    Cbioportal mixes Ensembl, Entrez, and HGNC gene definitions. Since only the gene names (HGNC)
    are used consistently throughout the variant files we must remove any genes that have conflicting
    definitions based solely on the gene name
    """
    genes_df = pandas.read_csv(small_mutations_filename, delimiter='\t')
    genes_df = genes_df[
        ['SYMBOL']
    ].copy()  # Gene in small mutations is ensembl ID which is not helpful
    genes_df = genes_df.rename(columns={'SYMBOL': GENE_NAME})

    for filename in [
        discrete_copy_variants_filename,
        continuous_copy_variants_filename,
        expression_filename,
    ]:
        df = pandas.read_csv(filename, delimiter='\t')
        df = df[[GENE_NAME, GENE_ID]].copy()
        genes_df = pandas.concat([genes_df, df])

    df = pandas.read_csv(fusions_filename, delimiter='\t')
    genes_df = pandas.concat([genes_df, df[[GENE_NAME, GENE_ID]].copy()])
    df[[GENE_NAME, 'Hugo_Symbol2']] = df.Fusion.str.split('-', 1, expand=True)
    genes_df = pandas.concat([genes_df, df[[GENE_NAME]].copy()])
    genes_df = pandas.concat(
        [genes_df, df[['Hugo_Symbol2']].rename(columns={'Hugo_Symbol2': GENE_NAME}).copy()]
    )
    genes_df = genes_df.drop_duplicates()
    genes_df = genes_df.dropna(subset=[GENE_ID])
    genes_df = genes_df.dropna(subset=[GENE_NAME])
    conflicts = genes_df[genes_df.duplicated(GENE_NAME)][GENE_NAME].tolist()
    return conflicts


def generate_reports(
    study_id,
    patients_filename,
    samples_filename,
    continuous_copy_variants_filename,
    discrete_copy_variants_filename,
    small_mutations_filename,
    expression_filename,
    fusions_filename,
    username,
    password,
    ipr_url,
    **kwargs,
):
    logger.info(f'generating study ({study_id}) reports')
    clinical_df = load_clinical_data(patients_filename, samples_filename,)

    gene_conflicts = find_conflicting_gene_names(
        continuous_copy_variants_filename,
        discrete_copy_variants_filename,
        small_mutations_filename,
        expression_filename,
        fusions_filename,
    )
    logger.warning(
        f'ignoring {len(gene_conflicts)} gene names with conflicting definitions: {", ".join(sorted(list(gene_conflicts)))}'
    )

    expression_df = load_zscore_data(expression_filename)

    small_mutations_df = load_small_mutations(small_mutations_filename)
    copy_variants_df = load_copy_variants(
        discrete_copy_variants_filename, continuous_copy_variants_filename,
    )
    fusions_df = load_fusions(fusions_filename)

    for index, row in clinical_df.iterrows():
        logger.info(f'creating a report for {row["patientId"]} {row["sample_id"]}')
        try:
            create_report(
                study_id,
                row['patientId'],
                row['sample_id'],
                clinical_df,
                expression_df,
                small_mutations_df,
                copy_variants_df,
                fusions_df,
                gene_conflicts,
                username,
                password,
                ipr_url,
                **kwargs,
            )
        except Exception as err:
            logger.error(err)
            logger.info('skipping to the next report')
