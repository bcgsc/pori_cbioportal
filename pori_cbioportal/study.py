"""
Process study download data files
"""
import re
from typing import Dict, Iterable, List, Optional

import pandas
from graphkb.match import INPUT_COPY_CATEGORIES, INPUT_EXPRESSION_CATEGORIES
from ipr import main
from ipr.connection import IprConnection

from .expression import load_zscore_data, upload_expression_density_plots
from .util import add_optional_columns, logger, read_csv

GENE_NAME = 'Hugo_Symbol'
GENE_ID = 'Entrez_Gene_Id'


def load_copy_variants(
    filename_discrete: str, filename_log2cna: Optional[str] = None
) -> pandas.DataFrame:
    discrete_df = read_csv(filename_discrete, dtype={GENE_NAME: 'string', GENE_ID: 'string'})
    discrete_df['type'] = 'discrete'

    if filename_log2cna:
        log2cna_df = read_csv(filename_log2cna, dtype={GENE_NAME: 'string', GENE_ID: 'string'})
        log2cna_df['type'] = 'log2cna'

        assert log2cna_df.columns.tolist() == discrete_df.columns.tolist()
        copy_varints_df = pandas.concat([log2cna_df, discrete_df])
    else:
        copy_varints_df = discrete_df
    copy_varints_df = copy_varints_df.rename(columns={GENE_NAME: 'gene', GENE_ID: 'gene_id'})
    return copy_varints_df


def load_small_mutations(filename, **kwargs) -> pandas.DataFrame:
    mutations_df = read_csv(
        filename,
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
        **kwargs,
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
            'Allele': 'altSeq',
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

    def hgvs_cds(row):
        if pandas.isnull(row.hgvsCds):
            return ''
        if ':' not in row.hgvsCds:
            return f'{row.transcript}:{row.hgvsCds}'
        return row.hgvsCds

    mutations_df['hgvsProtein'] = mutations_df.apply(hgvs_protein, axis=1)
    mutations_df['hgvsCds'] = mutations_df.apply(hgvs_cds, axis=1)
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
            return row.hgvsCds.split(':')[1] if ':' in row.hgvsCds else row.hgvsCds
        if row.hgvsGenomic:
            return row.hgvsGenomic.split(':')[1] if ':' in row.hgvsGenomic else row.hgvsGenomic
        return row.proteinChange

    mutations_df['proteinChange'] = mutations_df.apply(choose_main_variant, axis=1)

    return mutations_df.drop_duplicates()


def load_fusions(filename, **kwargs) -> pandas.DataFrame:
    mutations_df = read_csv(
        filename,
        dtype={
            'frame': 'string',
            'Tumor_Sample_Barcode': 'string',
            'Fusion': 'string',
            'DNA_support': 'string',
            'RNA_support': 'string',
        },
        **kwargs,
    )
    mutations_df = mutations_df.rename(
        columns={'Frame': 'frame', 'Tumor_Sample_Barcode': 'sample_id', 'Fusion': 'fusion'}
    )

    gene_pattern = r'^(?P<gene1>[A-Za-z0-9_]+(-(\d+|AS1))?)(-(?P<gene2>\S+))?(\s+(?P<eventType>\S+)(\s-\sArcher)?)?$'
    gene_df = mutations_df.fusion.str.extract(gene_pattern).fillna('').replace('nan', '')
    gene_df.loc[gene_df.eventType == '', ['eventType']] = 'fusion'
    # not certain what they mean by truncation, using the broader "fusion" type
    gene_df.loc[gene_df.eventType == 'truncation', ['eventType']] = 'fusion'

    gene_df['gene2'] = gene_df.apply(
        lambda row: row.gene2 if row.gene2 != 'intragenic' else row.gene1, axis=1
    )
    mutations_df = pandas.concat([mutations_df, gene_df], axis=1)

    def omic_support(detected_in):
        return bool('DNA' in detected_in and 'RNA' in detected_in)

    def detected_in(row):
        if pandas.isnull(row['DNA_support']) or pandas.isnull(row['RNA_support']):
            return ''

        if row['DNA_support'] == 'yes' and row['RNA_support'] == 'unknown':
            return 'DNA'
        if row['DNA_support'] == 'yes' and row['RNA_support'] == 'yes':
            return 'DNA/RNA'
        if row['DNA_support'] == 'unknown' and row['RNA_support'] == 'unknown':
            return ''
        if row['DNA_support'] == 'unknown' and row['RNA_support'] == 'yes':
            return 'RNA'
        print(row)
        raise NotImplementedError(
            'these values have never been not blank yet. need to implement handling'
        )

    # TODO: determine if these are ever not null and what format they take
    mutations_df['detectedIn'] = mutations_df.apply(detected_in, axis=1)
    mutations_df['omicSupport'] = mutations_df.detectedIn.apply(omic_support)

    mutations_df = mutations_df[['gene1', 'gene2', 'frame', 'sample_id', 'eventType']]
    mutations_df['exon1'] = ''
    mutations_df['exon2'] = ''
    mutations_df['breakpoint'] = ''
    return mutations_df.drop_duplicates()


def load_clinical_data(patients_filename, samples_filename) -> pandas.DataFrame:
    patients_df = read_csv(patients_filename)

    add_optional_columns(patients_df, ['OTHER_PATIENT_ID', 'SUBTYPE', 'SEX', 'DAYS_TO_BIRTH'], '')
    patients_df = patients_df.rename(
        columns={
            'SEX': 'gender',
            'PATIENT_ID': 'patientId',
            'SUBTYPE': 'subType',
            'OTHER_PATIENT_ID': 'alternateIdentifier',
        }
    )

    def days_to_age(days):
        if pandas.isnull(days) or days == '':
            return ''
        else:
            return str(abs(int(round(days / 365, 0))))

    patients_df['age'] = patients_df.DAYS_TO_BIRTH.apply(days_to_age)
    patients_df = patients_df.reset_index()
    samples_df = read_csv(samples_filename)
    add_optional_columns(samples_df, ['Tissue Source Site'], '')
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


def replace_values(rows: List[Dict], *fields: str, current_value='', replacement_value=None):
    for row in rows:
        for field in fields:
            if row[field] == current_value:
                row[field] = replacement_value
    return rows


def drop_fields_with_value(rows: List[Dict], *fields: str, value=''):
    for row in rows:
        for field in fields:
            if field in row and row[field] == value:
                del row[field]
    return rows


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
    ipr_project: str = 'TEST',
    zscore_threshold=2,
    percentile_threshold=97.5,
    copy_amplification_threshold=2,
    copy_homd_threshold=-2,
    study_size: Optional[int] = None,
    debugging_filename: Optional[str] = None,
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

    if expression_df is not None and sample_id in expression_df.columns:
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
        expression = expression.drop_duplicates(['gene', 'kbCategory']).to_dict('records')
    else:
        logger.warning(f'no expression data found for sample ({sample_id})')
        expression = []

    small_mutations = small_mutations_df[small_mutations_df['sample_id'] == sample_id].copy()
    small_mutations = small_mutations.drop(columns=['sample_id'])
    small_mutations = small_mutations[~small_mutations.gene.isin(gene_conflicts)]
    small_mutations = small_mutations.fillna('').to_dict('records')

    def categorize_copy_change(copy_change: int):
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
        copy_variants = drop_fields_with_value(
            copy_variants.fillna('').to_dict('records'), 'log2Cna'
        )
    else:
        copy_variants = []

    if fusions_df is not None:
        fusions = fusions_df[fusions_df.sample_id == sample_id].copy().drop(columns=['sample_id'])
        fusions = fusions[~fusions.gene1.isin(gene_conflicts)]
        fusions = fusions[~fusions.gene2.isin(gene_conflicts)]
        fusions = fusions.fillna('')
        fusions = fusions.to_dict('records')
    else:
        fusions = []

    clinical = clinical_df[clinical_df.sample_id == sample_id].to_dict('records')
    if len(clinical) != 1:
        raise NotImplementedError(
            f'unable to handle cases ({patient_id}) without or with multiple clinical records per sample ({sample_id}_'
        )
    clinical = clinical[0]
    content = main.create_report(
        generate_therapeutics=True,
        output_json_path=debugging_filename,
        content={
            'expressionVariants': replace_values(expression, 'kbCategory'),
            'copyVariants': replace_values(copy_variants, 'kbCategory'),
            'structuralVariants': replace_values(fusions, 'exon1', 'exon2'),
            'smallMutations': small_mutations,
            'comparators': [
                {'analysisRole': 'expression (disease)', 'name': study_id, 'size': study_size}
            ],
            'patientInformation': {
                'gender': clinical.get('gender'),
                'diagnosis': re.sub(r' \(NOS\)$', ', NOS', clinical['diagnosis']),
                'physician': 'not specified',
                'caseType': 'not specified',
                'biopsySite': clinical.get('biopsySite', ''),
                'age': clinical.get('age', ''),
            },
            'project': ipr_project,
            'kbDiseaseMatch': clinical['diagnosis'],
            'patientId': patient_id,
            'alternateIdentifier': clinical.get('alternateIdentifier'),
            'biopsyName': sample_id,
            'subtyping': clinical.get('subtyping', ''),
            'ploidy': clinical.get('ploidy', ''),
            'template': 'TCGA-cBioportal',
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
    genes_df = read_csv(small_mutations_filename)
    genes_df = genes_df[
        ['SYMBOL']
    ].copy()  # Gene in small mutations is ensembl ID which is not helpful
    genes_df = genes_df.rename(columns={'SYMBOL': GENE_NAME})

    for filename in [
        discrete_copy_variants_filename,
        continuous_copy_variants_filename,
        expression_filename,
    ]:
        if not filename:
            continue
        df = read_csv(filename)
        if GENE_ID in df.columns:
            df = df[[GENE_NAME, GENE_ID]].copy()
            genes_df = pandas.concat([genes_df, df])

    if fusions_filename:
        df = read_csv(fusions_filename)
        if GENE_ID in df.columns:
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
    study_id: str,
    patients_filename: str,
    samples_filename: str,
    continuous_copy_variants_filename: str,
    discrete_copy_variants_filename: str,
    small_mutations_filename: str,
    expression_filename: str,
    fusions_filename: str,
    username: str,
    password: str,
    ipr_url: str,
    ipr_project: str,
    patients_subset: Optional[List[str]] = [],
    strict: bool = False,
    **kwargs,
):
    logger.info(f'generating study ({study_id}) reports')
    clinical_df = load_clinical_data(patients_filename, samples_filename)

    gene_conflicts = find_conflicting_gene_names(
        continuous_copy_variants_filename,
        discrete_copy_variants_filename,
        small_mutations_filename,
        expression_filename,
        fusions_filename,
    )
    if gene_conflicts:
        logger.warning(
            f'ignoring {len(gene_conflicts)} gene names with conflicting definitions: {", ".join(sorted(list(gene_conflicts)))}'
        )

    if expression_filename:
        expression_df = load_zscore_data(expression_filename)
    else:
        expression_df = None

    small_mutations_df = load_small_mutations(small_mutations_filename)
    copy_variants_df = load_copy_variants(
        discrete_copy_variants_filename,
        continuous_copy_variants_filename,
    )
    if fusions_filename:
        fusions_df = load_fusions(fusions_filename)
    else:
        fusions_df = None

    patients_filter = {p.lower() for p in (patients_subset or [])}
    sample_count = clinical_df.sample_id.nunique()
    for _, row in clinical_df.iterrows():
        if patients_filter and row.patientId.lower() not in patients_filter:
            logger.warning(
                f'skipping patient {row.patientId} not in patients selected subset {patients_subset}'
            )
            continue
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
                study_size=sample_count,
                ipr_project=ipr_project,
                username=username,
                password=password,
                ipr_url=ipr_url,
                **kwargs,
            )
        except Exception as err:
            if strict:
                raise err
            logger.error(err)
            logger.info('skipping to the next report')
