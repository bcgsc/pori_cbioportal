import argparse
import logging
import os

from .study import generate_reports
from .util import LOG_LEVELS


def command_interface():
    parser = argparse.ArgumentParser()
    parser.add_argument('study_id', help='The study ID (ex. dlbc_tcga_pan_can_atlas_2018)')
    parser.add_argument(
        'patient_data',
        help='The clinical patient data for this study (ex. data_clinical_patient.txt)',
    )
    parser.add_argument(
        'sample_data',
        help='The clinical sample data for this study (ex. data_clinical_sample.txt)',
    )
    parser.add_argument(
        '--expression',
        '-e',
        help='The mRNA RNA Seq zscore data. (ex. data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt)',
        required=True,
    )
    parser.add_argument(
        '--small_mutations',
        '-m',
        help='The small mutation data (ex. data_mutations_extended.txt)',
        required=True,
    )
    parser.add_argument(
        '--discrete_cna',
        '-d',
        help='The discrete copy variant data (ex. data_CNA.txt)',
        required=True,
    )
    parser.add_argument(
        '--continuous_cna',
        '-c',
        help='The continuous (log2cna) copy variant data (ex. data_log2CNA.txt)',
        required=True,
    )
    parser.add_argument(
        '--ipr_project', help='The IPR project to upload this report to', default='TEST'
    )
    parser.add_argument(
        '--fusions', '-f', help='The fusion data (ex. data_fusions.txt)', required=True,
    )
    parser.add_argument(
        '--username',
        default=os.environ.get('USER'),
        required=bool(not os.environ.get('USER')),
        help='The username for logging in to GraphKB and IPR',
    )
    parser.add_argument(
        '--password', required=True, help='The password for logging in to GraphKB and IPR'
    )
    parser.add_argument(
        '--ipr_url', default='https://iprstaging-api.bcgsc.ca/api', help='The IPR API url'
    )
    parser.add_argument(
        '--graphkb_url',
        default='https://graphkbstaging-api.bcgsc.ca/api',
        help='The GraphKB API url',
    )
    parser.add_argument(
        '--patients_subset', nargs='+', help='Only create reports for a subset of patient IDs'
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=LOG_LEVELS['info'],
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        datefmt='%m-%d-%y %H:%M:%S',
    )

    generate_reports(
        args.study_id,
        args.patient_data,
        args.sample_data,
        continuous_copy_variants_filename=args.continuous_cna,
        discrete_copy_variants_filename=args.discrete_cna,
        small_mutations_filename=args.small_mutations,
        expression_filename=args.expression,
        fusions_filename=args.fusions,
        username=args.username,
        password=args.password,
        ipr_project=args.ipr_project,
        ipr_url=args.ipr_url,
        graphkb_url=args.graphkb_url,
        patients_subset=args.patients_subset,
    )
