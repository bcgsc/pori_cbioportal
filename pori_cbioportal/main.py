import argparse
import logging
import os

from .study import generate_reports
from .util import LOG_LEVELS


def command_interface():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--study_id', help='The study ID (ex. "BRCA TCGA PANCAN (2018)")')
    parser.add_argument(
        'data_folder', help='path to the folder containing the cbioportal export files'
    )
    parser.add_argument(
        '--strict',
        action='store_true',
        default=False,
        help='Flag to indicate the program should stop if any report fails to upload',
    )
    parser.add_argument(
        '--patient_data',
        default='data_clinical_patient.txt',
        help='The clinical patient data for this study',
    )
    parser.add_argument(
        '--sample_data',
        default='data_clinical_sample.txt',
        help='The clinical sample data for this study',
    )
    parser.add_argument(
        '--expression',
        default='data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt',
        help='The mRNA RNA Seq zscore data.',
    )
    parser.add_argument(
        '--small_mutations',
        default='data_mutations_extended.txt',
        help='The small mutation data',
    )
    parser.add_argument(
        '--discrete_cna',
        default='data_CNA.txt',
        help='The discrete copy variant data',
    )
    parser.add_argument(
        '--continuous_cna',
        help='The continuous (log2cna) copy variant data',
        default='data_log2CNA.txt',
    )
    parser.add_argument(
        '--ipr_project', help='The IPR project to upload this report to', default='TEST'
    )
    parser.add_argument('--fusions', '-f', help='The fusion data', default='data_fusions.txt')
    parser.add_argument(
        '--username',
        default=os.environ.get('USER'),
        required=bool(not os.environ.get('USER')),
        help='The username for logging in to GraphKB and IPR',
    )
    parser.add_argument(
        '--debugging_filename',
        help='path to write the JSON that was attempted to upload to IPR to on a failure',
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

    def resolve_data_path(filepath, required=True):
        """
        Join the path to the data folder if it is relative, don't if it is absolute
        """

        if filepath.startswith('/') and os.path.exists(os.path.join('/', filepath)):
            return filepath
        new_path = os.path.join(args.data_folder, filepath)
        if not os.path.exists(new_path):
            if required:
                raise argparse.ArgumentTypeError(f'file does not exist: {new_path} ({filepath})')
            return None
        return new_path

    logging.basicConfig(
        level=LOG_LEVELS['info'],
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        datefmt='%m-%d-%y %H:%M:%S',
    )

    generate_reports(
        args.study_id if args.study_id else os.path.basename(args.data_dir),
        resolve_data_path(args.patient_data),
        resolve_data_path(args.sample_data),
        continuous_copy_variants_filename=resolve_data_path(args.continuous_cna, False),
        discrete_copy_variants_filename=resolve_data_path(args.discrete_cna),
        small_mutations_filename=resolve_data_path(args.small_mutations),
        expression_filename=resolve_data_path(args.expression, False),
        fusions_filename=resolve_data_path(args.fusions, False),
        username=args.username,
        password=args.password,
        ipr_project=args.ipr_project,
        ipr_url=args.ipr_url,
        graphkb_url=args.graphkb_url,
        patients_subset=args.patients_subset,
        strict=args.strict,
        debugging_filename=args.debugging_filename,
    )
