import logging

import pandas

# name the logger after the package to make it simple to disable for packages using this one as a dependency
# https://stackoverflow.com/questions/11029717/how-do-i-disable-log-messages-from-the-requests-library
VERBOSE_ERROR_CODE = (logging.INFO + logging.DEBUG) // 2
logging.addLevelName(VERBOSE_ERROR_CODE, 'VERBOSE')
logger = logging.getLogger('pori_cbioportal')
# add shortbut for verbose logging
setattr(logger, 'verbose', lambda *pos, **kw: logger.log(VERBOSE_ERROR_CODE, *pos, **kw))
LOG_LEVELS = {
    'info': logging.INFO,
    'debug': logging.DEBUG,
    'warn': logging.WARN,
    'error': logging.ERROR,
    'verbose': VERBOSE_ERROR_CODE,
}


def add_optional_columns(df, columns, default_value=''):
    for column in columns:
        if column not in df.columns:
            df[column] = default_value


def read_csv(
    filename: str, comment: str = '#', delimiter: str = '\t', **kwargs
) -> pandas.DataFrame:
    """
    calls the pandas read_csv function with some comment defaults
    """
    logger.info(f'reading: {filename}')
    df = pandas.read_csv(filename, delimiter=delimiter, comment=comment, **kwargs)
    return df
