from setuptools import find_packages, setup

# Dependencies required to use your package
INSTALL_REQS = ['ipr >=3.0.0, <4.0.0', 'pandas', 'graphkb', 'matplotlib', 'seaborn']

# Dependencies required for development
DEV_REQS = ['flake8', 'black', 'mypy']

# Dependencies required only for running tests
TEST_REQS = ['pytest', 'pytest-cov']

# Dependencies required for deploying to an index server
DEPLOYMENT_REQS = ['twine', 'wheel', 'm2r']

long_description = ''
long_description_content_type = 'text/markdown'

try:
    import re

    import m2r

    long_description = m2r.parse_from_file('README.md')
    long_description = re.sub(
        r'.. code-block::.*', '.. code::', long_description
    )  # pyshop has issues with fenced code blocks
    long_description_content_type = 'text/rst'
except ImportError:
    with open('README.md', 'r') as fh:
        long_description = fh.read()


setup(
    name='pori_cbioportal',
    version='0.1.1',
    packages=find_packages(),
    install_requires=INSTALL_REQS,
    extras_require={
        'dev': TEST_REQS + DEPLOYMENT_REQS + DEV_REQS,
        'deploy': DEPLOYMENT_REQS,
        'test': TEST_REQS,
    },
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    python_requires='>=3',
    author='creisle',
    author_email='creisle@bcgsc.ca',
    maintainer='creisle',
    maintainer_email='creisle@bcgsc.ca',
    dependency_links=[],
    test_suite='tests',
    tests_require=TEST_REQS,
    entry_points={'console_scripts': ['pori_cbioportal=pori_cbioportal.main:command_interface']},
)
