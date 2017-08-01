import os
import re
from setuptools import find_packages

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

version = __import__('pbsv_polish').get_version()

_REQUIREMENTS_FILE = 'REQUIREMENTS.txt'
_README = 'README.md'


def _get_local_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


def _get_description(file_name):
    with open(file_name, 'r') as f:
        _long_description = f.read()
    return _long_description


def _get_requirements(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    rx = re.compile('^[A-z]')
    requirements = [l for l in lines if rx.match(l) is not None]
    if "READTHEDOCS" in os.environ:
        requirements = [r for r in requirements if not "pbcore" in r]
    return requirements

setup(
    name='pbsv_polish',
    version=version,
    package_dir={'pbsv_polish': 'pbsv_polish'},
    packages=find_packages('.'),
    license='BSD',
    author='awenger,yli',
    author_email='awenger@pacificbiosciences.com, yli@pacificbiosciences.com',
    description='PacBio structure variants polishing tool.',
    entry_points={'console_scripts': [
        'sv_pbdagcon = pbsv_polish.sv_pbdagcon:main',
        'polish_sv = pbsv_polish.polish_sv:main',
        'collect_polished_sv = pbsv_polish.collect_polished_sv:main',
        'substr_fasta = pbsv_polish.substr_fasta:main',
        'make_precise_sv_from_refs = pbsv_polish.make_precise_sv_from_refs:main'
    ]},
    install_requires=_get_requirements(_get_local_file(_REQUIREMENTS_FILE)),
    tests_require=['pytest'],
    long_description=_get_description(_get_local_file(_README)),
    classifiers=['Development Status :: 4 - Beta'],
    include_package_data=True,
    zip_safe=False
)
