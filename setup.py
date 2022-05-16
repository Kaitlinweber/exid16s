import io
import os
# import sys
# from shutil import rmtree

from setuptools import find_packages, setup, Command

NAME = 'exid16s'
DESCRIPTION = 'Package to extract and identify the 16S sequence.'
URL = 'https://github.com/Kaitlinweber/exid16s'
EMAIL = 'kaitlin@live.nl'
AUTHOR = 'Kaitlin Weber'
REQUIRES_PYTHON = '>=3.9'
VERSION = '0.1.0'


REQUIRED = [
    'biopython', 'pandas'
]


EXTRAS = {
    # 'fancy feature': ['django'],
}


EXCLUDE = [
    "test", "*.test", "*.test.*", "test.*"
]


here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION


setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=['exid16s'],
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='AGPL3',
    classifiers=[
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    entry_points={
        'console_scripts': [
            'exid16s = exid16s.__main__:main',
        ],
    },
)
