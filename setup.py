from setuptools import setup

requires = [
    'click',
    'biopython>=1.7',
    'pandas'
    ]

setup(
    name="tnseq2",
    version="1.0",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("Bioinformatic toolkit for barcode mapping and counting"),
    license="LICENSE",
    keywords="tnseq",
    install_requires=requires,
    #url = "https://github.com/SushiLab/TNSEQ_DEV",
    packages=['tnseq2'],
    entry_points={
        'console_scripts': ['tnseq2=tnseq2.main:main'],
    }
)
