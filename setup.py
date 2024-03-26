from setuptools import setup, find_packages

setup(
    name="tormentor",
    version='0.0.1',
    packages=find_packages(),
    author="Frederico Schmitt Kremer",
    author_email="fred.s.kremer@gmail.com",
    description="Tormentor",
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    keywords="bioinformatics",
    entry_points = {'console_scripts':[
        'tormentor = tormentor.main:main'
        ]},
    install_requires = [
        'pandas',
        'biopython',
        'bcbio-gff'
    ]
)