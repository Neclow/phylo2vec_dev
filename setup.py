from setuptools import find_packages, setup

setup(
    name="phylo2vec",
    version="0.1.0",
    description="Phylo2Vec: integer vector representation of binary (phylogenetic) trees",
    author="Neil Scheidwasser",
    author_email="neil.clow@sund.ku.dk",
    url="https://github.com/Neclow/phylo2vec",
    packages=find_packages,
    install_requires=[
        "numba==0.56.4",
        "numpy==1.23.5",
        "biopython==1.80.0",
        "python==3.10",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: Unix",
    ],
)