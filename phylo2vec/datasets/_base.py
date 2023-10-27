from importlib import resources

from Bio import SeqIO

DATA_MODULE = "phylo2vec.datasets.data"
DESCR_MODULE = "phylo2vec.datasets.descr"


def load_fasta(data_file_name, data_module=DATA_MODULE):
    return SeqIO.parse(_open_text(data_module, data_file_name), "fasta")


def load_descr(descr_file_name, descr_module=DESCR_MODULE):
    return _read_text(descr_module, descr_file_name)


# Largely inspired by https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/utils/fixes.py
def _open_text(module, file_name, encoding="utf-8"):
    return resources.files(module).joinpath(file_name).open("r", encoding=encoding)


# Largely inspired by https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/utils/fixes.py
def _read_text(module, file_name, encoding="utf-8"):
    return resources.files(module).joinpath(file_name).read_text(encoding=encoding)
