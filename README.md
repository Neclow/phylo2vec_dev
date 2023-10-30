# Phylo2Vec
This repository contains an implementation of Phylo2Vec. It is distributed under the GNU Lesser General Public License v3.0 (LGPL).

Link to the preprint: https://arxiv.org/abs/2304.12693

## Installation
### Dependencies:
* numba==0.56.4
* numpy==1.23.5
* biopython==1.80.0
* python==3.10

## Additional test dependencies
* ete3==3.1.3
* pytest==7.4.2
* six==1.16.0

### User installation:
TODO

## Development
### Testing
After installation, you can launch the test suite from outside the source directory (you will need to have pytest >= 7.4.2 installed):
```
pytest phylo2vec
```

## Citation and other work
```latex
@article{phylo2vec,
  title={Phylo2Vec: a vector representation for binary trees},
  author={Penn, Matthew J and Scheidwasser, Neil and Khurana, Mark P and Duch{\^e}ne, David A and Donnelly, Christl A and Bhatt, Samir},
  journal={arXiv preprint arXiv:2304.12693},
  year={2023}
}
```
* Preprint repository (core functions are deprecated): https://github.com/Neclow/phylo2vec_preprint
* C++ version (deprecated): https://github.com/Neclow/phylo2vec_cpp
* Related work: https://github.com/Neclow/GradME = phylo2vec + minimum evolution + gradient descent
* The file ```100trees.txt``` is an anonmyised subset of [TimeTree](http://www.timetree.org/) trees adapted from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.2fqz612vg. Kudos to [Mark Khurana](https://forskning.ku.dk/soeg/result/?pure=da%2Fpersons%2Fmark-poulsen-khurana(171ece7e-9567-4d48-8cf9-959b57de57c8).html) for the dataset
