---
title: 'phylo2vec: a library for vector-based phylogenetic tree manipulation'
tags:
  - Rust
  - Python
  - R
  - bioinformatics
  - phylogenetics
  - binary tree
authors:
  - name: Neil Scheidwasser
    orcid: 0000-0001-9922-0289
    affiliation: 1 # (Multiple affiliations must be quoted)
    equal-contrib: true
  - name: Ayush Nag
    affiliation: 2
    equal-contrib: true
  - name: Matthew J Penn
    orcid: 0000-0001-8682-5393
    affiliation: 4
  - name: Anthony MV Jakob
    affiliation: 5
  - name: Mark P Khurana
    affiliation: 1
  - name: Don Setiawan
    affiliation: 2
  - name: Madeline Gordon
    affiliation: 2
  - name: David A Duchêne
    orcid: 0000-0002-5479-1974
  - name: Samir Bhatt
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: "1, 4"
    orcid: 0000-0002-0891-4611
affiliations:
 - name: Section of Health Data Science and AI, University of Copenhagen, Copenhagen
   index: 1
 - name: eScience Institute, University of Washington, Seattle, United States
   index: 2
 - name: Department of Statistics, University of Oxford, Oxford, United Kingdom
   index: 3
 - name: MRC Centre for Global Infectious Disease Analysis, Imperial College London, London, United Kingdom
   index: 4
 - name: Independent researcher
   index: 5
date: 24 March 2024
geometry: margin=2cm
bibliography: paper.bib
---

# Summary

Phylogenetics is a cornerstone of many analysis frameworks in computational biology [@yang2014] and linguistics [@atkinson2005]. Recently, the advent of large-scale genomics and the SARS-CoV-2 pandemic has underscored the necessity to scale phylogenetic software to handle large datasets of genomes or phylogenetic trees [@kapli2020; @attwood2022; @khurana2024; @kraemer2025]. While significant efforts have focused on scaling phylogenetic inference [@turakhia2021; @sanderson2021; @demaio2023], visualization [@sanderson2022], and lineage identification [@mcbroome2024], a burgeoning body of research has been dedicated to efficient representations of data for genomes [@deorowicz2023] and phylogenetic trees such as phylo2vec [@penn2024], HOP [@chauve2025], and OLA [@richman2025]. Compared to traditional tree representations such as the Newick format [@felsenstein2004], which describes a phylogenetic tree as a string of nested parentheses enclosing pairs of leaves or subtrees, these modern representations utilize integer vectors to describe tree topology traversal. This approach offers several advantages, including increased easier manipulability, memory efficiency, and applicability to downstream tasks such as machine learning.

Here, we present the second version of \texttt{phylo2vec}, a high-performance software package for encoding, manipulating, and analyzing binary phylogenetic trees.  At its core, the package is based on the phylo2vec [@penn2024] representation of binary trees [@penn2024], which bijects any tree topology with $n$ leaves into an integer vector of size $n-1$. Compared to the traditional Newick format, phylo2vec was designed with fast sampling and rapid tree comparison in mind. This version features a core implementation in Rust, providing significant performance improvements and memory efficiency \textcolor{red}{figure 1: benchmark}, while remaining available in Python (superseding the version described in the original paper [@penn2024]) and R via dedicated wrappers, making it accessible to a broad audience in the bioinformatics community.

# Statement of need

The purpose of the \texttt{phylo2vec} library is threefold. First, the Rust-based core of the library aims at providing a robust phylogenetic tree manipulation library in Rust, complementing other efforts such as \texttt{light\_phylogeny} [@duchemin2018], which focuses on tree visualization and manipulation of reconciled phylogenies [@nakhleh2013], and \texttt{rust-bio} [@koster2016], a comprehensive bioinformatics library which does not yet cover phylogenetics. Second, \texttt{phylo2vec} aims at complementing existing phylogenetic libraries such as \texttt{ape} [@paradis2019] in R and \texttt{ete3} [@huerta2016] and \texttt{DendroPy} [@Moreno2024] in Python with fast tree sampling, fast tree comparison and efficient tree data compression [@penn2024]. Third, the inherent tree representation of phylo2vec offers a pathway to gradient-based optimization frameworks for phylogenetic inference. A notable example is GradME [@penn2023], which relaxes the vector representation of phylo2vec into a continuous space.

# Features

The presented version of \texttt{phylo2vec} addresses several limitations of [@penn2024]. In particular, it allows for branch length annotations, extending the vector representation of size $n-1$ to a matrix of size $(n-1) \times 3$, where $n$ denotes the number of leaves (or taxa) in a tree \textcolor{red}{Figure 2: example of phylo2vec with branch lengths}
and a **$\mathcal{O}(n \log n)$** implementation of vector-to-Newick conversion based on Adelson-Velsky and Landis (AVL; [@adelson1962]) trees. Moreover, the current release features several new additions, including fast cophenetic distance matrix calculation and a skeleton for Bayesian phylogenetic inference using Markov Chain Monte Carlo (MCMC). The inference framework leverages similarities between phylo2vec and \texttt{BEAGLE}'s inner representation of post-order traversal [@suchard2009; @ayres2012]. Lastly, user-friendliness is enhanced by a comprehensive documentation using \texttt{Jupyterbook} and step-by-step demos of the inner workings underlying phylo2vec's vector representation.

# Maintenance

\textcolor{red}{For VISS}: paragraph to write about implemented development practices. See <https://joss.theoj.org/papers/10.21105/joss.06943>

# Acknowledgements

SB acknowledges support from the Danish National Research Foundation via a chair grant (DNRF160) which also supports NS. SB acknowledges support from The Eric and Wendy Schmidt Fund For Strategic Innovation via the Schmidt Polymath Award (G-22-63345).

# References
