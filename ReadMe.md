This repository contains all the code to reproduce the results in [The dual PC algorithm and the role of Gaussianity for structure learning of Bayesian networks](https://www.sciencedirect.com/science/article/pii/S0888613X23001068).

- dualPC.R contains the R code for the main functions.
- dualPC_sim.R contains an example run with the different PC versions; it calls dualPC_algs.R which contains the individual runs.
- DAGfns.R contains functions for generating data and comparing the results.
- sims_collated/ contains the collected data from the simulations.
- plots/ contains the generated plots and the code used to create them.

Example
-------

```
require(datasets)
source("dualPC.R")
cor_mat <- cor(trees)
N <- nrow(trees)
dual_pc(cor_mat, N, alpha = 0.05) # Returns adjacency matrix of estimated CPDAG
```
Reference
---------

```
@inproceedings{NEURIPS2023_b146e7c8,
 author = {Giudice, Enrico and Kuipers, Jack and Moffa, Giusi},
 booktitle = {Advances in Neural Information Processing Systems},
 editor = {A. Oh and T. Neumann and A. Globerson and K. Saenko and M. Hardt and S. Levine},
 pages = {56602--56614},
 publisher = {Curran Associates, Inc.},
 title = {A Bayesian Take on Gaussian Process Networks},
 url = {https://proceedings.neurips.cc/paper_files/paper/2023/file/b146e7c87685fa208bd95ce4b08e330c-Paper-Conference.pdf},
 volume = {36},
 year = {2023}
}
```
