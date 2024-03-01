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
@article{GIUDICE2023108975,
  title = {The dual PC algorithm and the role of Gaussianity for structure learning of Bayesian networks},
  journal = {International Journal of Approximate Reasoning},
  volume = {161},
  pages = {108975},
  year = {2023},
  issn = {0888-613X},
  doi = {https://doi.org/10.1016/j.ijar.2023.108975},
  url = {https://www.sciencedirect.com/science/article/pii/S0888613X23001068},
  author = {Enrico Giudice and Jack Kuipers and Giusi Moffa}
}
```
