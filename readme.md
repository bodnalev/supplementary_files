# Supplementary files

This repository holds supplementary files for papers involving flag algebra calculations. Each calculation uses a modified version of sage, it is available from a separate GitHub repository [here](https://github.com/bodnalev/sage). Installation instructions and a quick guide is included. The package is still under development.

## Papers

1. `no_c5m/` folder contains calculations for [The Turan density of the tight 5-cycle minus one edge](https://arxiv.org/pdf/2412.21011),
2. `hypercube_ind_stats/` folder contains calculations for [Some exact values of the inducibility and statistics constants for
hypercubes](https://arxiv.org/pdf/2503.03408),
3. `no_cl/` folder contains calculations for [The Turan density of short tight cycles](https://arxiv.org/pdf/2506.03223),
4. `graph_ind/` folder contains calculations for [Some exact inducibility-type results for graphs via flag algebras](https://arxiv.org/abs/2507.01596).

## File structure

Each repository contains jupter notebooks showcasing the script that generated the results and the outputs. They can be easily inspected from a browser, without a sage installation. Additionally a certificates folder is included where for each claim a pickled object contains the relevant data to recover the proof and verify it. This data (where possible, see below), can be inspected with any python 3 instance.

## Contents of the certificates

Each certificate file contains a dictionary with the following information:

1. `result` - final bound the certificate proves.
2. `target size` - maximum size of flags used in the proof.
3. `target` - the vector the proof maximizes or minimizes.
4. `maximize` - a boolean value encoding if it is a maximization or minimization problem.
5. `positives` - a list of vectors encoding the positivity assumptions. 
6. `base flags` - list of untyped flags used in the proof, each with size equal to `target size`.
7. `typed flags` - a dictionary, with keys from `(n, tau, N)`, where `tau` is a type, `n, N` are numbers with `N` equal to `target size`. The corresponding value in the dictionary is a list of typed flags with size `n` and type `tau`.
8. `X matrices` - a list of flattened positive semi-definite matrices used in the sum of squares proof. The upper diagonal entries of the symmetric matrix are flattened to a row in reading order. The order of these matrices and the indexing of the columns agrees with the order of the types and typed flags in the `typed flags` entry.
9. `e vector` - a vector of positive coefficients indicating how each `positives` vector assumption was used.
10. `slack vector` - a vector of positive values indicating the difference between the bound and the element-wise proved values.
11. `phi vectors` - (optional) list of possible blow-up densities, either provided during the generation of the result to help the numerical optimization, or the construction found the the SDP solver.

## Data format

When the values in the matrices/vectors are python compatible, then they are stored in that format. This includes rational numbers in python's `Fractions` class, and floating points numbers in the default `float` class (for inexact results). In this case the certificate can be inspected and verified in pure python 3. 

However, some of the exact proofs work over field extensions of the rational numbers, in this case the certificate data contains sage values. These certificates can be inspected in any version of sage above 9.0.

The custom version of [sage](https://github.com/bodnalev/sage) used to generate these results includes a verifier for all the certificate types, however some of the results include a separate independently written verifier by [Jared Leon](https://github.com/WozMit) in pure python 3.
