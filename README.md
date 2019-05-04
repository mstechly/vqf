# Variation Quantum Factoring

## Introduction

This repository contains implementation of the algorithm presented in the article "Variational Quantum Factoring", by Eric R. Anschuetz, Jonathan P. Olson, Alán Aspuru-Guzik, Yudong Cao. It's available on [arxiv](https://arxiv.org/abs/1808.08927).

The notation in the code refers directly to the notation in the paper.

### Warning

For some reason, the preprocessing is not deterministic. Running the same case sporadically leads to getting different results. It seems to come from the fact, that there is some randomness inside sympy when it comes to ordering operations. Hence, preprocessing sometimes assigns `x=y` and sometimes `y=x`, which leads to different expressions after substition. From the mathematical point of view it doesn't matter, but from the practical - it does. Since set of implemented rules is incomplete, different substitution may occasionally lead to form which this algorithm cannot simplify. This effect diminished over time of development - i.e. improving rules and fixing bugs.

## Differences from the paper

Below you can find a list of points where I'm aware that my implementation differs from the one provided in the paper.

### General

- To perform simulations I decided to use pyQuil instead of QuTiP.

### Preprocessing
- Number of preprocessing rules is higher than what comes from the equations (5) in the paper, since some of them were not stated explicitly ("trivial relations"). 
- The results suggest that this script does more preprocessing than what has been done in the paper. 
- Algorithm can be ran with (as in paper, see footnote 38) or without prior knowledge about the length of the numbers p and q. 
- I have not implemented custom preprocessing rules for numbers 56153 and 291311, see footnote 40.

### QAOA
- The grid size used for the initialization of BFGS algorithm (Table I in the paper) seems to be chosen arbitrary. Therefore I have also used arbitrary grid size - it might be too low, depending on the number being factored.

## Tests 

To run tests please run `python -m pytest` from the main directory.