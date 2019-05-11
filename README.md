# Variation Quantum Factoring

## Introduction

This repository contains implementation of the algorithm presented in the article "Variational Quantum Factoring", by Eric R. Anschuetz, Jonathan P. Olson, Alán Aspuru-Guzik, Yudong Cao. It's available on [arxiv](https://arxiv.org/abs/1808.08927).

The notation in the code refers directly to the notation in the paper.

## Differences from the paper

Below you can find a list of points where I'm aware that my implementation differs from the one provided in the paper.

### General

- To perform simulations I decided to use pyQuil instead of QuTiP.

### Preprocessing
- Number of preprocessing rules is higher than what comes from the equations (5) in the paper, since some of them were not stated explicitly ("trivial relations"). 
- The results suggest that this script does more preprocessing than what has been done in the paper. 
- Algorithm can be ran with (as in paper, see footnote 38) or without prior knowledge about the length of the numbers p and q. 
- In order to use custom preprocessing rules for numbers 56153 and 291311 (see footnote 40), specific functions must be called. Otherwise, the regular preprocessing scheme will be applied.

### QAOA
- The grid size used for the initialization of BFGS algorithm (Table I in the paper) seems to be chosen arbitrary. Therefore I have also used arbitrary grid size - it might be too low, depending on the number being factored.
- Implementation of the BFGS algorithm is different from the one used in the original paper (see report in `research/2019_05_08_performance_checks`).


## Research

Research performed using this implementation can be found in the `research` directory. I follow the convention presented [here](https://github.com/BOHRTECHNOLOGY/public_research/blob/master/Guidelines/research_guidelines.md).


## Tests 

To run tests please run `python -m pytest` from the main directory.


## Issues

### Randomness in results

For some reason, the preprocessing is not deterministic. Running the same case sporadically leads to getting different results. It seems to come from the fact, that there is some randomness inside sympy when it comes to ordering operations. Hence, preprocessing sometimes assigns `x=y` and sometimes `y=x`, which leads to different expressions after substition. From the mathematical point of view it doesn't matter, but from the practical - it does. Since set of implemented rules is incomplete, different substitution may occasionally lead to form which this algorithm cannot simplify. This effect diminished over time of development - i.e. improving rules and fixing bugs.

### Known issues

I do not claim that the preprocessing part is perfect, though from manual inspection it seems to be working in most cases. Below are some known bugs.

- Current version of code doesn't produce correct results for number 1465 (factors: 293 and 5). This doesn't happen always (see note about randomness above).
- There are still some additional rules to add / cases to fix (see TODO in `preprocessing.py`).
- In cases exhibiting some form of symmetry (as described [here](https://arxiv.org/pdf/1411.6758.pdf)), procedure of calculating squared overlap might give wrong results. The fix for numbers 56153 and 291311, has been hardcoded, but it's far from being elegant and general solution.
