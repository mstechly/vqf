# Variation Quantum Factoring

## Introduction

This repository contains implementation of the algorithm presented in the article "Variational Quantum Factoring", by Eric R. Anschuetz, Jonathan P. Olson, Alán Aspuru-Guzik, Yudong Cao. It's available on [arxiv](https://arxiv.org/abs/1808.08927).

The notation in the code refers directly to the notation in the paper.

## Differences from the paper

Below you can find a list of points where I'm aware that my implementation differs from the one provided in the paper.

- To perform simulations I decided to use pyQuil instead of QuTiP.
- Number of preprocessing rules is higher than what comes from the equation (5) in the paper, since some of them were not stated explicitly ("trivial relations").
- The grid size used for the initialization of BFGS algorithm (Table I in the paper) seems to be chosen arbitrary. Therefore I have also used arbitrary grid size - it might be too low, depending on the number being factored.
