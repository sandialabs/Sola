# The Waanders Outer-Loop Framework

## Overview

WOLF is the **W**aanders **O**uter-**L**oop **F**ramework for iteratively solving large-scale problems such as optimization problems constrained by partial differential equations.
This project includes documentation and tutorials for two separate software ecosystems:

1. The **S**andbox for **A**djoint-**B**ased outer **L**oop analysis [(SABL)](./sabl/about), a [MATLAB](https://www.mathworks.com/products/matlab.html) framework designed for prototyping outer-loop problems.
2. The **M**ulti-**r**esolution **Hy**bridizable **D**ifferential **E**quations solver [(MrHyDE)](./mrhyde/about), a powerful [C++](https://cplusplus.com) package for solving partial differential equations at scale.

This project shows how to use SABL and MrHyDE to solve the following types of problems.

- [Constrained Optimization](./problems/optimization), including PDE-constrained optimization.
- [Bayesian Inversion](./problems/bayesinversion).
- [Optimal Experimental Design (OED)](./problems/oed).
- [Hyper-differential Sensitivity Analysis (HDSA)](./problems/hdsa).

:::{important}
This project is designed for advanced undergraduates and graduate students with some MATLAB experience and light familiarity with C++.
The goal is to lower the barrier to entry for large-scale high-performance computing tasks using sophisticated software libraries.
See the [Fundamentals](./appendix/fundamentals) page for the necessary mathematical background and references to materials for getting up to speed.
:::

## Contents

```{eval-rst}
.. toctree::
   :maxdepth: 1
   :caption: Outer-loop Problems

   problems/optimization
   problems/bayesinversion
   problems/oed
   problems/hdsa

.. toctree::
   :maxdepth: 1
   :caption: MATLAB: SABL

   sabl/about
   sabl/installation
   sabl/guides
   sabl/tutorials
   sabl/api

.. toctree::
   :maxdepth: 1
   :caption: C++: MrHyDE

   mrhyde/about
   mrhyde/installation
   mrhyde/primer
   mrhyde/guides
   mrhyde/tutorials

.. toctree::
   :maxdepth: 2
   :caption: Appendix

   appendix/fundamentals
   appendix/notation
   appendix/contributing
```
