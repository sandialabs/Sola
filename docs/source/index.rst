WOLF: Waanders Outer-Loop Framework
===================================

Overview
--------

WOLF is the **W**\ aanders **O**\ uter-\ **L**\ oop **F**\ ramework for iteratively solving large-scale problems such as optimization problems constrained by partial differential equations.
This project includes documentation and tutorials for two separate software ecosystems:

1. The **S**\ andbox for **A**\ djoint-\ **B**\ ased outer **L**\ oop analysis :doc:`SABL <./sabl/about>`, a `MATLAB <https://www.mathworks.com/products/matlab.html>`_ framework designed for prototyping outer-loop problems.
2. The **M**\ ulti-\ **r**\ esolution **Hy**\ bridizable **D**\ ifferential **E**\ quations solver :doc:`MrHyDE <./mrhyde/about>`, a powerful `C++ <https://cplusplus.com/>`_ package for solving partial differential equations at scale.

This project shows how to use SABL and MrHyDE to solve the following types of problems.

- :doc:`Constrained Optimization <./optimization>`, including PDE-constrained optimization.
- :doc:`Bayesian Inversion <./bayesinversion>`.
- :doc:`Optimal Experimental Design (OED) <./oed>`.
- :doc:`Hyper-differential Sensitivity Analysis (HDSA) <./hdsa>`.

.. important::
   This documentation is aimed at advanced undergraduates and graduate students with some MATLAB experience and light familiarity with C++.
   The goal is to lower the barrier to entry for large-scale high-performance computing tasks using sophisticated software libraries.
   See the [Fundamentals](./fundamentals) page for the necessary mathematical background and references to materials for getting up to speed.

Contents
--------

.. toctree::
   :maxdepth: 1
   :caption: Outer-loop Problems

   optimization
   bayesinversion
   oed
   hdsa

.. toctree::
   :maxdepth: 1
   :caption: MATLAB: SABL

   sabl/about
   sabl/installation
   sabl/anatomy
   sabl/example1

.. toctree::
   :maxdepth: 1
   :caption: C++: MrHyDE

   mrhyde/about
   mrhyde/installation
   mrhyde/primer
   mrhyde/anatomy
   mrhyde/example1

.. toctree::
   :maxdepth: 2
   :caption: Appendix

   fundamentals
   notation
