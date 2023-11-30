# API

This section documents the public API for the classes defined in SABL.

## Optimization

```{eval-rst}
.. toctree::
   :maxdepth: 1

   api/optimization_Objective
   api/optimization_Dynamic_Objective
   api/optimization_Constraint
   api/optimization_Dynamic_Constraint
   api/optimization_Reduced_Space_Optimization
```

### Automatic Differentiation

```{eval-rst}
.. toctree::
   :maxdepth: 1

   api/automatic_differentiation_Objective_AD.md
   api/automatic_differentiation_Dynamic_Objective_AD.md
   api/automatic_differentiation_Constraint_AD.md
   api/automatic_differentiation_Dynamic_Constraint_AD.md
```

### Model Reduction

```{eval-rst}
.. toctree::
   :maxdepth: 1

   api/model_reduction_Basis
   api/model_reduction_POD_Basis
   api/model_reduction_OpInf_Operator
   api/model_reduction_Constant_Operator
   api/model_reduction_Linear_Operator
   api/model_reduction_Quadratic_Operator
   api/model_reduction_Input_Operator
   api/model_reduction_OpInf_ROM_Constraint
   api/model_reduction_Reduced_Dynamic_Objective
```

## Bayesian Inversion

:::{admonition} TODO

- `Bayesian_Inversion`
- `Bayesian_Inversion_Objective`
- `Inf_Dim_Prior_Model`
- `Likelihood_Model`
- `Mass_Matrix_Sqrt`
- `Prior_Model`
:::

## Optimal Experimental Design (OED)

:::{admonition} TODO

- `Forward_Operator_GSVD`
- `Linear_OED`
- `Misfit_Hessian_GEVP`
:::

## Linear Algebra Tools

:::{admonition} TODO

- `Matrix_Sqrt`
- `Randomized_GEVP`
- `Randomized_GSVD`
:::

## Model Discrepancy (HDSA)

:::{admonition} TODO

- `Elliptic_GSVD.m`
- `HDSA_Bayes_Posterior_Data.m`
- `HDSA_MD_Interface.m`
- `HDSA_MD_Interface_Elliptic_Prior.m`
- `HDSA_MD_Prior_Sampling.m`
- `HDSA_MD_Update.m`
- `HDSA_Sabl_MD_Interface_Elliptic_Prior.m`
- `Hessian_GEVP.m`
- `M_z_Sqrt.m`
:::

## Model Discrepancy Continuation

:::{admonition} TODO

- `HDSA_MD_Continuation_Interface`
- `HDSA_MD_Continuation_Update`
- `HDSA_Sabl_MD_Continuation_Interface`
:::

---

This API documentation was generated using [`sphinxcontrib_matlabdomain`](https://github.com/sphinx-contrib/matlabdomain).
