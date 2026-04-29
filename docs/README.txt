The document "Optimization.pdf" describes the optimization problems that are solved in Sola and provides a guide for understanding the code contained in
-sola/src/optimization
-sola/src/automatic_differentation
-sola/src/optimization_under_uncertainty

The remaining sola modules contained in 
-sola/src/bayesian_inversion
-sola/src/linear_algebra_tools
-sola/src/model_discrepancy
-sola/src/optimal_experimental_design
-sola/src/pseudo_time_continuation
have not yet been documented.

The code in sola/src/python_adapter is designed to implement Sola pure virtual functions using models in Python. Matlab should be the "driver" and the pure virtual functions should be implemented in a Python class. This will allow Sola to execute its outer loop algorithms where the model constraining the analysis is implemented in Python.

Note that the Matlab and Python versions must be compatible. See https://www.mathworks.com/support/requirements/python-compatibility.html for a list of compatible versions. Developments thus have used Matlab 2022b and Python 3.9.16.

Also, note that the python_adapter may not be updated to support every functionality in Sola. Additional python interfaces can be implemented following the design of the code in sola/src/python_adapter.