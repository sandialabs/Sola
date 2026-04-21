This Python adapter is designed to implement Sola pure virtual functions using models in Python. Matlab should be the "driver" and the pure virtual functions should be implemented in a Python class. This will allow Sola to execute its outer loop algorithms where the model constraining the analysis is implemented in Python.

Note that the Matlab and Python versions must be compatible. See https://www.mathworks.com/support/requirements/python-compatibility.html for a list of compatible versions. Developments thus have have used Matlab 2022b and Python 3.9.16.
