# This is the retrieved script seen by MATLAB
from retriever import *

# Set port to call
PORT = 5001; 

# Access Functions
state_solve = lambda *x, **kwargs: call_remote_function('state_solve', PORT, *x, **kwargs)
