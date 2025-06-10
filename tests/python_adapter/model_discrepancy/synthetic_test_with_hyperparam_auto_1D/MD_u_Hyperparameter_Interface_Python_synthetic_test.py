import numpy as np
import sys
sys.path.append('../../../../src/python_adapter/model_discrepancy/interfaces/hyperparameter_interfaces')
from MD_u_Hyperparameter_Interface_Py import *

class MD_u_Hyperparameter_Interface_Python_synthetic_test(MD_u_Hyperparameter_Interface_Py):

    def __init__(self,m):
        m = int(m)
        self.x = np.linspace(0.0,1.0,m)
    
    def Load_Spatial_Node_Data_Py(self):
        spatial_nodes = []
        spatial_nodes.append(1.0+self.x)
        return spatial_nodes
