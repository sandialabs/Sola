import numpy as np
import sys
sys.path.append('../../../../src/python_adapter/model_discrepancy/')
from MD_Data_Interface_Py import *

class MD_Data_Interface_Python_Synthetic_Test(MD_Data_Interface_Py):

    def __init__(self,m):
        self.m = m
        self.x = np.linspace(0.0,1.0,m)
    
    def Load_Optimal_u_Py(self):
        return (1.0+self.x)**3

    def Load_Optimal_z_Py(self):
        return 1.0+self.x

    def Load_Z_Data_Py(self):
        Z = np.zeros((self.m,2))
        Z[:,0] = 1.0 + self.x
        Z[:,1] = self.x + self.x**2
        return Z

    def Load_d_Data_Py(self):
        Z = self.Load_Z_Data()
        D = 0.2*(Z**3)
        return D
