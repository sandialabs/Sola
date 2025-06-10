import numpy as np
import sys
sys.path.append('../../../../src/python_adapter/model_discrepancy/interfaces/u_prior_interfaces')
from MD_Elliptic_u_Prior_Interface_Py import *

class MD_Elliptic_u_Prior_Interface_Python_Synthetic_Test(MD_Elliptic_u_Prior_Interface_Py):

    def __init__(self,m):
        self.m = m
        self.x = np.linspace(0.0,1.0,m)
        h = self.x[1] - self.x[0]

        self.M = 4.0*np.eye((m))
        self.M[0,0] = 2.0
        self.M[m-1,m-1] = 2.0
        for i in range(0,m-1):
            self.M[i,i+1] = 1.0
            self.M[i+1,i] = 1.0
        self.M = (h/6.0)*self.M

        self.S = 2.0*np.eye((m))
        self.S[0,0] = 1.0
        self.S[m-1,m-1] = 1.0
        for i in range(0,m-1):
            self.S[i,i+1] = -1.0
            self.S[i+1,i] = -1.0
        self.S = (1.0/h)*self.S

        self.E_u = (5.e-2) * self.S + self.M

    def Apply_E_u_Inverse_Py(self, u_in):
        u_out = np.linalg.solve(self.E_u,u_in)
        return u_out

    def Apply_E_u_Inverse_Transpose_Py(self, u_in):
        u_out = np.linalg.solve(np.transpose(self.E_u),u_in)
        return u_out

    def Apply_M_u_Py(self, u_in):
        u_out = self.M@u_in
        return u_out

