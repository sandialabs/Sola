import numpy as np
import sys
sys.path.append('../../../../src/python_adapter/model_discrepancy/')
from MD_Opt_Prob_Interface_Py import *

class MD_Opt_Prob_Interface_Python_Synthetic_Test(MD_Opt_Prob_Interface_Py):

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

    def Apply_Solution_Operator_z_Jacobian_Transpose_Py(self, u_in, z):
        z_out = 3.0*np.diag(z**2)@u_in
        return z_out

    def Apply_RS_Hessian_Py(self, z_in, z):
        tmp1 = 3.0*np.diag(z**2)@z_in
        tmp2 = self.M@tmp1
        z_out = 3.0*np.diag(z**2)@tmp2
        return z_out

    def Misfit_Gradient_Py(self, u, z):
        tmp = u - (1.0+self.x)**3
        grad_u = self.M@tmp
        return grad_u

    def Apply_Misfit_Hessian_Py(self, u_in, u, z):
        u_out = self.M@u_in
        return u_out
