import numpy as np
import sys
sys.path.append('../../../../src/python_adapter/model_discrepancy/interfaces/z_prior_interfaces')
from MD_z_Prior_Interface_Py import *

class MD_z_Prior_Interface_Python_Synthetic_Test(MD_z_Prior_Interface_Py):

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

        self.E_z = (1.e2)* ( (1.e-2)*self.S + self.M )

        self.W_z = np.transpose(self.E_z)@np.linalg.solve(self.M,self.E_z)

    def Apply_M_z_Py(self, z_in):
        z_out = self.M@z_in
        return z_out

    def Apply_W_z_Inverse_Py(self, z_in):
        z_out = np.linalg.solve(self.W_z,z_in)
        return z_out

    def Apply_W_z_Py(self, z_in):
        z_out = self.W_z@z_in
        return z_out

    # Compute samples from a mean zero Gaussian with covariance W_z^{-1}
    def Sample_with_Covariance_W_z_Inverse(self, num_samples):
        num_samples = int(num_samples)
        Omega = np.random.standard_normal((self.m,num_samples)) 
        L = np.linalg.cholesky(self.W_z)
        z_out = np.transpose(L)@Omega
        return z_out
