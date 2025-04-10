import numpy as np
import sys

class Assemble_Mass_and_Stiffness():

    def __init__(self,m):
        m = int(m)
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
