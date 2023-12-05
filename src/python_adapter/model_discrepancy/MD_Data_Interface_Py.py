import numpy as np
from abc import ABC, abstractmethod

class MD_Data_Interface_Py():

    def __init__(self):
        pass
    
    # -------------------------------------------------------------------------------------------------------
    # Abstract methods which must be implemented
    
    @abstractmethod
    def Load_Optimal_u_Py(self):
        pass

    @abstractmethod
    def Load_Optimal_z_Py(self):
        pass

    # -------------------------------------------------------------------------------------------------------
    # Optional methods which must be implemented to enable some analyes

    def Load_Z_Data_Py(self):
        Z = []
        print('Load_Z_Data_Py must be implemented for any analyses except for optimal experimental design')
        return Z

    def Load_d_Data_Py(self):
        D = []
        print('Load_d_Data_Py must be implemented for any analyses except for optimal experimental design')
        return D

    # -------------------------------------------------------------------------------------------------------\
    # Matlab-Python interoperabilty methods

    def Load_Optimal_u(self):
        return self.Load_Optimal_u_Py()

    def Load_Optimal_z(self):
        return self.Load_Optimal_z_Py()

    def Load_Z_Data(self):
        Z = self.Load_Z_Data_Py()
        return Z

    def Load_d_Data(self):
        D = self.Load_d_Data_Py()
        return D
