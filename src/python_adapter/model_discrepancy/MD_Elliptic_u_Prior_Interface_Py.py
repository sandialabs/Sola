import numpy as np
from abc import ABC, abstractmethod

class MD_Elliptic_u_Prior_Interface_Py():
    
    def __init__(self):
        pass
    
    # -------------------------------------------------------------------------------------------------------
    # Abstract methods which must be implemented

    @abstractmethod
    def Apply_E_u_Inverse_Py(self, u_in):
        pass

    @abstractmethod
    def Apply_E_u_Inverse_Transpose_Py(self, u_in):
        pass

    @abstractmethod
    def Apply_M_u_Py(self, u_in):
        pass

    @abstractmethod
    def Apply_M_u_Inverse_Py(self, u_in):
        pass

    @abstractmethod
    def Apply_E_d_Py(self, u_in):
        pass

    @abstractmethod
    def Apply_E_d_Transpose_Py(self, u_in):
        pass

    # -------------------------------------------------------------------------------------------------------
    # Matlab-Python interoperabilty methods

    def Apply_E_u_Inverse(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_E_u_Inverse_Py(u_in)

    def Apply_E_u_Inverse_Transpose(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_E_u_Inverse_Transpose_Py(u_in)

    def Apply_M_u(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_M_u_Py(u_in)

    def Apply_M_u_Inverse(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_M_u_Inverse_Py(u_in)

    def Apply_E_d(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_E_d_Py(u_in)

    def Apply_E_d_Transpose(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_E_d_Transpose_Py(u_in)
