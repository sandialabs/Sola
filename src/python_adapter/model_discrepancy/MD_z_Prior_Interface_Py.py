import numpy as np
from abc import ABC, abstractmethod

class MD_z_Prior_Interface_Py():

    def __init__(self):
        pass

    # -------------------------------------------------------------------------------------------------------
    # Abstract methods which must be implemented

    @abstractmethod
    def Apply_W_z_Inverse_Py(self, z_in):
        pass

    # -------------------------------------------------------------------------------------------------------
    # Optional methods which must be implemented to enable some analyes

    # Factorize W_z^{-1} = F*F^T, function gives z_out = F*z_in
    # This function must be implemented to enable posterior update sampling
    def Apply_W_z_Inverse_Factor_Py(self, z_in):
        z_out = []
        print('Apply_W_z_Inverse_Factor_Py must be implemented to use sampling algorithms')
        return z_out

    # Apply W_z matrix
    # This function must be implemented to enable Hessian GEVP
    def Apply_W_z_Py(self, z_in):
        z_out = []
        print('Apply_W_z_Py must be implemented to use Hessian GEVP');
        return z_out

    # -------------------------------------------------------------------------------------------------------
    # Matlab-Python interoperabilty methods

    def Apply_W_z_Inverse(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_W_z_Inverse_Py(z_in)

    def Apply_W_z_Inverse_Factor(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_W_z_Inverse_Factor_Py(z_in)

    def Apply_W_z(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_W_z_Py(z_in)
