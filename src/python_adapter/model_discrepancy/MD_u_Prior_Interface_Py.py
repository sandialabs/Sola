import numpy as np
from abc import ABC, abstractmethod

class MD_u_Prior_Interface_Py():
    
    def __init__(self):
        
    
    # -------------------------------------------------------------------------------------------------------
    # Abstract methods which must be implemented

    @abstractmethod
    def Apply_W_u_Inverse_Py(self, u_in):
        pass

    @abstractmethod
    def Apply_W_d_Py(self, u_in):
        pass

    @abstractmethod
    def Apply_W_u_Plus_scalar_W_d_Inverse_Py(self, u_in, scalar):
        pass

    # -------------------------------------------------------------------------------------------------------
    # Optional methods which must be implemented to enable some analyes

    # Factorize W_u^{-1}=F*F^T, function gives u_out=F*u_in
    # This function must be implemented to enable posterior update sampling
    def Apply_W_u_Inverse_Factor_Py(self, u_in):
        u_out = [];
        print('Apply_W_u_Inverse_Factor_Py must be implemented to use sampling algorithms');
        return u_out    

    # Factorize (W_u+scalar*W_d)^{-1}=F*F^T, function gives u_out=F*u_in
    # This function must be implemented to enable posterior update sampling
    def Apply_W_u_Plus_scalar_W_d_Inverse_Factor_Py(self, u_in, scalar):
        u_out = []
        print('Apply_W_u_Plus_scalar_W_d_Inverse_Factor_Py must be implemented to use sampling algorithms')
        return u_out

    # -------------------------------------------------------------------------------------------------------
    # Matlab-Python interoperabilty methods

    def Apply_W_u_Inverse(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_W_u_Inverse_Py(u_in)

    def Apply_W_d(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_W_d_Py(u_in)

    def Apply_W_u_Plus_scalar_W_d_Inverse(self, u_in, scalar):
        u_in = np.array(u_in)
        return self.Apply_W_u_Plus_scalar_W_d_Inverse_Py(u_in,scalar)

    def Apply_W_u_Inverse_Factor(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_W_u_Inverse_Factor_Py(u_in)

    def Apply_W_u_Plus_scalar_W_d_Inverse_Factor(self, u_in, scalar):
        u_in = np.array(u_in)
        return self.Apply_W_u_Plus_scalar_W_d_Inverse_Factor_Py(u_in,scalar)
