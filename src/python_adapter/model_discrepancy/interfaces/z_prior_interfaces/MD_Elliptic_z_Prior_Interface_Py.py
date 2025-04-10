import numpy as np
from abc import ABC, abstractmethod

class MD_Elliptic_z_Prior_Interface_Py():
    
    def __init__(self):
        pass
    
    # -------------------------------------------------------------------------------------------------------
    # Abstract methods which must be implemented

    @abstractmethod
    def Apply_E_z_Inverse_Py(self, z_in):
        pass

    @abstractmethod
    def Apply_E_z_Inverse_Transpose_Py(self, z_in):
        pass

    @abstractmethod
    def Apply_M_z_Py(self, z_in):
        pass

    # -------------------------------------------------------------------------------------------------------
    # Optional methods which must be implemented to enable some analyes

    def Apply_E_z_Py(self, z_in):
        z_out = []
        print('Apply_E_z_Py must be implemented to use Hessian GEVP')
        return z_out

    def Apply_E_z_Transpose_Py(self, z_in):
        z_out = []
        print('Apply_E_z_Transpose_Py must be implemented to use Hessian GEVP')
        return z_out

    def Apply_M_z_Inverse_Py(self, z_in):
        z_out = []
        print('Apply_M_z_Inverse_Py must be implemented to use Hessian GEVP')
        return z_out

    # Compute samples from a mean zero Gaussian with covariance W_z^{-1}
    def Sample_with_Covariance_W_z_Inverse_Py(self, num_samples):
        z_out = [];
        print('Sample_with_Covariance_W_z_Inverse must be implemented to use sampling algorithms');
        return z_out

    # -------------------------------------------------------------------------------------------------------
    # Matlab-Python interoperabilty methods

    def Apply_E_z_Inverse(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_E_z_Inverse_Py(z_in)

    def Apply_E_z_Inverse_Transpose(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_E_z_Inverse_Transpose_Py(z_in)

    def Apply_M_z(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_M_z_Py(z_in)

    def Apply_E_z(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_E_z_Py(z_in)

    def Apply_E_z_Transpose(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_E_z_Transpose_Py(z_in)

    def Apply_M_z_Inverse(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_M_z_Inverse_Py(z_in)

    # Compute samples from a mean zero Gaussian with covariance W_z^{-1}
    def Sample_with_Covariance_W_z_Inverse(self,num_samples):
        return self.Sample_with_Covariance_W_z_Inverse_Py(num_samples)
