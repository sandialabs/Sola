import numpy as np
from abc import ABC, abstractmethod

class MD_u_Prior_Interface_Py():
    
    def __init__(self):
        pass
    
    # -------------------------------------------------------------------------------------------------------
    # Abstract methods which must be implemented

    @abstractmethod
    def Apply_W_u_Inverse_Py(self, u_in):
        pass

    @abstractmethod
    def Apply_M_u_Py(self, u_in):
        pass

    @abstractmethod
    def Apply_W_u_Plus_scalar_M_u_Inverse_Py(self, u_in, scalar):
        pass

    # -------------------------------------------------------------------------------------------------------
    # Optional methods which must be implemented to enable some analyes

    # Compute samples from a mean zero Gaussian with covariance W_u^{-1}
    def Sample_with_Covariance_W_u_Inverse_Py(self, num_samples):
        u_out = [];
        print('Sample_with_Covariance_W_u_Inverse_Py must be implemented to use sampling algorithms');
        return u_out

    # Compute samples from a mean zero Gaussian with covariance (W_u+scalar*M_u)^{-1}
    def Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse_Py(self, num_samples,scalar):
        u_out = [];
        disp('Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse_Py must be implemented to use sampling algorithms');
        return u_out

    # -------------------------------------------------------------------------------------------------------
    # Matlab-Python interoperabilty methods

    def Apply_W_u_Inverse(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_W_u_Inverse_Py(u_in)

    def Apply_M_u(self, u_in):
        u_in = np.array(u_in)
        return self.Apply_M_u_Py(u_in)

    def Apply_W_u_Plus_scalar_M_u_Inverse(self, u_in, scalar):
        u_in = np.array(u_in)
        return self.Apply_W_u_Plus_scalar_M_u_Inverse_Py(u_in,scalar)

    # Compute samples from a mean zero Gaussian with covariance W_u^{-1}
    def Sample_with_Covariance_W_u_Inverse(self, num_samples):
        return self.Sample_with_Covariance_W_u_Inverse_Py(num_samples)

    # Compute samples from a mean zero Gaussian with covariance (W_u+scalar*M_u)^{-1}
    def Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(self, num_samples, scalar):
        return self.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse_Py(num_samples,scalar)
