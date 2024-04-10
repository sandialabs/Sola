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

    # Compute samples from a mean zero Gaussian with covariance W_z^{-1}
    def Sample_with_Covariance_W_z_Inverse_Py(this,num_samples)
        z_out = [];
        print('Sample_with_Covariance_W_z_Inverse must be implemented to use sampling algorithms');
        return z_out
    end

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

    % Compute samples from a mean zero Gaussian with covariance W_z^{-1}
    def Sample_with_Covariance_W_z_Inverse(this,num_samples):
        return self.Sample_with_Covariance_W_z_Inverse_Py(num_samples)

    def Apply_W_z(self, z_in):
        z_in = np.array(z_in)
        return self.Apply_W_z_Py(z_in)
