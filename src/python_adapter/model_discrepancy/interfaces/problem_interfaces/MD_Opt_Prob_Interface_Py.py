import numpy as np
from abc import ABC, abstractmethod

class MD_Opt_Prob_Interface_Py():

    def __init__(self):
        pass
    
    # -------------------------------------------------------------------------------------------------------
    # Abstract methods which must be implemented

    @abstractmethod
    def Apply_Solution_Operator_z_Jacobian_Transpose_Py(self, u_in, z):
        pass

    @abstractmethod
    def Apply_RS_Hessian_Py(self, z_in, z):
        pass

    @abstractmethod
    def Misfit_Gradient_Py(self, u, z):
        pass

    @abstractmethod
    def Apply_Misfit_Hessian_Py(self, u_in, u, z):
        pass

    # -------------------------------------------------------------------------------------------------------
    # Optional methods which must be implemented to enable some analyes

    def State_Solve_Py(self, z):
        u = []
        print('State_Solve_Py must be implemented to use continuation algorithm')
        return u

    def Apply_Solution_Operator_z_Jacobian_Py(self, z_in, z):
        u_out = []
        print('Apply_Solution_Operator_z_Jacobian_Py must be implemented to use continuation algorithm')
        return u_out

    # -------------------------------------------------------------------------------------------------------
    # Matlab-Python interoperabilty methods

    def Apply_Solution_Operator_z_Jacobian_Transpose(self, u_in, z):
        u_in = np.array(u_in)
        z = np.array(z)
        return self.Apply_Solution_Operator_z_Jacobian_Transpose_Py(u_in, z)

    def Apply_RS_Hessian(self, z_in, z):
        z_in = np.array(z_in)
        z = np.array(z)
        return self.Apply_RS_Hessian_Py(z_in, z)

    def Misfit_Gradient(self, u, z):
        u = np.array(u)
        z = np.array(z)
        return self.Misfit_Gradient_Py(u, z)

    def Apply_Misfit_Hessian(self, u_in, u, z):
        u_in = np.array(u_in)
        u = np.array(u)
        z = np.array(z)
        return self.Apply_Misfit_Hessian_Py(u_in, u, z)

    def State_Solve(self, z):
        z = np.array(z)
        return self.State_Solve_Py(z)

    def Apply_Solution_Operator_z_Jacobian(self, z_in, z):
        z_in = np.array(z_in)
        z = np.array(z)
        return self.Apply_Solution_Operator_z_Jacobian_Py(z_in, z)
