import numpy as np
from abc import ABC, abstractmethod

class MD_u_Hyperparameter_Interface_Py():

    def __init__(self):
        pass
    
    # -------------------------------------------------------------------------------------------------------
    # Optional methods which must be implemented to enable some analyes

    def Load_Spatial_Node_Data_Py(self):
        spatial_nodes = []
        print('Load_Spatial_Node_Data is required for hyperparameter algorithm-based initialization')
        return spatial_nodes

    def Load_Time_Node_Data_Py(self):
        time_nodes = []
        print('Load_Time_Node_Data is required for hyperparameter algorithm-based initialization')
        return time_nodes

    # -------------------------------------------------------------------------------------------------------\
    # Matlab-Python interoperabilty methods

    def Load_Spatial_Node_Data(self):
        spatial_nodes = self.Load_Spatial_Node_Data_Py()
        return spatial_nodes

    def Load_Time_Node_Data(self):
        time_nodes = self.Load_Time_Node_Data_Py()
        return time_nodes
