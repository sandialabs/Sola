import numpy as np
from typing import Literal
from array import array
from dolfin import *
from dolfin_adjoint import *
from pyadjoint.enlisting import Enlist

# Handy Helpers to pause verbosity and annotation in one go:
class stop_annotate_verbose(object):
    def __enter__(self): pause_annotation(); set_log_active(False)
    def __exit__(self, *args): continue_annotation(); set_log_active(True)
class stop_verbose(object):
    def __enter__(self): set_log_active(False)
    def __exit__(self, *args): set_log_active(True)


def compute_jacobian_action(J, m, m_dot):
    tape = get_working_tape()
    tape.reset_tlm_values()
    m_dot = Enlist(m_dot); m = Enlist(m)

    for i in range(len(m_dot)):
        m[i].tlm_value = m_dot[i]

    with stop_annotating():
        with tape.marked_nodes(m):
            tape.evaluate_tlm()
            if isinstance(J.block_variable.tlm_value, np.ndarray):
                raise NotImplementedError("J.block_variable.tlm_value returns array.")
            else:
                Jmdots = J.block_variable.tlm_value
    
    return Jmdots



# Fix Boundary Measure
class MeasureDiff(object):
    def __init__(self, measure_1, measure_2):
        self.measure_1 = measure_1
        self.measure_2 = measure_2
    def __rmul__(self, other):
        return other * self.measure_1 - other * self.measure_2
    def __str__(self):
        return "{\n    " + str(self.measure_1) + "\n  - " + str(self.measure_2) + "\n}"
setattr(Measure, "__sub__", lambda self, other: MeasureDiff(self, other))
# Corrects an issue with FEniCS that incorrectly assigns fenics.ds
class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0)
class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 1)

def generate_ds(int_mesh):
    bdry = MeshFunction("size_t", int_mesh, 0, 0)
    Gamma_left = LeftBoundary(); Gamma_left.mark(bdry, 0)
    Gamma_right = RightBoundary(); Gamma_right.mark(bdry, 1)
    ds_full = Measure("ds", domain=int_mesh, subdomain_data=bdry)
    ds = ds_full(1) - ds_full(0)
    return ds

# Fenics conversion helper functions:
def fenics_convert(fenics_input, return_type: Literal["vertex", "vector", "petsc", "function"], fun_space=None):
    if return_type == "function":
        return convert_to_function(fenics_input, fun_space)
    elif return_type == "vector":
        return convert_to_numpy_vector(fenics_input, fun_space)
    elif return_type == "petsc":
        return convert_to_petsc_vector(fenics_input, fun_space)
    elif return_type == "vertex":
        return convert_to_vertex(fenics_input, fun_space)
    else:
        raise ValueError("Unknown return_type received (must be one of vector/petsc/function).")


def convert_to_function(fenics_input, fun_space=None):
    if isinstance(fenics_input, Function):
        return fenics_input;
    elif isinstance(fenics_input, Expression):
        return project(fenics_input, fun_space);
    elif isinstance(fenics_input, str):
        return project(Expression(fenics_input, degree=1), fun_space);
    elif isinstance(fenics_input, PETScVector):
        output = Function(fun_space);
        output.vector()[:] = fenics_input[:];
        return output
    elif isinstance(fenics_input, np.ndarray):
        output = Function(fun_space);
        output.vector()[:] = fenics_input.flatten();
        return output
    elif isinstance(fenics_input, array):
        output = Function(fun_space);
        # output.vector()[:] = np.array(fenics_input).flatten();
        output.vector()[:] = fenics_input;
        return output
    elif isinstance(fenics_input, memoryview):
        output = Function(fun_space);
        output.vector()[:] = np.array(fenics_input).flatten();
        return output
    else:
        raise TypeError(f"Input can only be of type(s): str, fenics.Expression, fenics.Function, numpy.ndarray, array.array, memoryview \n Recieved input of type: {type(fenics_input)}")

def convert_to_numpy_vector(fenics_input, fun_space=None):
    if isinstance(fenics_input, np.ndarray) or isinstance(fenics_input, memoryview) or isinstance(fenics_input, array):
        return np.array(fenics_input);
    elif isinstance(fenics_input, PETScVector):
        return fenics_input[:]
    fenics_function = convert_to_function(fenics_input, fun_space);
    return fenics_function.vector()[:]

def convert_to_petsc_vector(fenics_input, fun_space=None):
    if isinstance(fenics_input, PETScVector):
        return fenics_input
    fenics_function = convert_to_function(fenics_input, fun_space);
    return fenics_function.vector()

def convert_to_vertex(fenics_input, fun_space=None):
    fenics_function = convert_to_function(fenics_input, fun_space);
    return fenics_function.compute_vertex_values()
