# Server Imports
from flask import Flask, request, jsonify
import inspect
import numpy as np
import base64
import json
from __main__ import *

app = Flask(__name__)

# Helper functions to serialize and deserialize NumPy arrays
def serialize_array(arr):
    encoded = base64.b64encode(arr.tobytes()).decode('utf-8')
    return {
        'data': encoded,
        'shape': arr.shape,
        'dtype': str(arr.dtype)
    }

def deserialize_array(serialized):
    decoded = base64.b64decode(serialized['data'])
    return np.frombuffer(decoded, dtype=serialized['dtype']).reshape(serialized['shape'])

# Function to dynamically call any function
@app.route('/call_function', methods=['POST'])
def call_function():
    data = request.json
    # print(f"Received request: {data}")
    function_name = data['function_name']
    params = data['params']
    kwargs = data.get('kwargs', {})
    
    # Deserialize any arrays in the parameters
    for i, param in enumerate(params):
        if isinstance(param, dict) and param.get('type') == 'ndarray':
            params[i] = deserialize_array(param['data'])
    
    # Deserialize any arrays in the kwargs
    for key, value in kwargs.items():
        if isinstance(value, dict) and value.get('type') == 'ndarray':
            kwargs[key] = deserialize_array(value['data'])
    
    # Get the function from the global scope
    func = globals().get(function_name)
    
    if not func:
        return jsonify(error=f"Function {function_name} not found"), 404
    
    # Check if the function is callable
    if not callable(func):
        return jsonify(error=f"{function_name} is not callable"), 400
    
    # Get the function signature
    sig = inspect.signature(func)
    
    # Bind the parameters to the function signature
    try:
        bound_args = sig.bind(*params, **kwargs)
    except TypeError as e:
        return jsonify(error=str(e)), 400
    
    # Call the function with the bound arguments
    result = func(*bound_args.args, **bound_args.kwargs)
    
    # Serialize the result if it's a NumPy array
    if isinstance(result, np.ndarray):
        result = {
            'type': 'ndarray',
            'data': serialize_array(result)
        }
    
    return jsonify(result=result)

# Function to get the value of a variable
@app.route('/get_variable', methods=['GET'])
def get_variable():
    variable_name = request.args.get('variable_name')
    # print(f"Received request for variable: {variable_name}")
    
    # Get the variable from the global scope
    variable = globals().get(variable_name)
    
    if variable is None:
        return jsonify(error=f"Variable {variable_name} not found"), 404
    
    # Serialize the variable if it's a NumPy array
    if isinstance(variable, np.ndarray):
        variable = {
            'type': 'ndarray',
            'data': serialize_array(variable)
        }
    
    return jsonify(variable=variable)