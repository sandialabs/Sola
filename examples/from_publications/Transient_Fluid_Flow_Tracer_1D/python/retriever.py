import requests
import numpy as np
import base64
import json
from array import array


def convert_ifvec_numpy(input):
    if isinstance(input, np.ndarray) or isinstance(input, memoryview) or isinstance(input, array):
        return np.array(input);
    else:
        return input;

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

def call_remote_function(function_name, port, *params, **kwargs):
    url = f'http://localhost:{port}/call_function'  # Ensure this URL is correct
    serialized_params = []
    
    for param in params:
        param = convert_ifvec_numpy(param)
        if isinstance(param, tuple):
            print(param)
        if isinstance(param, np.ndarray):
            serialized_params.append({
                'type': 'ndarray',
                'data': serialize_array(param)
            })
        else:
            serialized_params.append(param)
    
    serialized_kwargs = {}
    for key, value in kwargs.items():
        if isinstance(value, tuple):
            print(key, value)
        value = convert_ifvec_numpy(value)
        if isinstance(value, np.ndarray):
            serialized_kwargs[key] = {
                'type': 'ndarray',
                'data': serialize_array(value)
            }
        else:
            serialized_kwargs[key] = value
    
    data = {
        'function_name': function_name,
        'params': serialized_params,
        'kwargs': serialized_kwargs
    }
    # print(f"Sending request to {url} with data: {data}")
    response = requests.post(url, json=data)
    
    # print(f"Received response: {response.text}")
    if response.status_code == 200:
        result = response.json()['result']
        if isinstance(result, dict) and result.get('type') == 'ndarray':
            return deserialize_array(result['data'])
        return result
    else:
        raise Exception(response.json().get('error', 'Unknown error'))

def get_remote_variable(variable_name, port):
    url = f'http://localhost:{port}/get_variable?variable_name={variable_name}'  # Ensure this URL is correct
    # print(f"Sending request to {url}")
    response = requests.get(url)
    
    # print(f"Received response: {response.text}")
    if response.status_code == 200:
        variable = response.json()['variable']
        if isinstance(variable, dict) and variable.get('type') == 'ndarray':
            return deserialize_array(variable['data'])
        return variable
    else:
        raise Exception(response.json().get('error', 'Unknown error'))