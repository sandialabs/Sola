# Define your functions here

# Modify this.
from fluid_flow_1d_hifi_eval import *

from start_server import *
input("Warning: This script executes functions over localhost. Are you sure you want to continue (press enter)?");
app.run(host='127.0.0.1', port=5001, debug=True)
