# Define your functions here

# Modify this.
from fluid_flow_1d_hifi_eval import *

from start_server import *
app.run(host='0.0.0.0', port=5001, debug=True)

"""
conda activate FenicsEnvNew
python /Users/mmadhav/GitRepos/sabl/examples/Transient_Fluid_Flow_Tracer_1D/python/serverside/run_server_hifi.py
python /Users/mmadhav/GitRepos/sabl/examples/Transient_Fluid_Flow_Tracer_1D/python/serverside/run_server_lofi_linear.py
""";