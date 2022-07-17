"""
© 2022. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.
"""

import sys
sys.path.append("../")
from utils import *
from create_telecloning_circuits import *
from qiskit import *
import copy
from qiskit.transpiler import CouplingMap
import numpy as np

OPT_LEVEL = 3

device_name = "ibmq_montreal"

#subg = [0, 1, 4, 7]
#subg = [5, 8, 11, 14]
#subg = [10, 12, 15, 18]
subg = [16, 19, 22, 25]

n = len(subg)

IBMQ.load_account()
provider = IBMQ.get_provider(hub='', group='', project='')
backend = provider.get_backend(device_name)
coupling_map = CouplingMap(backend.configuration().coupling_map)

N_shots = 30000
N_Angle_slices = 17

#Create a name for the qubit set used when storing results
subgraph_name = "_"
for a in subg:
	subgraph_name += str(a)+"_"
subgraph_name += "_"

print("Initial layout:", subgraph_name)

def get_transpiled_circuits(angle_ry, angle_rz):
	qc_tc, clone_indices, _, post_select_indices = construct_PCC_postselect_00(angle_ry, angle_rz)
	tomography_circuits_00 = parallel_qubit_state_tomography_general(clone_indices, copy.deepcopy(qc_tc), angle_ry, angle_rz)
	tomography_circuits_00 = add_measurements_for_postselect_circuits_for_tomography_circuits(copy.deepcopy(tomography_circuits_00), post_select_indices)
	transpiled_00 = transpile(copy.deepcopy(tomography_circuits_00), coupling_map=coupling_map, optimization_level=OPT_LEVEL, basis_gates=["x", "sx", "cx", "rz"], initial_layout=subg)
	
	qc_tc, clone_indices, _, post_select_indices = construct_PCC_postselect_11(angle_ry, angle_rz)
	tomography_circuits_11 = parallel_qubit_state_tomography_general(clone_indices, copy.deepcopy(qc_tc), angle_ry, angle_rz)
	tomography_circuits_11 = add_measurements_for_postselect_circuits_for_tomography_circuits(copy.deepcopy(tomography_circuits_11), post_select_indices)
	transpiled_11 = transpile(copy.deepcopy(tomography_circuits_11), coupling_map=coupling_map, optimization_level=OPT_LEVEL, basis_gates=["x", "sx", "cx", "rz"], initial_layout=subg)
	
	qc_tc, clone_indices, _, post_select_indices = construct_PCC_postselect_10(angle_ry, angle_rz)
	tomography_circuits_10 = parallel_qubit_state_tomography_general(clone_indices, copy.deepcopy(qc_tc), angle_ry, angle_rz)
	tomography_circuits_10 = add_measurements_for_postselect_circuits_for_tomography_circuits(copy.deepcopy(tomography_circuits_10), post_select_indices)
	transpiled_10 = transpile(copy.deepcopy(tomography_circuits_10), coupling_map=coupling_map, optimization_level=OPT_LEVEL, basis_gates=["x", "sx", "cx", "rz"], initial_layout=subg)
	
	qc_tc, clone_indices, _, post_select_indices = construct_PCC_postselect_01(angle_ry, angle_rz)
	tomography_circuits_01 = parallel_qubit_state_tomography_general(clone_indices, copy.deepcopy(qc_tc), angle_ry, angle_rz)
	tomography_circuits_01 = add_measurements_for_postselect_circuits_for_tomography_circuits(copy.deepcopy(tomography_circuits_01), post_select_indices)
	transpiled_01 = transpile(copy.deepcopy(tomography_circuits_01), coupling_map=coupling_map, optimization_level=OPT_LEVEL, basis_gates=["x", "sx", "cx", "rz"], initial_layout=subg)
	
	return transpiled_00, transpiled_01, transpiled_10, transpiled_11

angles_ry = np.linspace(0, math.pi, num=N_Angle_slices)
angles_rz = np.linspace(0, 2*math.pi, num=N_Angle_slices)

compiled_circuits_00 = []
compiled_circuits_01 = []
compiled_circuits_10 = []
compiled_circuits_11 = []

for angle_ry in angles_ry:
	for angle_rz in angles_rz:
		transpiled_00, transpiled_01, transpiled_10, transpiled_11 = get_transpiled_circuits(float(angle_ry), float(angle_rz))
		compiled_circuits_00 += transpiled_00
		compiled_circuits_01 += transpiled_01
		compiled_circuits_10 += transpiled_10
		compiled_circuits_11 += transpiled_11

print("total circuits before error mitigation =", len(compiled_circuits_11))
qubits_to_measure = extract_qubit_indices_which_are_measured(compiled_circuits_11)

measurement_err_mitigation_circuits = generate_transpiled_measurement_error_mitigation_circuits(qubits_to_measure, coupling_map)

compiled_circuits_11 = measurement_err_mitigation_circuits+compiled_circuits_11
compiled_circuits_00 = measurement_err_mitigation_circuits+compiled_circuits_00
compiled_circuits_01 = measurement_err_mitigation_circuits+compiled_circuits_01
compiled_circuits_10 = measurement_err_mitigation_circuits+compiled_circuits_10


job = backend.run(compiled_circuits_00, shots=N_shots)
id = job.job_id()
print(id)
job_filename = "PCC_"+device_name+subgraph_name+"_"+str(N_shots)+"_"+str(N_Angle_slices)+"_00_OPT_"+str(OPT_LEVEL)+".txt"
file = open("job_ids_postselect/"+job_filename, "w")
file.write(str(id))
file.close()

job = backend.run(compiled_circuits_01, shots=N_shots)
id = job.job_id()
print(id)
job_filename = "PCC_"+device_name+subgraph_name+"_"+str(N_shots)+"_"+str(N_Angle_slices)+"_01_OPT_"+str(OPT_LEVEL)+".txt"
file = open("job_ids_postselect/"+job_filename, "w")
file.write(str(id))
file.close()

job = backend.run(compiled_circuits_10, shots=N_shots)
id = job.job_id()
print(id)
job_filename = "PCC_"+device_name+subgraph_name+"_"+str(N_shots)+"_"+str(N_Angle_slices)+"_10_OPT_"+str(OPT_LEVEL)+".txt"
file = open("job_ids_postselect/"+job_filename, "w")
file.write(str(id))
file.close()

job = backend.run(compiled_circuits_11, shots=N_shots)
id = job.job_id()
print(id)
job_filename = "PCC_"+device_name+subgraph_name+"_"+str(N_shots)+"_"+str(N_Angle_slices)+"_11_OPT_"+str(OPT_LEVEL)+".txt"
file = open("job_ids_postselect/"+job_filename, "w")
file.write(str(id))
file.close()

