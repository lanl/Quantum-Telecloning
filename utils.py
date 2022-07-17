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

from qiskit import QuantumCircuit, transpile
import copy
import re
from qiskit.ignis.mitigation.measurement import complete_meas_cal
from qiskit.ignis.verification import marginal_counts
import seaborn as sns
import matplotlib.pyplot as plt

def count_postselection(counts, post_select_dictionary):
	total = sum(list(counts.values()))
	single_qubit_counts = {}
	valid_counts = 0
	for bitstring in counts:
		bitstring_reversed = bitstring[::-1]
		not_valid = False
		for qubit_index in post_select_dictionary:
			target_value = post_select_dictionary[qubit_index]
			if bitstring_reversed[qubit_index] != str(target_value):
				not_valid = True
		if not_valid == True:
			continue
		valid_counts += counts[bitstring]
	return valid_counts / float(total)
def decrement_keys_by_n(dictionary, n):
	"""
	decrements the largest key value by n
	"""
	out = {}
	for k in dictionary:
		if k == max(list(dictionary.keys())):
			out[k-n] = dictionary[k]
		else:
			out[k] = dictionary[k]
	return out

def parallel_qubit_state_tomography_general(qubit_indices, circuit, angle_ry, angle_rz):
	pauli_y_basis = circuit.copy("Y_"+str(angle_ry)+"_"+str(angle_rz))
	pauli_y_basis.barrier()
	for idx in qubit_indices:
		pauli_y_basis.sdg(idx)
		pauli_y_basis.h(idx)
	pauli_y_basis.barrier()
	for idx in qubit_indices:
		pauli_y_basis.measure(idx, idx)
	
	pauli_x_basis = circuit.copy("X_"+str(angle_ry)+"_"+str(angle_rz))
	pauli_x_basis.barrier()
	for idx in qubit_indices:
		pauli_x_basis.h(idx)
	pauli_x_basis.barrier()
	for idx in qubit_indices:
		pauli_x_basis.measure(idx, idx)
	
	pauli_z_basis = circuit.copy("Z_"+str(angle_ry)+"_"+str(angle_rz))
	pauli_z_basis.barrier()
	for idx in qubit_indices:
		pauli_z_basis.measure(idx, idx)
	
	return [pauli_y_basis, pauli_x_basis, pauli_z_basis]
def decrement_keys_by_1(dictionary):
	"""
	decrements the largest key value by 1
	"""
	out = {}
	for k in dictionary:
		if k == max(list(dictionary.keys())):
			out[k-1] = dictionary[k]
		else:
			out[k] = dictionary[k]
	return out
def save_matrix(matrix, state_labels, filename):
	sns.heatmap(matrix, yticklabels=state_labels, xticklabels=state_labels, cmap="gist_gray_r")
	plt.xlabel("Prepared state", fontsize=15)
	plt.ylabel("Measured state", fontsize=15)
	fig = plt.gcf()
	fig.set_size_inches(8, 6)
	plt.tight_layout()
	#plt.savefig("measurement_error_calibration_matrices/deferred_measurement_"+TYPE+"_"+device_name+subgraph_name+"_"+str(N_shots)+"_"+str(N_Angle_slices)+".pdf")
	plt.savefig(filename)
	plt.close(fig)
def pull_only_measured_indices(counts, indices):
	return marginal_counts(counts, indices)
def extract_qubit_indices_which_are_measured(compiled_circuits):
	openqasm = compiled_circuits[0].qasm()
	measurement_qubits = []
	for line in openqasm.splitlines():
		if "measure" in line:
			res = re.findall(r"\[\s*\+?(-?\d+)\s*\]", line)
			measurement_qubits.append(int(res[0]))
	return measurement_qubits
def generate_transpiled_measurement_error_mitigation_circuits(qubits_to_measure, coupling_map):
	assert type(qubits_to_measure) is list
	meas_calibs, state_labels = complete_meas_cal(qubit_list=[i for i in range(len(qubits_to_measure))], circlabel='mcal')
	transpiled = transpile(copy.deepcopy(meas_calibs), coupling_map=coupling_map, optimization_level=0, basis_gates=["x", "sx", "cx", "rz"], initial_layout=qubits_to_measure)
	return transpiled
def add_measurements_for_postselect_circuits(qc, post_select_indices):
	qc.barrier()
	for i in list(post_select_indices.keys()):
		qc.measure(i, i)
	return qc
def add_measurements_for_postselect_circuits_for_tomography_circuits(circuits, post_select_indices):
	out = []
	for qc in circuits:
		for i in list(post_select_indices.keys()):
			qc.measure(i, i)
		out.append(qc)
	return out
def state_tomography_clone1(qubit_indices, circuit):
	pauli_y_basis = circuit.copy("Y")
	pauli_y_basis.barrier()
	pauli_y_basis.sdg(qubit_indices[0])
	pauli_y_basis.h(qubit_indices[0])
	pauli_y_basis.measure(qubit_indices[0], qubit_indices[0])
	
	pauli_x_basis = circuit.copy("X")
	pauli_x_basis.barrier()
	pauli_x_basis.h(qubit_indices[0])
	pauli_x_basis.measure(qubit_indices[0], qubit_indices[0])
	
	pauli_z_basis = circuit.copy("Z")
	pauli_z_basis.barrier()
	pauli_z_basis.measure(qubit_indices[0], qubit_indices[0])
	
	return [pauli_y_basis, pauli_x_basis, pauli_z_basis]
def state_tomography_clone2(qubit_indices, circuit):
	pauli_y_basis = circuit.copy("Y")
	pauli_y_basis.barrier()
	pauli_y_basis.sdg(qubit_indices[1])
	pauli_y_basis.h(qubit_indices[1])
	pauli_y_basis.measure(qubit_indices[1], qubit_indices[1])
	
	pauli_x_basis = circuit.copy("X")
	pauli_x_basis.barrier()
	pauli_x_basis.h(qubit_indices[1])
	pauli_x_basis.measure(qubit_indices[1], qubit_indices[1])
	
	pauli_z_basis = circuit.copy("Z")
	pauli_z_basis.barrier()
	pauli_z_basis.measure(qubit_indices[1], qubit_indices[1])
	
	return [pauli_y_basis, pauli_x_basis, pauli_z_basis]
def get_single_qubit_measurements_from_parallel_results(counts, index):
	single_qubit_counts = {}
	ones = 0
	zeros = 0
	for bitstring in counts:
		bitstring_reversed = bitstring[::-1]
		if bitstring_reversed[index] == "1":
			ones += counts[bitstring]
		elif bitstring_reversed[index] == "0":
			zeros += counts[bitstring]
	single_qubit_counts = {"0": zeros, "1": ones}
	return single_qubit_counts
def get_single_qubit_measurements_from_parallel_results_with_postselection(counts, index, post_select_dictionary):
	"""
	index is the clone index
	"""
	total = sum(list(counts.values()))
	single_qubit_counts = {}
	ones = 0
	zeros = 0
	for bitstring in counts:
		bitstring_reversed = bitstring[::-1]
		not_valid = False
		for qubit_index in post_select_dictionary:
			target_value = post_select_dictionary[qubit_index]
			if bitstring_reversed[qubit_index] != str(target_value):
				not_valid = True
		if not_valid == True:
			continue
		if bitstring_reversed[index] == "1":
			ones += counts[bitstring]
		elif bitstring_reversed[index] == "0":
			zeros += counts[bitstring]
	#print("Proportion kept = ", (zeros+ones) / float(total))
	single_qubit_counts = {"0": zeros, "1": ones}
	return single_qubit_counts
def append_parallel_qubit_tomography_to_Honeywell_QASM(qubit_indices, circuit):
	pauli_y_basis = copy.deepcopy(circuit)
	pauli_x_basis = copy.deepcopy(circuit)
	pauli_z_basis = copy.deepcopy(circuit)
	
	for qubit in qubit_indices:
		pauli_y_basis += "sdg q["+str(qubit)+"];\n"
		pauli_y_basis += "h q["+str(qubit)+"];\n"
		pauli_y_basis += "measure q["+str(qubit)+"] -> c["+str(qubit)+"];\n"
	
	for qubit in qubit_indices:
		pauli_x_basis += "h q["+str(qubit)+"];\n"
		pauli_x_basis += "measure q["+str(qubit)+"] -> c["+str(qubit)+"];\n"
	
	for qubit in qubit_indices:
		pauli_z_basis += "measure q["+str(qubit)+"] -> c["+str(qubit)+"];\n"
	
	return pauli_y_basis, pauli_x_basis, pauli_z_basis
def remove_measurements_from_qasm_string(QASM):
	split_str = QASM.split("\n")
	original = []
	split_str.reverse()
	barrier_removed_count = 0
	measure_removed_count = 0
	for line in split_str:
		if "barrier" in line:
			if barrier_removed_count < 1:
				barrier_removed_count += 1
				continue
		if "measure" in line:
			if measure_removed_count < 2:
				measure_removed_count += 1
				continue
		original.append(line)
	original.reverse()
	output_qasm = ""
	for a in original:
		output_qasm += a+"\n"
	return output_qasm
def get_state_tomography_clone1_for_post_selection(qubit_indices, circuit):
	qasm_reprentation = circuit.qasm()
	qasm_reprentation = remove_measurements_from_qasm_string(qasm_reprentation)
	circuit = QuantumCircuit.from_qasm_str(qasm_reprentation)
	
	pauli_y_basis = circuit.copy("Y")
	pauli_y_basis.barrier()
	pauli_y_basis.sdg(qubit_indices[0])
	pauli_y_basis.h(qubit_indices[0])
	pauli_y_basis.measure(qubit_indices[0], qubit_indices[0])
	
	pauli_x_basis = circuit.copy("X")
	pauli_x_basis.barrier()
	pauli_x_basis.h(qubit_indices[0])
	pauli_x_basis.measure(qubit_indices[0], qubit_indices[0])
	
	pauli_z_basis = circuit.copy("Z")
	pauli_z_basis.barrier()
	pauli_z_basis.measure(qubit_indices[0], qubit_indices[0])
	
	return [pauli_y_basis, pauli_x_basis, pauli_z_basis]
def get_state_tomography_clone2_for_post_selection(qubit_indices, circuit):
	#assert len(qubit_indices) == 2
	qasm_reprentation = circuit.qasm()
	qasm_reprentation = remove_measurements_from_qasm_string(qasm_reprentation)
	circuit = QuantumCircuit.from_qasm_str(qasm_reprentation)
	
	pauli_y_basis = circuit.copy("Y")
	pauli_y_basis.barrier()
	pauli_y_basis.sdg(qubit_indices[1])
	pauli_y_basis.h(qubit_indices[1])
	pauli_y_basis.measure(qubit_indices[1], qubit_indices[1])
	
	pauli_x_basis = circuit.copy("X")
	pauli_x_basis.barrier()
	pauli_x_basis.h(qubit_indices[1])
	pauli_x_basis.measure(qubit_indices[1], qubit_indices[1])
	
	pauli_z_basis = circuit.copy("Z")
	pauli_z_basis.barrier()
	pauli_z_basis.measure(qubit_indices[1], qubit_indices[1])
	
	return [pauli_y_basis, pauli_x_basis, pauli_z_basis]
def merge_dictionaries(d1, d2):
	combined = copy.deepcopy(d1)
	for k in d1:
		combined[k] = d1[k]+d2[k]
	return combined
