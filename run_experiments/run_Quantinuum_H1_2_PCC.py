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
import numpy as np
import math
from utils import *
from qtuum.api_wrappers import QuantinuumAPI as QAPI
from create_telecloning_circuits import *

machine = 'H1-2'
qapi = QAPI(machine=machine)
status = qapi.status()
print(machine, ":", status)

def run(angle_ry, angle_rz):
	qc, clone_indices, rho, HQS_QASM = construct_PCC_Honeywell_QASM_circuit_Qiskit_version(angle_ry, angle_rz, False)
	QASM = qc.qasm()+HQS_QASM
	QASM = QASM.replace('qelib1.inc', 'hqslib1.inc')
	tomography_circuit_Y, tomography_circuit_X, tomography_circuit_Z = append_parallel_qubit_tomography_to_Honeywell_QASM(clone_indices, QASM)
	
	job_id = qapi.submit_job(tomography_circuit_Y, shots=300, machine='H1-2')
	print(job_id)
	file = open("Quantinuum_H1_2_results/PCC_job_"+str(angle_ry)+"_"+str(angle_rz)+"_Y.txt", "w")
	file.write(job_id)
	file.close()
	
	job_id = qapi.submit_job(tomography_circuit_X, shots=300, machine='H1-2')
	print(job_id)
	file = open("Quantinuum_H1_2_results/PCC_job_"+str(angle_ry)+"_"+str(angle_rz)+"_X.txt", "w")
	file.write(job_id)
	file.close()
	
	job_id = qapi.submit_job(tomography_circuit_Z, shots=300, machine='H1-2')
	print(job_id)
	file = open("Quantinuum_H1_2_results/PCC_job_"+str(angle_ry)+"_"+str(angle_rz)+"_Z.txt", "w")
	file.write(job_id)
	file.close()

N_Angle_slices = 6

angles_ry = np.linspace(0, math.pi, num=N_Angle_slices)
print(angles_ry)

for angle_ry in angles_ry:
	print(angle_ry)
	run(angle_ry, math.pi/2.0)
