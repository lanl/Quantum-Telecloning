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

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
import qiskit.quantum_info as qi
import numpy
import time
import math
from copy import deepcopy

''' message creation functions
    - message_ry does real amplitudes
    - message_ry adds complex amplitudes
'''
def message_ry(angle):
    qr_message = QuantumRegister(1, name='qM')
    qc_message = QuantumCircuit(qr_message)
    qc_message.ry(angle, 0)
    return qc_message
def message_ry_rz(angle_ry, angle_rz):
    qr_message = QuantumRegister(1, name='qM')
    qc_message = QuantumCircuit(qr_message)
    qc_message.ry(angle_ry, 0)
    qc_message.rz(angle_rz, 0)
    return qc_message

''' fracAngle(l,n) 
    is an angle theta, such that Ry(theta)|0⟩ = sqrt(l/n)|0⟩ + sqrt(1-l/n)|1⟩'''
def fracAngle(l,n):
    return 2*numpy.arccos(numpy.sqrt(l/n))


''' Split & Cyclic Shift Unitary SCS
    creates a circuit such that SCS(m) |0^(n-l)1^l⟩ = sqrt(l/n)|1 0^(n-l) 1^(l-1)⟩ + sqrt(1-l/n)|0^(n-l)1^l⟩
'''
def SCS(m=2):
    qr_SCS = QuantumRegister(m, name='qSCS{}{}'.format(m,m))
    qc_SCS = QuantumCircuit(qr_SCS)
    # circuit, little endian
    # blocki
    qc_SCS.ry(1.0*(numpy.pi)/2, 1)
    qc_SCS.cx(1,0)
    qc_SCS.ry(0.5*fracAngle(m-1,m),[0,1])
    qc_SCS.cx(1,0)
    qc_SCS.ry(-1.0*(numpy.pi)/2,1)
    # blockii
    for i in range(1,m-1):
        qc_SCS.cx(i,i+1)
        qc_SCS.ry(-0.25*fracAngle(i+1,m),i)
        qc_SCS.cx(i-1, i)
        qc_SCS.ry(0.25*fracAngle(i+1,m),i)
        qc_SCS.cx(i+1, i)
        qc_SCS.ry(-0.25*fracAngle(i+1,m),i)
        qc_SCS.cx(i-1, i)
        qc_SCS.ry(0.25*fracAngle(i+1,m),i)
        qc_SCS.cx(i,i+1)
    # return
    return qc_SCS


''' Dicke State Unitary DSU
    creates a circuit such that DSU(m) |0^(n-l)1^l⟩ = |D^n_l⟩
'''
def DSU(m=2):
    qr_DSU = QuantumRegister(m, name='qU{}{}'.format(m,m))
    qc_DSU = QuantumCircuit(qr_DSU)
    # circuit, little endian
    for l in range(m,1,-1):
        qc_DSU.compose(SCS(l), qubits=[*range(0,l)], inplace=True)
    # return
    return qc_DSU
''' old version '''
def DSU_old(m=2):
    qr_DSU = QuantumRegister(m, name='qU{}{}'.format(m,m))
    qc_DSU = QuantumCircuit(qr_DSU)
    # circuit, little endian
    for l in range(m,1,-1):
        # blocki
        qc_DSU.ry(1.0*(numpy.pi)/2, 1)
        qc_DSU.cx(1,0)
        qc_DSU.ry(0.5*fracAngle(l-1,l),[0,1])
        qc_DSU.cx(1,0)
        qc_DSU.ry(-1.0*(numpy.pi)/2,1)
        # blockii
        for i in range(1,l-1):
            qc_DSU.cx(i,i+1)
            qc_DSU.ry(-0.25*fracAngle(i+1,l),i)
            qc_DSU.cx(i-1, i)
            qc_DSU.ry(0.25*fracAngle(i+1,l),i)
            qc_DSU.cx(i+1, i)
            qc_DSU.ry(-0.25*fracAngle(i+1,l),i)
            qc_DSU.cx(i-1, i)
            qc_DSU.ry(0.25*fracAngle(i+1,l),i)
            qc_DSU.cx(i,i+1)
    # return
    return qc_DSU


''' inputHW
    - creates the entangled Hamming Weight Superpositions to be fed into Dicke State Unitaries
    - general method for Ancilla TeleCloning State, for both LNN and Full Connectivity
    - methods for m=2 and m=3 without ancillas, for both LNN and Full connectivity
'''
def inputHW(m=2,ancilla=True,topology='LNN'):
    qr_HW = QuantumRegister(2*m if ancilla else m+1, name='qHW')
    qc_HW = QuantumCircuit(qr_HW)
    # TCstate A^(m-1)PC^m
    if (ancilla):
        # entangling Hamming weights on LNN connectivity
        if (topology == 'LNN'):
            qc_HW.ry(fracAngle(1,m+1),m-1)
            qc_HW.cx(m-1,m)
            for i in range(1,m):
                # sifts
                qc_HW.cx([*range(m,2*m-i)], [*range(m+1,2*m-i+1)])
                qc_HW.cx([*range(m+1,2*m-i+1)], [*range(m,2*m-i)])
                # increase HW
                qc_HW.ry(0.5*fracAngle(m-i,m-i+1), m-2)
                qc_HW.cx(m-2,m-1)
                qc_HW.cx(m-1,m-2)
                qc_HW.ry(-0.5*fracAngle(m-i,m-i+1), m-1)
                qc_HW.cx(m-1,m)
                # sifts
                if (i<=m-2):
                    qc_HW.cx([*range(m-2,i-1,-1)], [*range(m-3,i-2,-1)])
                    qc_HW.cx([*range(m-3,i-2,-1)],[*range(m-2,i-1,-1)])
        # entangling Hamming weights on full connectivity
        else: # topology == 'Full'
            qc_HW.ry(fracAngle(1,m+1),0)
            for i in range(1,m):
                qc_HW.ry(0.5*fracAngle(m-i,m-i+1), i)
                qc_HW.cx(i-1,i)
                qc_HW.ry(-0.5*fracAngle(m-i,m-i+1), i)
            qc_HW.cx([*range(0,m)],[*range(2*m-1,m-1,-1)]) 
    # TCstate PC^m
    else:
        assert m<=3, "no telecloning circuit without ancilla known for m>3"
        if (m == 2):
            qc_HW.ry(fracAngle(2,3),0)
            qc_HW.cx(0,1)
        elif (m == 3):
            if (topology == 'LNN'):
                qc_HW.ry(fracAngle(2,3),0)
                qc_HW.cx(0,1)
                qc_HW.ry(0.5*fracAngle(1,2),2)         
                qc_HW.cx(1,2)
                qc_HW.ry(-0.5*fracAngle(1,2),2)
                qc_HW.x(0)
                qc_HW.ry(0.5*fracAngle(1,4),1)         
                qc_HW.cx(0,1)
                qc_HW.ry(-0.5*fracAngle(1,4),1)
                qc_HW.x(0)                
            else:
                qc_HW.ry(fracAngle(2,3),0)
                qc_HW.x(1)
                qc_HW.ry(0.5*fracAngle(1,2),2)         
                qc_HW.cx(0,2)
                qc_HW.ry(-0.5*fracAngle(1,2),2)
                qc_HW.x(0)
                qc_HW.ry(-0.5*fracAngle(3,4),1)         
                qc_HW.cx(0,1)
                qc_HW.ry(0.5*fracAngle(3,4),1)
                qc_HW.x(0)   
                
    # return
    return qc_HW

''' TCstate
    creates TeleCloning State, depending on availability of ancillas:
    - A^(m-1) P C^m: m-1 Ancillas, 1 Port qubits, m Clone qubits
    - P C^m: 1 Port qubit, m Clone qubits, only for m<=3    
    - anc_opt:  Optimizes DSU(m) to SCS(m) on the ancilla+port qubits
                False   for legacy reasons / backwards compatibility to QCE paper
                True    for newer Quantinuum experiments 
'''
def TCstate(m=2,ancilla=True,topology='LNN',anc_opt=False):
    qr_TCstate = QuantumRegister(2*m if ancilla else m+1, name='qTCstate')
    qc_TCstate = QuantumCircuit(qr_TCstate)

    # set input Hamming weight
    qc_TCstate.compose(inputHW(m,ancilla,topology), inplace=True)
    # feed into dicke state unitary(ies)
    if (ancilla):
        if (anc_opt):
            qc_TCstate.compose(SCS(m), qubits=[*range(0,m)], inplace=True)
        else:
            qc_TCstate.compose(DSU(m), qubits=[*range(0,m)], inplace=True)
        qc_TCstate.compose(DSU(m), qubits=[*range(2*m-1,m-1,-1)], inplace=True)
    else:
        qc_TCstate.compose(DSU(m), qubits=[*range(1,m+1)], inplace=True)
    
    # return
    return qc_TCstate


''' LOCC
    different implementations of Local Operations and Classical Communication:
    1) Bell measurement of Message and Port Qubit
    2) X, Z operations on Clone Qubits depending on measurement
    3) Returns a dictionary of Qubit Positions
    
    locc can be done through:
    - 'dfm':  deferred measurement / quantum operations
    - 'psXZ': postselection to apply X,Z 
    - 'locc': classical feed forward
'''
def LOCC(m=2,ancilla=True,topology='LNN',locc='dfm', Qiskit_if_statements=True):
    HQS_QASM = ""
    # circuit
    if (locc == 'locc'):
        qr_LOCC = QuantumRegister(2*m+1 if ancilla else m+2, name='qLOCC')    
        cr_LOCC = ClassicalRegister(2*m+1 if ancilla else m+2, name='cLOCC')
        qc_LOCC = QuantumCircuit(qr_LOCC, cr_LOCC)
    else:
        qr_LOCC = QuantumRegister(2*m+1 if ancilla else m+2, name='qLOCC')    
        qc_LOCC = QuantumCircuit(qr_LOCC)
        
    # qubit positions
    if (ancilla):
        qpos = {'message': 0, 'ancillas': [*range(1,m)], 'port': m, 'clones': [*range(m+1,2*m+1)]}
    else:
        qpos = {'message': 0, 'ancillas': [], 'port': 1, 'clones': [*range(2,m+2)]}
    
    # implement locc
    if (locc == 'locc'):
        qc_LOCC.cx(qpos['message'], qpos['port']) 
        qc_LOCC.h(qpos['message'])
        qc_LOCC.measure(qpos['message'], qpos['message'])
        qc_LOCC.measure(qpos['port'], qpos['port'])
        for clone in qpos['clones']:
            if Qiskit_if_statements:
                qc_LOCC.x(clone).c_if(cr_LOCC[qpos['port']], 1)
                qc_LOCC.z(clone).c_if(cr_LOCC[qpos['message']], 1)         
            else:
                HQS_QASM += "if(c["+str(qpos['port'])+"]==1) x q["+str(clone)+"];\n"
                HQS_QASM += "if(c["+str(qpos['message'])+"]==1) z q["+str(clone)+"];\n"
    # implement deferred measurement
    elif (locc == 'dfm'):
        if (topology == 'LNN'):
            # bell measurement cnot
            qc_LOCC.cx(qpos['message'], qpos['port']) 
            # port-controlled clone-Xgates, combined with swaps
            qc_LOCC.cx(qpos['clones'], [qpos['port'],*qpos['clones'][0:m-1]])
            qc_LOCC.cx([qpos['port'],*qpos['clones'][0:m-1]], qpos['clones'])
            qpos['clones'], qpos['port'] = [qpos['port'],*qpos['clones'][0:m-1]], qpos['clones'][m-1]
            # message-controlled clone-Zgates, combined with swaps
            qc_LOCC.cx([qpos['message'],*qpos['clones'][0:m-2]], qpos['clones'][0:m-1])
            qc_LOCC.cx(qpos['clones'][0:m], [qpos['message'],*qpos['clones'][0:m-1]]) 
            qpos['clones'][0:m-1], qpos['message'] = [qpos['message'],*qpos['clones'][0:m-2]], qpos['clones'][m-2]
            #qc_LOCC.h(qpos['message']) # bell measurement h, not necessary for dfm
        else:
            qc_LOCC.cx(qpos['message'], qpos['port'])       # bell measurement cnot
            qc_LOCC.cx([qpos['port']]*m, qpos['clones'])    # port-controlled clone-Xgates
            qc_LOCC.cx(qpos['clones'], [qpos['message']]*m) # message-controlled clone-Zgates
            #qc_LOCC.h(qpos['message']) # bell measurement h, not necessary for dfm


    # implement postselection
    elif (locc[0:2] == 'ps'):
        qc_LOCC.cx(qpos['message'], qpos['port']) # bell measurement cnot
        qc_LOCC.h(qpos['message'])                      # bell measurement h
        if (locc[2] == '1'): qc_LOCC.x(qpos['clones'])
        if (locc[3] == '1'): qc_LOCC.z(qpos['clones'])
        #print(locc[2], locc[3], locc)

    
    # return
    return qc_LOCC, qpos, HQS_QASM

''' construct_circuit
    Main Function for all circuits, with options:
    - ryangle: message angle y-rotation, in [0,pi]
    - rzangle: message angle z-rotation, in [0,2*pi]
    - m: number of clones > 1
    - ancilla:  True/False depending on whether ancilla are allowed for the TCstate
                must be True if n>3
    - topology: 'LNN'  if linear nearest neighbor between TCstate + edge between message & port
                'full' otherwise (only used/defined as an ELSE statement to IF 'LNN')
    - locc:     'locc' for local operations and classical (feed-forward) communication (Honeywell only)
                'dfm'  for deferred measurement / quantum operations
                'ps00' for postselecting 00, no X-gate, no Z-gate
                'ps01' for postselecting 01, Z-gate
                'ps10' for postselecting 10, X-gate
                'ps11' for postselecting 11, X-gate followed by Z-gate
    - anc_opt:  Optimizes DSU(m) to SCS(m) on the ancilla+port qubits
                False   for legacy reasons / backwards compatibility to QCE paper
                True    for newer Quantinuum experiments 
'''
def construct_circuit(ryangle, rzangle, m, ancilla, topology, locc, Qiskit_version_if_statements, anc_opt=False):
    # setup circuit
    width = 2*m+1 if ancilla else m+2
    qc = QuantumCircuit(width, width)
    # setup message, TCstate, locc
    qc_message = message_ry_rz(ryangle, rzangle)
    rho = qi.DensityMatrix.from_instruction(qc_message)
    qc_TCstate = TCstate(m,ancilla,topology,anc_opt)
    qc_LOCC, qpos, HQS_QASM = LOCC(m,ancilla,topology,locc, Qiskit_version_if_statements)
    # compose circuit
    qc.compose(qc_message, qubits=[0], inplace=True)
    qc.compose(qc_TCstate, qubits=[*range(1,width)], inplace=True)
    qc.barrier()
    qc.compose(qc_LOCC, qubits=[*range(0,width)], inplace=True)
    
    # return circuit, clone positions, message density matrix, postselect positions & values (if any)
    if (locc[0:2] == 'ps'):
        return qc, qpos['clones'], rho, {qpos['port']: int(locc[2]), qpos['message']: int(locc[3])}
    if Qiskit_version_if_statements:
        return qc, qpos['clones'], rho
    else:
        return qc, qpos['clones'], rho, HQS_QASM    

''' IBM
    - 2 clones
    - no ancillas
'''
def construct_deferred_measurement_PCC(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=False, topology='LNN', locc='dfm', Qiskit_version_if_statements=True)

def construct_PCC_postselect_00(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=False, topology='LNN', locc='ps00', Qiskit_version_if_statements=True)

def construct_PCC_postselect_01(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=False, topology='LNN', locc='ps01', Qiskit_version_if_statements=True)

def construct_PCC_postselect_10(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=False, topology='LNN', locc='ps10', Qiskit_version_if_statements=True)

def construct_PCC_postselect_11(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=False, topology='LNN', locc='ps11', Qiskit_version_if_statements=True)

''' IBM
    - 2 clones
    - with ancillas
'''
def construct_deferred_measurement_APCC(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=True, topology='LNN', locc='dfm', Qiskit_version_if_statements=True)

def construct_APCC_postselect_00(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=True, topology='LNN', locc='ps00', Qiskit_version_if_statements=True)

def construct_APCC_postselect_01(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=True, topology='LNN', locc='ps01', Qiskit_version_if_statements=True)

def construct_APCC_postselect_10(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=True, topology='LNN', locc='ps10', Qiskit_version_if_statements=True)

def construct_APCC_postselect_11(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=True, topology='LNN', locc='ps11', Qiskit_version_if_statements=True)

''' IBM
    - 3 clones
    - no ancillas
'''
def construct_deferred_measurement_PCCC(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=False, topology='LNN', locc='dfm', Qiskit_version_if_statements=True)

def construct_PCCC_postselect_00(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=False, topology='LNN', locc='ps00', Qiskit_version_if_statements=True)

def construct_PCCC_postselect_01(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=False, topology='LNN', locc='ps01', Qiskit_version_if_statements=True)

def construct_PCCC_postselect_10(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=False, topology='LNN', locc='ps10', Qiskit_version_if_statements=True)

def construct_PCCC_postselect_11(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=False, topology='LNN', locc='ps11', Qiskit_version_if_statements=True)

''' IBM
    - 3 clones
    - with ancillas
'''
def construct_deferred_measurement_AAPCCC(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=True, topology='LNN', locc='dfm', Qiskit_version_if_statements=True)

def construct_AAPCCC_postselect_00(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=True, topology='LNN', locc='ps00', Qiskit_version_if_statements=True)

def construct_AAPCCC_postselect_01(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=True, topology='LNN', locc='ps01', Qiskit_version_if_statements=True)

def construct_AAPCCC_postselect_10(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=True, topology='LNN', locc='ps10', Qiskit_version_if_statements=True)

def construct_AAPCCC_postselect_11(ryangle, rzangle):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=True, topology='LNN', locc='ps11', Qiskit_version_if_statements=True)


''' Honeywell
    - 2 clones
'''
def construct_PCC_Honeywell_QASM_circuit_Qiskit_version(ryangle, rzangle, Qiskit_version_if_statements):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=False, topology='full', locc='locc', Qiskit_version_if_statements=Qiskit_version_if_statements)

def construct_APCC_Honeywell_QASM_circuit_Qiskit_version(ryangle, rzangle, Qiskit_version_if_statements):
    return construct_circuit(ryangle, rzangle, m=2, ancilla=True, topology='full', locc='locc', Qiskit_version_if_statements=Qiskit_version_if_statements)

''' Honeywell
    - 3 clones
'''
def construct_PCCC_Honeywell_QASM_circuit_Qiskit_version(ryangle, rzangle, Qiskit_version_if_statements):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=False, topology='full', locc='locc', Qiskit_version_if_statements=Qiskit_version_if_statements)

def construct_AAPCCC_Honeywell_QASM_circuit_Qiskit_version(ryangle, rzangle, Qiskit_version_if_statements):
    return construct_circuit(ryangle, rzangle, m=3, ancilla=True, topology='full', locc='locc', Qiskit_version_if_statements=Qiskit_version_if_statements)

