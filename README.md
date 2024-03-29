# Quantum-Telecloning

This repository contains python code and data associated with the papers [Quantum Telecloning on NISQ computers](https://arxiv.org/abs/2205.00125) and [Optimized Telecloning Circuits: Theory and Practice of Nine NISQ Clones](https://arxiv.org/abs/2210.10164)

The goal of this project is to implement several quantum telecloning algorithms on Quantinuum (formerly Honeywell) and IBMQ NISQ computers. 

For interacting with the IBMQ devices, one needs to set up an IBMQ account and load their credentials locally so that the Qiskit `IBMQ.load_account()` will load your account credentials. 

In order to execute the Quantinuum experiments one needs to first install the required Quantinuum API source code to interact with the backend; this should be in the form of a directory called `qtuum`

Then, for the purposes of computing the single state fidelity of the clones using Qiskit Ignis, a slight modification to the Qiskit Ignis source code (at least for `qiskit-ignis==0.7.1`) is required. This is explained in `modified_qiskit_ignis_code/`

The code that executes these quantum telecloning circuits on IBMQ and Quantinuum backends (`run_experiments/`) functions by saving job ids to local storage, which can then be used to retrieve the results after all circuits have been executed. 

However, these circuits are small enough to be easily simulated locally. The directory `classical_simulation_code` contains all of the relevant code to execute these local simulations. These local simulations are written entirely to be Qiskit compatible (with the slight modifications to Ignis source code), meaning they do not require any account setup from IBMQ or Quantinuum. In order to run the source code you will need to install python3 along with the required libraries in `requirements.txt`. You can install these requirements directly using the command `python3 -m pip install -r requirements.txt`

### Run experiments on IBMQ or Quantinuum devices
The python scripts that have the capability to send jobs to the various NISQ backends are in the directory `run_experiments`

### Fidelity figures
The directories `figures_Quantinuum`, `figures_IBMQ_post_select`, `figures_IBMQ_deferred_measurement` contain figures which show clone fidelities as a function of varying message states when the telecloning circuits are executed on NISQ devices. 

### Circuit drawings
The directory `circuit_drawings/` contains Qiskit circuit drawings for all telecloning circuit variants, including mid-circuit measurement with classical condition operations. 

Note that these circuit drawings do not include any state tomography operations on the clone qubits. 

## How to Cite?
bibtex:
```latex
@InProceedings{pelofske2022telecloning,
  author        = {Pelofske, Elijah and B{\"{a}}rtschi, Andreas and Garcia, Bryan and Kiefer, Boris and Eidenbenz, Stephan},
  booktitle     = {IEEE International Conference on Quantum Computing \& Engineering QCE'22},
  title         = {{Quantum Telecloning on NISQ Computers}},
  year          = {2022},
  doi           = {10.1109/QCE53715.2022.00083}
  archiveprefix = {arXiv},
  eprint        = {2205.00125},
}

@article{pelofske2022_M9_telecloning,
  doi = {10.48550/ARXIV.2210.10164},
  url = {https://arxiv.org/abs/2210.10164},
  author = {Pelofske, Elijah and Bärtschi, Andreas and Eidenbenz, Stephan},
  title = {Optimized Telecloning Circuits: Theory and Practice of Nine NISQ Clones},
  archiveprefix = {arXiv},
  year = {2022},
  eprint = {2210.10164},
}
```

Alternatively:
```
Elijah Pelofske, Andreas Bärtschi, Bryan Garcia, Boris Kiefer, and Stephan Eidenbenz. Quantum Telecloning on NISQ Computers. In IEEE International Conference on Quantum Computing & Engineering QCE’22, 2022. To appear. arXiv:2205.00125.
```

## Authors
- [Elijah Pelofske](mailto:epelofske@lanl.gov): Information Sciences, Los Alamos National Laboratory
- [Andreas Bärtschi](mailto:baertschi@lanl.gov): Information Sciences, Los Alamos National Laboratory
- Bryan Garcia: Department of Physics, New Mexico State University
- Boris Kiefer: Department of Physics, New Mexico State University
- Stephan Eidenbenz: Information Sciences, Los Alamos National Laboratory


## Copyright Notice:
© 2022. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.

**LANL C Number: C22038**

## License:
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
 
2.Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
 
3.Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
