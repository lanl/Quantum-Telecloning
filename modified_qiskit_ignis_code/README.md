# Qiskit Ignis modification

Qiskit Ignis is now deprecated. However, for reconstructing the clone fidelities we use a modified version of Qiskit Ignis. Specifically, a modified version of `qiskit-ignis==0.7.1`. 

The relevant source code change lets us submit tomography data of the form {"X": counts, "Y": counts, "Z": counts}, where counts is of the form {"0": int, "1": int}. Importantly the circuit names need to be X, Y and Z. 

The code that needs to be changed is `base_fitter.py`; [https://github.com/Qiskit/qiskit-ignis/blob/master/qiskit/ignis/verification/tomography/fitters/base_fitter.py](https://github.com/Qiskit/qiskit-ignis/blob/master/qiskit/ignis/verification/tomography/fitters/base_fitter.py)

You will need to install qiskit ignis (`qiskit-ignis==0.7.1`) and then make the following local changes:

- Lines 72, 73, 74, 75. Originally they should look like this:

```
72:        self._data = {}
73:        if isinstance(result, Result):
74:            result = [result]  # unify results handling
75:        self.add_data(result, circuits)
```

Modify these lines to be:

```
        self._data = result
        #if isinstance(result, Result):
        #    result = [result]  # unify results handling
        #self.add_data(result, circuits)
```
