# Ignis modification

`base_fitter.py` in Qiskit Ignis has been modified very slightly to accept strictly a dictionary of the form {"X": counts, "Y": counts, "Z": counts}, where counts is of the form {"0": int, "1": int}. Importantly the circuit names need to be X, Y and Z

The original `base_fitter.py` code from Qiskit Ignis can be found here: [https://github.com/Qiskit/qiskit-ignis/blob/master/qiskit/ignis/verification/tomography/fitters/base_fitter.py](https://github.com/Qiskit/qiskit-ignis/blob/master/qiskit/ignis/verification/tomography/fitters/base_fitter.py)

The relevant changes are on lines 73, 74, 75; those lines should be commented out. And 72 should be modified. 

```
# Add initial data
self._data = result
#if isinstance(result, Result):
#    result = [result]  # unify results handling
#self.add_data(result, circuits)
```

Lastly, line 46 need to be changed from `result: Union[Result, List[Result]],` to `result,`


