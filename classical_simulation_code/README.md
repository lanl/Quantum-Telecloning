# Local Simulation

Assuming that you have your local Qiskit Ignis modifications in place, then you will be able to execute any of these python3 scripts in order to simulate the telecloning circuits. 

These scripts simply execute noiseless simulations of the telecloning circuits with state tomography with different message qubit states; the computed fidelity of the resulting clones is then printed. 

As an example, if you execute `python3 run_postselect_local_PCC.py`, the output will look like this:

```
Ry angle = 0.0
Rz angle = 0.0
Clone 1 fidelity: 0.8317095068666248
Clone 2 fidelity: 0.8341369334619103
***
Ry angle = 0.0
Rz angle = 0.06346651825433926
Clone 1 fidelity: 0.8338608660515862
Clone 2 fidelity: 0.8338608660515862
***
Ry angle = 0.0
Rz angle = 0.12693303650867852
Clone 1 fidelity: 0.8335387587861028
Clone 2 fidelity: 0.8330723874879237
***
Ry angle = 0.0
Rz angle = 0.1903995547630178
Clone 1 fidelity: 0.8308090872857722
Clone 2 fidelity: 0.8357247243257603
***
Ry angle = 0.0
Rz angle = 0.25386607301735703
Clone 1 fidelity: 0.8312324929971991
Clone 2 fidelity: 0.8349673202614382
***
Ry angle = 0.0
Rz angle = 0.3173325912716963
Clone 1 fidelity: 0.8391994707244451
Clone 2 fidelity: 0.8351306649024136
***
Ry angle = 0.0
Rz angle = 0.3807991095260356
Clone 1 fidelity: 0.8358945108039191
Clone 2 fidelity: 0.8321030734129652
```

