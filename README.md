# Code to accompany ["Escaping the Shadow of Bell's Theorem in Network Nonlocality"](https://arxiv.org/abs/2406.15587)

## Maria Ciudad Alañón,  Emanuel-Cristian Boghiu, Paolo Abiuso and Elie Wolfe

This repository contains the codes used to obtain all the results in "Escaping the Shadow of Bell's Theorem in Network Nonlocality". Maria Ciudad Alañón, Emanuel-Cristian Boghiu, Paolo Abiuso and Elie Wolfe.

All the code is written in Python. They use the solver Gurobi to solve linear and bilinear problems.

The files used for the results of the paper are:

- [test_one_classical_source_AB.py](test_one_classical_source_AB.py): this file contains a function that given a probability distribution in the bilocality scenario (with the minimal configuration), it tells you whether is compatible with one classical source between A and B and one nonclassical source between B and C.
- [test_one_classical_source_BC.py](test_one_classical_source_BC.py): this file contains a function that given a probability distribution in the bilocality scenario (with the minimal configuration), it tells you whether is compatible with one nonclassical source between A and B and one classical source between B and C.
- [test_bilocality.py](test_bilocality.py): this file contains a function that given a probability distribution in the bilocality scenario, it tells you whether is compatible with the two sources being classical.
- [test_bilocal_minimal_configuration.py](test_bilocal_minimal_configuration.py): this file is an alternative formulation of test_bilocality.py but with the restriction that it can only be used for the minimal configuration of the bilocality scenario.

As explained in the paper, with minimal configuration we mean binary outputs for the three parties and binary inputs for the extreme parties (while no input for the middle-party).
