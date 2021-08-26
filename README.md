# Overview
This directory has source codes which realize neural network simulation based on Hodgkin-Huxley model.

# Requirements
* Clang
* Python

# How to use
## Parameters
If you need to change parameters which are related to Hodgkin-Huxley model, you can change them in the file, `src/*.h`.

## Run simulation
You can simulate neural network acrivity by the following command:
```
$ cd src
$ ./run.sh
```

## Plot results
If you want figures of neural network activity sucn as membrane potential, you can get them by the following command:
```
$ python plot.py
```

# References
* A. L. Hodgkin and A. F. Huxley, [A quantitative description of membrance current and its application to conduction and excitation in nerve](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/), J Physiol, 1952.
