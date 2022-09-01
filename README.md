# acc-lanczos
OpenACC accelerated lambda-lanczos

This is an OpenACC accelerated version of [lambda-lanczos](https://github.com/mrcdr/lambda-lanczos)

Contains 4 versions after compiling the whole project:

```
make full
```

<br />
lanczos_std:&emsp;&ensp;Sample program for testing the Grid logic for calculating Laplace operator (2D and 3D available).<br />
<br />
lanczos_cpu:&emsp; Lanczos algorithm for gauge-covariant lattice Laplace operator in three dimensions.<br />
<br />
lanczos_gpu:&emsp; OpenACC accelerated Lanczos algorithm for gauge-covariant lattice Laplace operator in three dimensions (operator only).<br />
<br />
lanczos_acc:&emsp;&ensp;OpenACC accelerated Lanczos algorithm for gauge-covariant lattice Laplace operator in three dimensions (operator and gram-schmidt orthogonalization).<br />

---

Prerequisites:<br />
<br />
* g++ (with support for C++11)<br />
* nvc++ (as part of NVIDIA HPC SDK)<br />

Requires syetem with Nvidia CUDA-enabled GPU. If multiple GPUs are available, set ```CUDA_VISIBLE_DEVICES``` to one of them.

---
Benchmarking results for SU3 input data for 3D grid sizes: 10<sup>3</sup>, 20<sup>3</sup>, 30<sup>3</sup>, 40<sup>3</sup> and 50<sup>3</sup>

![Time taken by serial, partial-parallelized and fullly parallel versions](https://github.com/Beck-919/acc-lanczos/blob/master/stats/line_time.png?raw=true)

![Speedup of partial-parallelized and fullly parallel versions with serial baseline](https://github.com/Beck-919/acc-lanczos/blob/master/stats/bar_speedup.png?raw=true)
