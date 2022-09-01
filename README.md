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
lanczos_acc:&emsp;&ensp;OpenACC accelerated Lanczos algorithm for gauge-covariant lattice Laplace operator in three dimensions (operator and<br />
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; gram-schmidt orthogonalization).<br />
