Data Format:

```

vlen

site : direction
x1 x2 x3
x4 x5 x6
x7 x8 x9

.
.
.
.
.

<vlen * 3 times>

```

Where,&emsp; vlen		- 3D grid dimension times 3 <br />
&emsp;&emsp;&emsp;&emsp; direction	- 0/1/2 for right/bottom/hind <br />
&emsp;&emsp;&emsp;&emsp; x1 x2 x3 &emsp; --¬ <br />
&emsp;&emsp;&emsp;&emsp; x4 x5 x6 &emsp;&emsp; |-- 3 x 3 complex numbers <br />
&emsp;&emsp;&emsp;&emsp; x7 x8 x9 &emsp; __|<br />

Sample files' naming conventions: <br />

sample_edgedata_\<vlen\>_\<l\>_\<u\>txt	- gluon matrices are complex numbers distributed in U(l,u).<br />
sample_su3_\<vlen\>.txt			- gluon matrices are SU3 simulated data.<br />

NOTE: <br />
Find the larger files at: https://drive.google.com/drive/folders/1wI8U1A1Rwb2TMxi2UUTEJMoe6G2iubAY?usp=sharing
