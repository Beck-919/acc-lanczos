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

Where,	vlen		- 3D grid dimension times 3 <br />
	direction	- 0/1/2 for right/bottom/hind <br />
	x1 x2 x3	--¬ <br />
	x4 x5 x6	  |-- 3 x 3 complex numbers <br />
	x7 x8 x9	__|<br />


Sample files filenames:

sample_edgedata_<vlen>_<l>_<u>txt	- gluon matrices are complex numbers distributed in U(l,u).
sample_su3_<vlen>.txt			- gluon matrices are SU3 simulated data.

NOTE: Find the larger files at: https://drive.google.com/drive/folders/1wI8U1A1Rwb2TMxi2UUTEJMoe6G2iubAY?usp=sharing
