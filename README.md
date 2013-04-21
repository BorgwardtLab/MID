MID (Mutual Information Dimension)
==================================

MID for measuring statistical dependence between two random variables.


Summary
-------

An estimation algorithm for MID (Mutual Information Dimension), which measures statistical dependence between two random variables.
This algorithm has the following advantages:

* **Nonlinear dependences** (and also linear dependences) can be measured,
* **Scalable**; the average-case time complexity is *O*(*n*log*n*), where *n* is the number of data points, and
* **Parameter-free**.

Please see the following article for detailed information and refer it in your published research:

* Sugiyama, M., Borgwardt, K. M.: **Measuring Statistical Dependence via the Mutual Information Dimension**,
	*Proceedings of the 23rd International Joint Conference on Artificial Intelligence (IJCAI 2013)*, to appear.


Installation
------------

The code consists of only one C file "MID.c".
Thus you can use by compiling it, for example, type into your terminal:

	$ gcc -O3 MID.c -o MID


Usage
-----

To calculate MID between two variables, type:

	$ ./MID <input_file>
	
`<input_file>` is a comma-separated text file with two columns without row and column names.
Columns correspond to respective variables.
For example,	

	0.921,0.930
	0.491,0.492
	0.990,0.993
	0.775,0.777
	...
	0.577,0.561

The followings are shown at standard output.

* **dimX** (the information dimension of the first variable)
* **dimY** (the information dimension of the second variable)
* **dimXY** (the information dimension of X and Y)
* **MID** (equivalent to dimX + dimY - dimXY)

Example
-------

	$ gcc -O3 MID.c -o MID
	$ ./MID ./sampledata/linear.csv
	idim_x:  0.994690
	idim_y:  0.994690
	idim_xy: 0.994690
	MID:     0.994690
	$ ./MID ./sampledata/noise.csv
	idim_x:  0.995130
	idim_y:  0.996233
	idim_xy: 1.755107
	MID:     0.236256


Information
-----------

* Author: Mahito Sugiyama
* Affiliation: Machine Learning & Computational Biology Research Group, MPIs Tübingen, Germany
* Mail: mahito.sugiyama@tuebingen.mpg.de