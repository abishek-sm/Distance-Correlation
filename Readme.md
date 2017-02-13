# Distance-Correlation

Added code to the Distance Correlation R package to handle multiple computations
# Description:

Distance Correlation is a relatively new method for computing relationships between non-linear vectors. I used this in my thesis research project to compute a correlation coefficient.( https://cran.r-project.org/web/packages/energy/energy.pdf )
Problem with the existing package: The existing R package can only compute the Distance Correlation (DC) value for 2 vectors. I needed the package to compute 2-by-2 DC values for a lot of vectors. i.e. an example computation had about 20,000 vectors where I had to compute C(20,000,2) computations. Doing it sequentially did not make sense because there was a lot of recomputation of intermediate distance matrices (https://arxiv.org/pdf/1010.0297.pdf)
# Optimizations

To reduce the number of computations, I cached intermediate results computed to be stored in memory and made extensive use of pointers to 2-D and 3-D arrays and deallocated them as needed. With my machine I used a 8 Gigs of memory for a significant speedup in the process.
# Result

I put up the entire code I had which now accepts a matrix of values and computes the distance correlation matrix and stores a square matrix of distance correlation values for each pair in the matrix. For example, if you gave it 20,000 columns with 500 values each, then the result is a 20,000 x 20,000 symmetric matrix with each element representing the value of the distance correlation. 
