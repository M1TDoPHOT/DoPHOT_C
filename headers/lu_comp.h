/*
 * This code was (mostly) written by Ken Turkowski, who said:
 *
 * Oh, that. I wrote it in college the first time. It's open source - I think I
 * posted it after seeing so many people solve equations by inverting matrices
 * by computing minors naïvely.
 * -Ken
 *
 * The views represented here are mine and are not necessarily shared by
 * my employer.
   	Ken Turkowski			turk@apple.com
	Immersive Media Technologist 	http://www.worldserver.com/turk/
	Apple Computer, Inc.
	1 Infinite Loop, MS 302-3VR
	Cupertino, CA 95014
 */



/* This module solves linear equations in several variables (Ax = b) using
 * LU decomposition with partial pivoting and row equilibration.  Although
 * slightly more work than Gaussian elimination, it is faster for solving
 * several equations using the same coefficient matrix.  It is
 * particularly useful for matrix inversion, by sequentially solving the
 * equations with the columns of the unit matrix.
 *
 * lu_decompose() decomposes the coefficient matrix into the LU matrix,
 * and lu_solve() solves the series of matrix equations using the
 * previous LU decomposition.
 *
 *	Ken Turkowski (apple!turk)
 *	written 3/2/79, revised and enhanced 8/9/83.
 */

/* lu_decompose() decomposes the coefficient matrix A into upper and lower
 * triangular matrices, the composite being the LU matrix.
 *
 * The arguments are:
 *
 *	a - the (n x n) coefficient matrix
 *	n - the order of the matrix
 *
 *  1 is returned if the decomposition was successful,
 *  and 0 is returned if the coefficient matrix is singular.
 */

int lu_decompose(float **a, int n, int *ps, float **lu);

/* lu_solve() solves the linear equation (Ax = b) after the matrix A has
 * been decomposed with lu_decompose() into the lower and upper triangular
 * matrices L and U.
 *
 * The arguments are:
 *
 *	x - the solution vector
 *	b - the constant vector
 *	n - the order of the equation
*/

void lu_solve(float **lu, int n, int *ps, float *b, float *x);
