#include "gaussjordan_3d.h"
/******************************************************************************/

/* Perform Gauss-Jordan elimination with row-pivoting to obtain the solution to

* the system of linear equations

* A X = B

*

* Arguments:

*                 lhs                -        left-hand side of the equation, matrix A

*                 rhs                -        right-hand side of the equation, matrix B

*                 nrows        -        number of rows in the arrays lhs and rhs

*                 ncolsrhs-        number of columns in the array rhs

*

* The function uses Gauss-Jordan elimination with pivoting.  The solution X to

* the linear system winds up stored in the array rhs; create a copy to pass to

* the function if you wish to retain the original RHS array.

*

* Passing the identity matrix as the rhs argument results in the inverse of

* matrix A, if it exists.

*

* No library or header dependencies, but requires the function swaprows, which

* is included here.

*/

namespace VecMath
{

	//  swaprows - exchanges the contents of row0 and row1 in a 2d array

	void swaprows3d(double** arr, long row0, long row1) {

		double* temp;

		temp = arr[row0];

		arr[row0] = arr[row1];

		arr[row1] = temp;

	}



	//        gjelim

	void gjelim3d(S3DMatrix* lhs, S3DMatrix* rhs, long nrows, long ncolsrhs) {



		//        augment lhs array with rhs array and store in arr2

		double** arr2 = new double* [nrows];

		for (long row = 0; row < nrows; ++row)

			arr2[row] = new double[nrows + ncolsrhs];



		for (long row = 0; row < nrows; ++row) {

			for (long col = 0; col < nrows; ++col) {

				arr2[row][col] = lhs->a[row][col];

			}

			for (long col = nrows; col < nrows + ncolsrhs; ++col) {

				arr2[row][col] = rhs->a[row][col - nrows];

			}

		}



		//        perform forward elimination to get arr2 in row-echelon form

		for (long dindex = 0; dindex < nrows; ++dindex) {

			//        run along diagonal, swapping rows to move zeros in working position

			//        (along the diagonal) downwards

			if ((dindex == (nrows - 1)) && (arr2[dindex][dindex] == 0)) {

				return; //  no solution

			}
			else if (arr2[dindex][dindex] == 0) {

				swaprows3d(arr2, dindex, dindex + 1);

			}

			//        divide working row by value of working position to get a 1 on the

			//        diagonal

			if (arr2[dindex][dindex] == 0.0) {

				return;

			}
			else {

				double tempval = arr2[dindex][dindex];

				for (long col = 0; col < nrows + ncolsrhs; ++col) {

					arr2[dindex][col] /= tempval;

				}

			}



			//        eliminate value below working position by subtracting a multiple of

			//        the current row

			for (long row = dindex + 1; row < nrows; ++row) {

				double wval = arr2[row][dindex];

				for (long col = 0; col < nrows + ncolsrhs; ++col) {

					arr2[row][col] -= wval * arr2[dindex][col];

				}

			}

		}



		//        backward substitution steps

		for (long dindex = nrows - 1; dindex >= 0; --dindex) {

			//        eliminate value above working position by subtracting a multiple of

			//        the current row

			for (long row = dindex - 1; row >= 0; --row) {

				double wval = arr2[row][dindex];

				for (long col = 0; col < nrows + ncolsrhs; ++col) {

					arr2[row][col] -= wval * arr2[dindex][col];

				}

			}

		}



		//        assign result to replace rhs

		for (long row = 0; row < nrows; ++row) {

			for (long col = 0; col < ncolsrhs; ++col) {

				rhs->a[row][col] = arr2[row][col + nrows];

			}

		}



		for (long row = 0; row < nrows; ++row)

			delete[] arr2[row];

		delete[] arr2;

	}
}