#include "fonction_qr_mpfr.h"

// Memory allocation function of n by n matrix
mpfr_t* allocateMatrix(int n) {
	return malloc(n * n * sizeof(mpfr_t));
}

// Memory allocation function of a vector (length n)
mpfr_t* allocateVector(int n) {
	return malloc(n * sizeof(mpfr_t));
}

//// The MPFR library is to be used with initialised mpfr_t variable 
//// which will contain multiprecision floats.
// Initialise a mpfr_t square Matrix (dimension n)
void initMatrix(mpfr_t* A, int n, mpfr_prec_t prec) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_init2(A[i*n + j], prec);
}

// Initialise a mpfr_t Vector (dimension n)
void initVector(mpfr_t* v, int n, mpfr_prec_t prec) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_init2(v[i], prec);
}
////

// Free the memory of a mpfr_t square Matrix (dimension n)
void clearMatrix(mpfr_t* A, int n) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_clear(A[i*n + j]);
	free(A);
}

// Free the memory of a mpfr_t Vector (dimension n)
void clearVector(mpfr_t* v, int n) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_clear(v[i]);
	free(v);
}

// Compute the 2-norm of a Vector of length n
// The output is stored in renv
// the Vector x is stored like this :
// x_0 = x[incx * 0]
// x_1 = x[incx * 1]
// ...
// x_k = x[incx * k]
void norm2(mpfr_t renv, int n, mpfr_t* x, int incx) {
	mpfr_set_d(renv, 0, MPFR_RNDN);
	int i, tmp_indx;
	for(i = 0, tmp_indx = 0 ; i < n ; i++, tmp_indx += incx)
		mpfr_hypot(renv, x[tmp_indx], renv, MPFR_RNDN);
}

// Affect to all the element of the Matrix A the value 0
void set0Matrix(mpfr_t* A, int n) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_set_d(A[i*n + j], 0, MPFR_RNDN);
}

// Affect to all the element of the Vector v the value 0
void set0Vector(mpfr_t* v, int n) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_set_d(v[i], 0, MPFR_RNDN);
}

// Affect randomly the element of the Matrix A such as every line form unitary vector
void setAleaMatrix(mpfr_t* A, int n, mpfr_t tmp) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_set_d(A[i*n + j], ((double) rand()) / RAND_MAX, MPFR_RNDN);
	for(i = 0 ; i < n ; i++) {
		norm2(tmp, n, A + i, 1);
		for(j = 0 ; j < n ; j++)
			mpfr_div(A[i*n + j], A[i*n + j], tmp, MPFR_RNDN); // A[i*n + j] /= tmp
	}
}

// Affect randomly the element of the vectot v such as v is an unitary vector
void setAleaVector(mpfr_t* v, int n, mpfr_t tmp) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_set_d(v[i], ((double) rand()) / RAND_MAX, MPFR_RNDN);
	norm2(tmp, n, v, 1);
	for(i = 0 ; i < n ; i++)
		mpfr_div(v[i], v[i], tmp, MPFR_RNDN); // v[i] /= tmp;
}

// A := In (unit Matrix of size nxn)
void identityMatrix(mpfr_t* A, int n) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_set_d(A[i*n + j], 0, MPFR_RNDN);
	for(i = 0 ; i < n ; i++)
		mpfr_set_d(A[i*n + i], 1, MPFR_RNDN);
}

void copyVector(mpfr_t* dest, mpfr_t* src, int n) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_set(dest[i], src[i], MPFR_RNDN);
}

void copyMatrix(mpfr_t* dest, mpfr_t* src, int n) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_set(dest[i*n + j], src[i*n + j], MPFR_RNDN);
}

// A := A^^t (transpose)
void trans(mpfr_t* A, int n, mpfr_t tmp) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < i ; j++) {
			mpfr_set(tmp, A[i*n + j], MPFR_RNDN);
			mpfr_set(A[i*n + j], A[j*n + i], MPFR_RNDN);
			mpfr_set(A[j*n + i], tmp, MPFR_RNDN);
		}
}

// Compute the 2-norm of a vector x = xre + i * xim (complex vector)
// The result is stored in renv
void norm2Complex(mpfr_t renv, int n, mpfr_t* xre, mpfr_t* xim, mpfr_t tmp) {
	int i;
	mpfr_set_d(renv, 0, MPFR_RNDN);
	for(i = 0 ; i < n ; i++) {
		mpfr_mul(tmp, xre[i], xre[i], MPFR_RNDN);
		mpfr_add(renv, renv, tmp, MPFR_RNDN);
		mpfr_mul(tmp, xim[i], xim[i], MPFR_RNDN);
		mpfr_add(renv, renv, tmp, MPFR_RNDN);
		// renv += xre[i] * xre[i] + xim[i] * xim[i];
	}
	mpfr_sqrt(renv, renv, MPFR_RNDN);
}

// Compute the infinite norm of a vector v = vre + i * vim (complex vector)
// The result is stored in renv
void normInfVectComplex(mpfr_t renv, mpfr_t* vre, mpfr_t* vim, int n, mpfr_t tmp) {
	int i;
	mpfr_set_d(renv, -1, MPFR_RNDN);
	for(i = 0 ; i < n ; i++) {
		mpfr_hypot(tmp, vre[i], vim[i], MPFR_RNDN);
		if(mpfr_cmp(renv, tmp) < 0)
			mpfr_set(renv, tmp, MPFR_RNDN);
	}
}

// Compute the 2-norm of a matrix A (Real matrix)
// The result is stored in renv
void normInfMatrix(mpfr_t renv, mpfr_t* A, int n, mpfr_t tmp1, mpfr_t tmp2) {
	int i, j;
	mpfr_set_d(renv, -1, MPFR_RNDN);
	for(i = 0 ; i < n ; i++) {
		mpfr_set_d(tmp1, 0, MPFR_RNDN);
		for(j = 0 ; j < n ; j++) {
			mpfr_abs(tmp2, A[i*n + j], MPFR_RNDN);
			mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
			// tmp += fabs(A[i*n + j]);
		}
		if(mpfr_cmp(renv, tmp1) < 0)
			mpfr_set(renv, tmp1, MPFR_RNDN);
	}
}

// The Hessenberg reduction of the matrix A by successive multiplications of Householder matrix.
// The function overwrite A by his reduction
// Computation of the Orthogonal matrix U such as A := U^t * A * U
void toHessenbergForm(int n, mpfr_t* A, mpfr_t* U, mpfr_prec_t prec) {
	
	mpfr_t alpha, r, tmp_calc;
	int i, j, k, taille;
	mpfr_t* x = allocateVector(n); // Vector of the Householder Matrix
	mpfr_t* tmp_A = allocateVector(n); // work vector
	mpfr_t* tmp_U = allocateVector(n); // work vector
	
	mpfr_init2(alpha, prec);
	mpfr_init2(r, prec);
	mpfr_init2(tmp_calc, prec);
	
	initVector(x, n, prec);
	initVector(tmp_A, n, prec);
	initVector(tmp_U, n, prec);
	
	// Reduction of the differents columns
	for(i = 0 ; i < n - 2 ; i++) { 
		taille = n - i;
		
		// we're going to store in x the value of the wanted Vector of the Householder Matrix
		for(j = i ; j < n ; j++) // Copy of the part of the column that we want to reduce
			mpfr_set(x[j - i], A[j*n + i], MPFR_RNDN);

		// Set to zero the first element
		mpfr_set_d(*x, 0, MPFR_RNDN);
		
		// New element sub-diagonal
		norm2(alpha, taille, x, 1);
		mpfr_copysign(alpha, alpha, x[1], MPFR_RNDN);
		mpfr_mul_d(alpha, alpha, -1, MPFR_RNDN);
		//// alpha = -copysign(norm2(taille, x, 1), x[1])
		
		mpfr_mul(r, alpha, alpha, MPFR_RNDN);
		mpfr_mul(tmp_calc, alpha, x[1], MPFR_RNDN);
		mpfr_sub(r, r, tmp_calc, MPFR_RNDN);
		mpfr_div_d(r, r, 2, MPFR_RNDN);
		mpfr_sqrt(r, r, MPFR_RNDN);
		//// r = sqrt((alpha * alpha - x[1] * alpha) / 2)
		
		if(mpfr_cmp_ui(r, 0) == 0) continue; // if r == 0, we don't need to proceed (there is nothing to do)
		
		mpfr_mul_d(tmp_calc, r, 2, MPFR_RNDN); // tmp_calc := 2*r
		
		mpfr_sub(x[1], x[1], alpha, MPFR_RNDN);
		mpfr_div(x[1], x[1], tmp_calc, MPFR_RNDN);
		
		for(j = 2 ; j < taille ; j++) // Computation of the Vector of the Householder Matrix
			mpfr_div(x[j], x[j], tmp_calc, MPFR_RNDN);
		
		//// Computation of HA
		for(j = 0 ; j < taille ; j++) { // Computation of tmp_A := x^^t * A
			mpfr_set_d(tmp_A[j], 0, MPFR_RNDN);
			for(k = 0 ; k < taille ; k++) {
				mpfr_mul(tmp_calc, x[k], A[(i + k)*n + i + j], MPFR_RNDN);
				mpfr_add(tmp_A[j], tmp_A[j], tmp_calc, MPFR_RNDN);
			}
		}
		for(j = 0 ; j < taille ; j++) // Computation of A - 2*x*tmp_A = A - 2*x*(x^^t * A)
			for(k = 0 ; k < taille ; k++) {
				mpfr_mul(tmp_calc, tmp_A[k], x[j], MPFR_RNDN);
				mpfr_mul_d(tmp_calc, tmp_calc, 2, MPFR_RNDN);
				mpfr_sub(A[(i + j)*n + i + k], A[(i + j)*n + i + k], tmp_calc, MPFR_RNDN);
			}
		////
		//// Computation of (HA)H
		for(j = 0 ; j < n ; j++) { // Computation of tmp_A := (HA) * x
			mpfr_set_d(tmp_A[j], 0, MPFR_RNDN);
			mpfr_set_d(tmp_U[j], 0, MPFR_RNDN);
			for(k = 0 ; k < taille ; k++) {
				mpfr_mul(tmp_calc, x[k], A[j*n + i + k], MPFR_RNDN);
				mpfr_add(tmp_A[j], tmp_A[j], tmp_calc, MPFR_RNDN);
				// tmp_A[j] += x[k] * A[j*n + i + k];
				mpfr_mul(tmp_calc, x[k], U[j*n + i + k], MPFR_RNDN);
				mpfr_add(tmp_U[j], tmp_U[j], tmp_calc, MPFR_RNDN);
				// tmp_U[j] += x[k] * U[j*n + i + k];
			}
		}
		for(j = 0 ; j < n ; j++) // Computation of HAH = HA - 2 * tmp_A * x^^t = HA - 2 * ((HA) * x) * x^^t
			for(k = 0 ; k < taille ; k++) {
				mpfr_mul(tmp_calc, tmp_A[j], x[k], MPFR_RNDN);
				mpfr_mul_d(tmp_calc, tmp_calc, 2, MPFR_RNDN);
				mpfr_sub(A[j*n + i + k], A[j*n + i + k], tmp_calc, MPFR_RNDN);
				// A[j*n + i + k] -= 2 * x[k] * tmp_A[j];
				
				mpfr_mul(tmp_calc, tmp_U[j], x[k], MPFR_RNDN);
				mpfr_mul_d(tmp_calc, tmp_calc, 2, MPFR_RNDN);
				mpfr_sub(U[j*n + i + k], U[j*n + i + k], tmp_calc, MPFR_RNDN);
				// U[j*n + i + k] -= 2 * x[k] * tmp_U[j];
			}
		////
	}
	
	// Set to zero the elements under the sub-diagonal
	for(i = 2 ; i < n ; i++)
		for(j = 0 ; j < i - 1 ; j++)
			mpfr_set_d(A[i*n + j], 0, MPFR_RNDN);
	
	// Free the work variable
	mpfr_clear(alpha);
	mpfr_clear(r);
	mpfr_clear(tmp_calc);
	
	clearVector(x, n);
	clearVector(tmp_A, n);
	clearVector(tmp_U, n);
}

// Compute the coefficients of the Givens Matrix G which apply the following transformation :
// G * (A_ik) -> (x)
//     (A_jk)    (0)
void GivensCoefficients(int n, mpfr_t* A, int i, int j, int k, mpfr_t c, mpfr_t s, mpfr_t tmp) {
	mpfr_hypot(tmp, A[i*n + k], A[j*n + k], MPFR_RNDN);
	mpfr_div(c, A[i*n + k], tmp, MPFR_RNDN);
	mpfr_div(s, A[j*n + k], tmp, MPFR_RNDN);
	mpfr_mul_d(s, s, -1, MPFR_RNDN);
}

// Left multiplication by a Givens matrix of coefficient c and s to a Hessenberg Matrix A
// ! A must be a Hessenberg matrix
void GivensRotationGH(int n, mpfr_t* A, int i, int j, mpfr_t c, mpfr_t s, mpfr_t tmp1, mpfr_t tmp2) {
	int ind;
	for(ind = i ; ind < n ; ind++) {
		mpfr_set(tmp1, A[i*n + ind], MPFR_RNDN);
		mpfr_mul(tmp2, A[j*n + ind], s, MPFR_RNDN);
		mpfr_mul(A[i*n + ind], A[i*n + ind], c, MPFR_RNDN);
		mpfr_sub(A[i*n + ind], A[i*n + ind], tmp2, MPFR_RNDN);
		//// A[i*n + ind] = A[i*n + ind] * c - A[j*n + ind] * s
		
		mpfr_mul(tmp2, tmp1, s, MPFR_RNDN);
		mpfr_mul(A[j*n + ind], A[j*n + ind], c, MPFR_RNDN);
		mpfr_add(A[j*n + ind], A[j*n + ind], tmp2, MPFR_RNDN);
		//// A[j*n + ind] = tmp1 * s + A[j*n + ind] * c
	}
}

// Rigth multiplication by a Givens matrix of coefficient c and s to a Hessenberg Matrix A
// ! A must be a Hessenberg matrix
void GivensRotationDH(int n, mpfr_t* A, int i, int j, mpfr_t c, mpfr_t s, mpfr_t tmp1, mpfr_t tmp2) {
	int ind;
	for(ind = 0 ; ind < i + 2 ; ind++) {
		mpfr_set(tmp1, A[ind*n + i], MPFR_RNDN);
		mpfr_mul(tmp2, A[ind*n + j], s, MPFR_RNDN);
		mpfr_mul(A[ind*n + i], A[ind*n + i], c, MPFR_RNDN);
		mpfr_add(A[ind*n + i], A[ind*n + i], tmp2, MPFR_RNDN);
		//// A[ind*n + i] = A[ind*n + i] * c + A[ind*n + j] * s
		
		mpfr_mul(tmp2, tmp1, s, MPFR_RNDN);
		mpfr_mul(A[ind*n + j], A[ind*n + j], c, MPFR_RNDN);
		mpfr_sub(A[ind*n + j], A[ind*n + j], tmp2, MPFR_RNDN);
		//// A[ind*n + j] = -1 * tmp1 * s + A[ind*n + j] * c
	}
}

// Update of the sub-diagonal
// If an element of the sub-diagonal is small enough (|A_(i+1)i| <= tol * (|A_ii| + |A_(i + 1)(i + 1)|))
// he is explicitly set to zero
void subDiagonalUpdate(int n, mpfr_t* A, mpfr_t tol, mpfr_t tmp1, mpfr_t tmp2) {
	int i;
	for(i = 0 ; i < n - 1 ; i++) {
		mpfr_abs(tmp2, A[i*n + i], MPFR_RNDN);
		mpfr_abs(tmp1, A[(i + 1)*n + (i + 1)], MPFR_RNDN);
		mpfr_add(tmp2, tmp2, tmp1, MPFR_RNDN);
		mpfr_mul(tmp2, tmp2, tol, MPFR_RNDN);
		mpfr_abs(tmp1, A[(i+1)*n + i], MPFR_RNDN);
		if(mpfr_lessequal_p(tmp1, tmp2))
			mpfr_set_d(A[(i+1)*n + i], 0, MPFR_RNDN);
	}
}

// Compute p and q such as A is 2x2 block diagonal from the index 0 to p and the index q to n
// The iterations are to be applied only on the square Matrix from p to q
void deflation(int n, mpfr_t* A, int* p, int* q) {
	int i;
	for(i = *p ; i < n - 1 ; i++, (*p)++)
		if(mpfr_cmp_ui(A[(i+1)*n + i], 0) != 0 && (i == n - 2 || mpfr_cmp_ui(A[(i+2)*n + i + 1], 0) != 0))
			break;
	for(i = *q ; i > 0 ; i--, (*q)--)
		if(mpfr_cmp_ui(A[i*n + (i - 1)], 0) != 0 && (i == 1 || mpfr_cmp_ui(A[(i - 1)*n + (i - 2)], 0) != 0))
			break;
}

// One step of the QR algorithm with a simple shift of a Hessenberg Matrix A
void QRstepH(int n, mpfr_t* A, mpfr_t* rot, int p, int q, mpfr_t c, mpfr_t s, mpfr_t tmp1, mpfr_t tmp2, mpfr_t mu) {
	int i;
	
	mpfr_set(mu, A[q*n + q], MPFR_RNDN);
	
	// Computation of A - mu*I
	for(i = p ; i <= q ; i++)
		mpfr_sub(A[i*n + i], A[i*n + i], mu, MPFR_RNDN);
	
	// Computation of the QR decomposition
	// R is stored in A
	// Q is stored as the successives givens coefficients in rot
	for(i = p ; i <= q - 1 ; i++) {
		GivensCoefficients(n, A, i, i + 1, i, c, s, tmp1);
		GivensRotationGH(n, A, i, i + 1, c, s, tmp1, tmp2);
		mpfr_set(rot[2*i], c, MPFR_RNDN);
		mpfr_set(rot[2*i + 1], s, MPFR_RNDN);
	}
	// Computation of RQ
	for(i = p ; i <= q - 1 ; i++) {
		mpfr_mul_d(rot[2*i + 1], rot[2*i + 1], -1, MPFR_RNDN);
		GivensRotationDH(n, A, i, i + 1, rot[2*i], rot[2*i + 1], tmp1, tmp2);
	}
	
	// A + mu*I
	for(i = p ; i <= q ; i++)
		mpfr_add(A[i*n + i], A[i*n + i], mu, MPFR_RNDN);
}

// Iteration of the QR algorithm to a Hessenberg matrix A
// A is overwrited
// A should be at the end a 2x2 block diaganol matrix
void iterQR(int n, mpfr_t* A, mpfr_prec_t prec) {
	mpfr_t c, s, tmp1, tmp2, tol, mu;
	int i = 0;
	int p = 0, q = n-1;
	int tmpp = p, tmpq = q;
	mpfr_t* rot = allocateVector(2*n);
	
	mpfr_init2(c, prec);
	mpfr_init2(s, prec);
	mpfr_init2(tmp1, prec);
	mpfr_init2(tmp2, prec);
	mpfr_init2(tol, prec);
	mpfr_init2(mu, prec);
	
	// Tolerance before setting to zero
	// The smallest number x such as 1 + x != 1 (Compute according to the precision)
	mpfr_set_d(tol, 2, MPFR_RNDN);
	mpfr_pow_si(tol, tol, -1 * (prec - 1), MPFR_RNDN);
	
	initVector(rot, 2*n, prec);
	
	subDiagonalUpdate(n, A, tol, tmp1, tmp2);
	deflation(n, A, &p, &q); // Updating the p and q index
	while(p < q) {
		i++;
		QRstepH(n, A, rot, p, q, c, s, tmp1, tmp2, mu); // QR step
		subDiagonalUpdate(n, A, tol, tmp1, tmp2); // Updating sub-diagonal
		deflation(n, A, &p, &q); // Updating the p and q index
	}
	
	mpfr_clear(c);
	mpfr_clear(s);
	mpfr_clear(tmp1);
	mpfr_clear(tmp2);
	mpfr_clear(tol);
	mpfr_clear(mu);
	
	clearVector(rot, 2*n);
}

//// FUNCTION TO DEAL WITH COMPLEX MATRIX

//// Decomposition and resolution are done on complex matrix
//// Re is the matrix of the real parts and Im is the matrix of the imaginary parts

// !! For the addition and substraction, the output variable can be the same as the input variable
void addComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2) {
	mpfr_add(renvRe, re1, re2, MPFR_RNDN);
	mpfr_add(renvIm, im1, im2, MPFR_RNDN);
}

void subComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2) {
	mpfr_sub(renvRe, re1, re2, MPFR_RNDN);
	mpfr_sub(renvIm, im1, im2, MPFR_RNDN);
}

void mulComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2, mpfr_t tmp) {
	mpfr_mul(renvRe, re1, re2, MPFR_RNDN);
	mpfr_mul(tmp, im1, im2, MPFR_RNDN);
	mpfr_sub(renvRe, renvRe, tmp, MPFR_RNDN);
	// renvRe = re1 * re2  - im1 * im2;
	
	mpfr_mul(renvIm, re1, im2, MPFR_RNDN);
	mpfr_mul(tmp, re2, im1, MPFR_RNDN);
	mpfr_add(renvIm, renvIm, tmp, MPFR_RNDN);
	// renvIm = re1 * im2 + re2 * im1;
}

void divComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp) {
	mpfr_mul(tmp, re2, re2, MPFR_RNDN);
	mpfr_mul(tmpre, im2, im2, MPFR_RNDN);
	mpfr_add(tmp, tmp, tmpre, MPFR_RNDN);
	// tmp = re2 * re2 + im2 * im2;
	
	mpfr_div(tmpre, re2, tmp, MPFR_RNDN);
	mpfr_div(tmpim, im2, tmp, MPFR_RNDN);
	mpfr_mul_d(tmpim, tmpim, -1, MPFR_RNDN);
	// tmpre = re2 / tmp;
	// tmpim = -1 * im2 / tmp;
	
	mpfr_mul(renvRe, re1, tmpre, MPFR_RNDN);
	mpfr_mul(tmp, im1, tmpim, MPFR_RNDN);
	mpfr_sub(renvRe, renvRe, tmp, MPFR_RNDN);
	// renvRe = re1 * tmpre  - im1 * tmpim;
	
	mpfr_mul(renvIm, re1, tmpim, MPFR_RNDN);
	mpfr_mul(tmp, tmpre, im1, MPFR_RNDN);
	mpfr_add(renvIm, renvIm, tmp, MPFR_RNDN);
	// renvIm = re1 * tmpim + tmpre * im1;
}

// Compute the LU decomposition of a complex matrix H = Re + i * Im
// !!! H must be a Hessemberg matrix
void luDecompositionH(mpfr_t* Re, mpfr_t* Im, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp1, mpfr_t tmp2, mpfr_t tmp3) {
	int i, k;
	for(i = 0 ; i < n - 1 ; i++) {
		divComplex(tmpre, tmpim, Re[(i+1)*n + i], Im[(i+1)*n + i], Re[i*n + i], Im[i*n + i], tmp1, tmp2, tmp3);
		mpfr_set(Re[(i+1)*n + i], tmpre, MPFR_RNDN);
		mpfr_set(Im[(i+1)*n + i], tmpim, MPFR_RNDN);
		for(k = i + 1 ; k < n ; k++) {
			mulComplex(tmpre, tmpim, Re[(i+1)*n + i], Im[(i+1)*n + i], Re[i*n + k], Im[i*n + k], tmp1);
			subComplex(Re[(i+1)*n + k], Im[(i+1)*n + k], Re[(i+1)*n + k], Im[(i+1)*n + k], tmpre, tmpim);
		}
	}
}

// Resolve the system Ux = b with U = Re + i * Im an upper triangular complex matrix
// The solution x is stored in the complex vector b = bre + i * bim
void upperTriangularSolve(mpfr_t* Re, mpfr_t* Im, mpfr_t* bre, mpfr_t* bim, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp1, mpfr_t tmp2, mpfr_t tmp3) {
	int i, j;
	for(i = n - 1 ; i >= 0 ; i--) {
		for(j = n - 1 ; j > i ; j--) {
			mulComplex(tmpre, tmpim, Re[i*n + j], Im[i*n + j], bre[j], bim[j], tmp1);
			subComplex(bre[i], bim[i], bre[i], bim[i], tmpre, tmpim);
		}
		divComplex(tmpre, tmpim, bre[i], bim[i], Re[i*n + i], Im[i*n + i], tmp1, tmp2, tmp3);
		mpfr_set(bre[i], tmpre, MPFR_RNDN);
		mpfr_set(bim[i], tmpim, MPFR_RNDN);
	}
}

// Resolve the system Lx = b with L = Re + i * Im a lower triangular complex matrix
// The solution x is stored in the complex vector b = bre + i * bim
void lowerTriangularSolve(mpfr_t* Re, mpfr_t* Im, mpfr_t* bre, mpfr_t* bim, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < i ; j++) {
			mulComplex(tmpre, tmpim, Re[i*n + j], Im[i*n + j], bre[j], bim[j], tmp);
			subComplex(bre[i], bim[i], bre[i], bim[i], tmpre, tmpim);
		}
}

// Resolve the system Ax = b with A = LU. We must have already applied the luDecompositionH to A
// The solution x is stored in the complex vector b = bre + i * bim
void systemSolveLU(mpfr_t* Re, mpfr_t* Im, mpfr_t* bre, mpfr_t* bim, int n, mpfr_t tmp1, mpfr_t tmp2, mpfr_t tmp3, mpfr_t tmp4, mpfr_t tmp5) {
	lowerTriangularSolve(Re, Im, bre, bim, n, tmp1, tmp2, tmp3);
	upperTriangularSolve(Re, Im, bre, bim, n, tmp1, tmp2, tmp3, tmp4, tmp5);
}

// Computation of C = A*B
void mulMatrixComplexReal(mpfr_t* Cre, mpfr_t* Cim, mpfr_t* Are, mpfr_t* Aim, mpfr_t* Bre, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp1, mpfr_t tmp2) {
	int i, j, k;
	mpfr_set_d(tmp1, 0, MPFR_RNDN);
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++) {
			mpfr_set_d(Cre[i*n + j], 0, MPFR_RNDN);
			mpfr_set_d(Cim[i*n + j], 0, MPFR_RNDN);
			for(k = 0 ; k < n ; k++) {
				mulComplex(tmpre, tmpim, Are[i*n + k], Aim[i*n + k], Bre[k*n + j], tmp1, tmp2);
				addComplex(Cre[i*n + j], Cim[i*n + j], Cre[i*n + j], Cim[i*n + j], tmpre, tmpim);
			}
		}
}
////

// Find the Eigen values of A. A must have been through the QR algorithm
// re contains the real part of the Eigen values and im the imaginary part 
void findEigenValues(int n, mpfr_t* A, mpfr_t* re, mpfr_t* im, mpfr_prec_t prec) {
	int i = 0;
	mpfr_t a, b, c, d, tmp, disc;
	
	mpfr_init2(a, prec);
	mpfr_init2(b, prec);
	mpfr_init2(c, prec);
	mpfr_init2(d, prec);
	mpfr_init2(tmp, prec);
	mpfr_init2(disc, prec);
	
	while(i < n) {
		if(i != n - 1 && mpfr_cmp_ui(A[(i+1)*n + i], 0) != 0) { // if we find a 2x2 diagonal block, we compute the roots of his characteristic polynomial
			mpfr_set(a, A[i*n + i], MPFR_RNDN);
			mpfr_set(b, A[i*n + i + 1], MPFR_RNDN);
			mpfr_set(c, A[(i+1)*n + i], MPFR_RNDN);
			mpfr_set(d, A[(i+1)*n + i + 1], MPFR_RNDN);
			
			mpfr_mul(disc, c, b, MPFR_RNDN);
			mpfr_mul(tmp, a, d, MPFR_RNDN);
			mpfr_sub(disc, tmp, disc, MPFR_RNDN);
			mpfr_mul_d(disc, disc, 4, MPFR_RNDN);
			mpfr_add(tmp, a, d, MPFR_RNDN);
			mpfr_pow_ui(tmp, tmp, 2, MPFR_RNDN);
			mpfr_sub(tmp, tmp, disc, MPFR_RNDN);
			
			if(mpfr_cmp_ui(tmp, 0) < 0) {
				mpfr_mul_d(tmp, tmp, -1, MPFR_RNDN);
				mpfr_sqrt(disc, tmp, MPFR_RNDN);
				
				mpfr_add(re[i], a, d, MPFR_RNDN);
				mpfr_div_d(re[i], re[i], 2, MPFR_RNDN);
				
				mpfr_set(re[i + 1], re[i], MPFR_RNDN);
				
				mpfr_div_d(im[i], disc, 2, MPFR_RNDN);
				mpfr_mul_d(im[i + 1], im[i], -1, MPFR_RNDN);

				i += 2;
			}
			else {
				mpfr_sqrt(disc, tmp, MPFR_RNDN);
				
				mpfr_add(re[i], a, d, MPFR_RNDN);
				mpfr_set(re[i + 1], re[i], MPFR_RNDN);
				
				mpfr_add(re[i], re[i], disc, MPFR_RNDN);
				mpfr_sub(re[i + 1], re[i + 1], disc, MPFR_RNDN);
				
				mpfr_div_d(re[i], re[i], 2, MPFR_RNDN);
				mpfr_div_d(re[i + 1], re[i + 1], 2, MPFR_RNDN);

				mpfr_set_d(im[i], 0, MPFR_RNDN);
				mpfr_set_d(im[i + 1], 0, MPFR_RNDN);

				i += 2;
			}
		}
		else { // we read the eigen value on the diagonal coefficient
			mpfr_set(re[i], A[i*n + i], MPFR_RNDN);
			mpfr_set_d(im[i], 0, MPFR_RNDN);
			i++;
		}
	}
	
	mpfr_clear(a);
	mpfr_clear(b);
	mpfr_clear(c);
	mpfr_clear(d);
	mpfr_clear(tmp);
	mpfr_clear(disc);
}

// Computation of the Eigen vectors associated with the complex Eigen values val = re + i * im
void findEigenVectors(int n, mpfr_t* A, mpfr_t* eigVecre, mpfr_t* eigVecim, mpfr_t* re, mpfr_t* im, mpfr_prec_t prec) {
	int i, j, compt, essai;
	
	mpfr_t tol, tmp1, tmp2, tmp3, tmp4, tmp5;
	
	mpfr_t* Are = allocateMatrix(n);
	mpfr_t* Aim = allocateMatrix(n);
	
	mpfr_t* tmpre = allocateVector(n);
	mpfr_t* tmpim = allocateVector(n);
	
	mpfr_init2(tol, prec);
	mpfr_init2(tmp1, prec);
	mpfr_init2(tmp2, prec);
	mpfr_init2(tmp3, prec);
	mpfr_init2(tmp4, prec);
	mpfr_init2(tmp5, prec);
	
	initMatrix(Are, n, prec);
	initMatrix(Aim, n, prec);
	
	initVector(tmpre, n, prec);
	initVector(tmpim, n, prec);
	
	// Computation of the tolerance of the computation
	mpfr_set_d(tol, 2, MPFR_RNDN);
	mpfr_pow_si(tol, tol, -1 * (prec - 3), MPFR_RNDN);
	normInfMatrix(tmp1, A, n, tmp2, tmp3);
	mpfr_mul(tol, tol, tmp1, MPFR_RNDN);
	// tol = 2^(3 - prec) * normInfMatrix(A, n);
	
	// Computation of each eigen vector by an inverse iteration
	// We start with a random vector
	// We apply to it an inverse iteration
	// If the iteration converges, it requires only few iterations (1 or 2)
	// Else, we try again with an another random vector until it converges
	for(i = 0 ; i < n ; i++) {
		essai = 0;
		
		// We copy A to proceed the computation on it
		copyMatrix(Are, A, n);
		set0Matrix(Aim, n);
		for(j = 0 ; j < n ; j++) {
			mpfr_sub(Are[j*n + j], Are[j*n + j], re[i], MPFR_RNDN);
			// Are[j*n + j] -= re[i];
			mpfr_sub(Aim[j*n + j], Aim[j*n + j], im[i], MPFR_RNDN);
			// Aim[j*n + j] -= im[i];
		}

		// LU decomposition
		luDecompositionH(Are, Aim, n, tmp1, tmp2, tmp3, tmp4, tmp5);
		do {
			compt = 0;
			essai++;
			
			// We choose a complex or real random vector depending on wether the eigen value is complex or real
			setAleaVector(eigVecre + i*n, n, tmp1);
			if(mpfr_cmp_d(im[i], 0) == 0)
				set0Vector(eigVecim + i*n, n);
			else
				setAleaVector(eigVecim + i*n, n, tmp1);
			
			// Iteration
			do {
				compt++;
				copyVector(tmpre, eigVecre + i*n, n);
				copyVector(tmpim, eigVecim + i*n, n);

				systemSolveLU(Are, Aim, eigVecre + i*n, eigVecim + i*n, n, tmp1, tmp2, tmp3, tmp4, tmp5);
				
				norm2Complex(tmp1, n, eigVecre + i*n, eigVecim + i*n, tmp2);

				for(j = 0 ; j < n ; j++) {
					mpfr_div(eigVecre[i*n + j], eigVecre[i*n + j], tmp1, MPFR_RNDN);
					// eigVecre[i*n + j] /= tmp;
					mpfr_div(eigVecim[i*n + j], eigVecim[i*n + j], tmp1, MPFR_RNDN);
					// eigVecim[i*n + j] /= tmp;
					mpfr_div(tmpre[j], tmpre[j], tmp1, MPFR_RNDN);
					// tmpre[j] /= tmp;
					mpfr_div(tmpim[j], tmpim[j], tmp1, MPFR_RNDN);
					// tmpim[j] /= tmp;
				}
				normInfVectComplex(tmp1, tmpre, tmpim, n, tmp2);
			} while(mpfr_cmp(tmp1, tol) > 0 && compt <= 3); // End if the eigen vector has already converged or if there are too many iterations
		}  while(compt > 3 && essai <= n*n); // End if there are not too many iterations or if we made too many try (to avoid the program to block)
		if(mpfr_cmp_d(im[i], 0) != 0) { // If we were computing an complex eigen vector, create also its conjugate eigen vector for the conjugate eigen value
			copyVector(eigVecre + (i + 1)*n, eigVecre + i*n, n);
			for(j = 0 ; j < n ; j++)
				mpfr_mul_d(eigVecim[(i + 1)*n + j], eigVecim[i*n + j], -1, MPFR_RNDN);
			i++;
		}
		
		if(essai > n*n) {
			if(mpfr_cmp_d(im[i], 0) != 0)
				i += 2;
			else
				i++;
//			printf("Calcul du Vecteur propre de la valeur propre numéro %d n'a pas aboutie au bout de %d essais, sauté\n", i, essai);
		}
//		printf("Vecteur propre de la valeur propre numéro %d calculé en %d essais\n", i, essai);
	}
	
	mpfr_clear(tol);
	mpfr_clear(tmp1);
	mpfr_clear(tmp2);
	mpfr_clear(tmp3);
	mpfr_clear(tmp4);
	mpfr_clear(tmp5);
	
	clearMatrix(Are, n);
	clearMatrix(Aim, n);
	
	clearVector(tmpre, n);
	clearVector(tmpim, n);
}

// Computation of the eigen values of matrix and, if ask_for_Evector is not equal to 0, the eigen vectors
void spectrum(int n, mpfr_t* A, mpfr_t* re, mpfr_t* im, int ask_for_Evector, mpfr_t* evre, mpfr_t* evim, mpfr_prec_t prec) {
	srand(time(NULL));
	
	// Temporary variables
	mpfr_t tmp1, tmp2, tmp3, tmp4;
	mpfr_init2(tmp1, prec);
	mpfr_init2(tmp2, prec);
	mpfr_init2(tmp3, prec);
	mpfr_init2(tmp4, prec);
	//
	
	// Matrix of the transformation of the Hessenberg reduction
	mpfr_t* U = allocateMatrix(n);
	initMatrix(U, n, prec);
	identityMatrix(U, n);

	mpfr_t* tmpevre = allocateMatrix(n);
	mpfr_t* tmpevim = allocateMatrix(n);
	mpfr_t* H = allocateMatrix(n);
	
	initMatrix(tmpevre, n, prec);
	initMatrix(tmpevim, n, prec);
	initMatrix(H, n, prec);
	
	toHessenbergForm(n, A, U, prec);
	trans(U, n, tmp1);
	
	copyMatrix(H, A, n); // We store the Hessenberg reduction so that we could use it later
	
	iterQR(n, A, prec); // Apply QR algorithm to A
	
	findEigenValues(n, A, re, im, prec); // Reading of the eigen values
	
	if(ask_for_Evector) { // If we asked for eigen vectors
		findEigenVectors(n, H, tmpevre, tmpevim, re, im, prec); // Computation of the eigen vectors of H (Hessenberg reduction of A)
		mulMatrixComplexReal(evre, evim, tmpevre, tmpevim, U, n, tmp1, tmp2, tmp3, tmp4); // Multiplication by U to find the eigen vectors of A
		
		// The vectors are organized by lines
		// We are going to organize them in columns
		trans(evre, n, tmp1);
		trans(evim, n, tmp1);
	}
	
	mpfr_clear(tmp1);
	mpfr_clear(tmp2);
	mpfr_clear(tmp3);
	mpfr_clear(tmp4);
	
	clearMatrix(U, n);
	clearMatrix(tmpevre, n);
	clearMatrix(tmpevim, n);
	clearMatrix(H, n);
}














