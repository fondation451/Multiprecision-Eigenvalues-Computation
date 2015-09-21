#ifndef H_FONCTION_QR_H
#define H_FONCTION_QR_H

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <mpfr.h>


////////
//// Compute the Eigen values and, if asked, the Eigen vectors of an nxn square Matrix A.
//// A should be a linearised nxn matrix of mpfr_t floating point.
//// re and im should be vector of size n with mpfr_t variables fully initialized.
//// evre et evim should be nxn matrix with mpfr_t variables fully initialized.
//// ask_for_Evector should be 0 or 1. If it is 0, the eigen vectors are not be compute.
//// prec is the precision of the computation. It should the same as the precision of the different matrix/vectors.

//// At the end, the n eigen values are stored in re and im such as re is the real part and im is the imaginary part.
//// If asked, the eigen vectors are stored in evre and evim by column. The first column correspond to the first eigen value in re and im...
//// evre is the real part and evim is the imaginary part.
//// !!!! A is overwrited !!!!
void spectrum(int n, mpfr_t* A, mpfr_t* re, mpfr_t* im, int ask_for_Evector, mpfr_t* evre, mpfr_t* evim, mpfr_prec_t prec);

#endif






















