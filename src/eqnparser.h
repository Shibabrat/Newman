// eqnparser.h
//
// Header file for eqnparser module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#ifndef _EQNPARSER_H_
#define _EQNPARSER_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
using namespace std;
#include "macros.h"

typedef double type;  


extern char AnalyticEq[ND][200];

//Function prototypes
void next();
void syntax();
void unknown(char *s);
type variable(double t, double *X);
type constant(double t, double *X);
type function(double t, double *X);
type term(double t, double *X);
type H(double t, double *X);
type G(double t, double *X);
type F(double t, double *X);
type E(double t, double *X);
void GetEquationVelocity(double t, double *X, double *dXdt);

static inline type factorial(type v)
{
	type i, r = 1;
	for(i=2;i<=v;i++) r *= i;
	return r;
}

#endif


