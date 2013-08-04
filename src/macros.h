//macros.h
#ifndef _MACROS_H_
#define _MACROS_H_

#define MAXSTP 10000

#define TINY_FTLE (1.0e-6)
#define TINY_DATA (1.0e-10)
#define TRUE 1
#define FALSE 0
#define FAILURE 0
#define SUCCESS 1
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHORTSTRING 100
#define LONGSTRING  200


#endif

#ifndef M_PI
#define M_PI (3.141592653589793238)
#endif
#ifndef HUGE_VAL
#define HUGE_VAL (99999999.0)
#endif

