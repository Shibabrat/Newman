// eqnparser.cpp
//
// eqnparser module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include "eqnparser.h"

char c;
char *p;

void GetEquationVelocity(double t, double *X, double *dXdt)
{
  for(int d=0;d<ND;++d) {
    p=AnalyticEq[d];
    next();
    dXdt[d] = E(t,X);
    if(c != 0) 
      syntax();
  }
}

void next()
{
	do { 
    c = *p++;
  }
	while(c == ' ');
}

void syntax()
{
	cout << "Syntax error in analytic equation" << endl;		
	exit(1);
}


void unknown(char *s)
{
	cout << s << " : unknown variable or function" << endl;
	exit(1);	
}

type variable(double t, double *X)
{
  if(c == 't') { 
    next();
    return t;
  }
  
#if ND==2
  if(c == 'x') { 
    next();
    return X[0];
  }
  if(c == 'y') { 
    next();
    return X[1];
  }
#elif ND==3
  if(c == 'x') { 
    next();
    return X[0];
  }
  if(c == 'y') { 
    next();
    return X[1];
  }
  if(c == 'z') { 
    next();
    return X[2];
  }
#elif ND==4
   if(c == 'x') { 
    next();
    return X[0];
  }
  if(c == 'y') { 
    next();
    return X[1];
  }
  if(c == 'z') { 
    next();
    return X[2];
  }
  if(c == 'w') { 
    next();
    return X[3];
  }
#endif		
  return 0.0;
}


type constant(double t, double *X)
{
	type r = 0;
	
	while(c >= '0' && c <= '9') 
		{ r = 10*r + (c-'0'); next(); }
	
	if(c == '.')
		{
    type p = 0.1;
    next();
    
    while(c >= '0' && c <= '9') 
				{ r += p * (c-'0'); p /= 10; next(); }
		}
	
	if(c == 'e' || c == 'E')
		{ 
    type m = 1;			
    
    next();
		  if(c == '-') { m = -m; next(); } else if(c == '+') next();
      
			r *= pow(10,m*term(t,X));
		}
		
	return r;
}


type function(double t, double *X)
{
	char f[20], *p;
	type v;
	
	p = f;
	while(p-f < 19 && c >= 'a' && c <= 'z') { *p++ = c; next(); }
	
	*p = 0;
	
	if(!strcmp(f,"pi")) return M_PI;
	if(!strcmp(f,"e" )) return M_E;
  
	v = term(t,X);
	
#define mathfunc(a,b) if(!strcmp(f,a)) return b;
  
	mathfunc("abs"   , fabs(v));
	mathfunc("fabs"  , fabs(v));
	mathfunc("floor" , floor(v));
	mathfunc("ceil"  , ceil(v));
	mathfunc("sqrt"  , sqrt(v));
	mathfunc("exp"   , exp(v));
  
	mathfunc("sin"   , sin(v));
	mathfunc("cos"   , cos(v));
	mathfunc("tan"   , tan(v));
	mathfunc("asin"  , asin(v));
	mathfunc("acos"  , acos(v));
	mathfunc("atan"  , atan(v));
  
	mathfunc("sinh"  , sinh(v));
	mathfunc("cosh"  , cosh(v));
	mathfunc("tanh"  , tanh(v));
	mathfunc("asinh" , asinh(v));
	mathfunc("acosh" , acosh(v));
	mathfunc("atanh" , atanh(v));
  
  unknown(f);
  return 0;
}

type term(double t, double *X)
{
	const int m = 1;
  
	if(c=='(' || c=='[') {
    type r;
    
    next();
    r = E(t,X);
    if(c != ')' && c !=']') 
      syntax();
    
    next();
    return m * r;
  }
	
	else if((c >= '0' && c <= '9') || c == '.') {
    return m * constant(t,X);
  
  }
  else if( c=='t' || c=='x' || c=='y' || c=='z' )
    return m * variable(t,X);
  
	else if(c >= 'a' && c <= 'z')
    return m * function(t,X);
  
  return 0.0;
}


type H(double t,double *X)
{
	type r = term(t,X);
	
	if(c == '!') { next(); r = factorial(r); }
	return r;
}

type G(double t, double *X)
{
	type q, r = H(t,X);
	
	while(c == '^')
		{ next(); q = G(t,X); r = pow(r,q); }
	
	return r;
}

type F(double t, double *X)
{
	type r = G(t,X);
  
	while(1)
		{
    if(c=='*') { next(); r *= G(t,X); }
    else
			if(c=='/') { next(); r /= G(t,X); }
    else
			if(c=='%') { next(); r = fmod(r,G(t,X)); }
    else break;
		}
	
	return r;
}

type E(double t, double *X)
{
	type r = F(t,X);
	while(1) {
    if(c=='+') { next(); r += F(t,X); }
    else
			if(c=='-') { next(); r -= F(t,X); }
    else break;
  }
	
	return r;
}


