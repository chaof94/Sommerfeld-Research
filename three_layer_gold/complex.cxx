#include <math.h>
#include "complex.h"


complex::complex(const double v1,const double v2)
{
    re = v1;
    im = v2;
}


complex &complex::operator *= (const complex& b) {
    double d = re*b.re - im*b.im;
	im = re*b.im + im*b.re;
	re = d;
    return (*this);
} 

complex &complex::operator /= (const complex& b) {
    complex c = Conj(b) * (*this);
    double d = b.Norm();
    re = c.re/d;
    im = c.im/d;
    return (*this);
}

complex operator / (const complex&a, const complex& b) {
	complex c = a;
	c /= b;
	return c;
} 

complex operator / (const complex&a, const double b) {
    complex c;
    c.re = a.re/b;
    c.im = a.im/b;
    return c;
} 

complex operator / (const double a,  const complex& b) {
	complex c = 1.0;
	c /= b;
	return c*a;
}

complex operator * (const complex&a, const complex& b) {
	complex c = a;
	c *= b;
	return c;
} 

complex log( const complex& x) {
    double mag = x.Norm();
    double phi = x.Phase();
    mag = log(mag)*0.5;
    complex b = complex(mag,phi);
    return b;
}

complex log10( const complex& x) {
    double mag = x.Norm();
    double phi = x.Phase();
    mag = log10(mag)*0.5;
    phi = phi*0.4342944; //log10(e)
    complex b = complex(mag,phi);
    return b;
}


complex sqrt( const complex& x) {
    double mag = x.Mag();
    double phi = 0.5 * x.Phase();
    mag = sqrt(mag);
    complex b = complex(cos(phi),sin(phi));
    return b*mag;
} 

complex pow( const complex& x, double d) {
    double mag = x.Norm();
    double phi = d * x.Phase();
    mag = pow(mag,d*0.5);
    complex b = complex(cos(phi),sin(phi));
    return b*mag;
}


