#ifndef COMPLEX_H
#define COMPLEX_H

#include <math.h>
#include "defs.h"
//#define PI 3.141592653589793238462643

enum complex_format { REAL_IMAG, MAG_PHASE_DEG, MAG_PHASE_RAD };
class complex
{
    public: 

        double re, im;
        inline complex( void );
        //  complex( const double v1, const double v2 = 0.0,
        //          complex_format the_format = REAL_IMAG ); 

        complex( const double v1, const double v2 = 0.0);
        inline complex( const complex &z );
        inline double SetReal( const double the_real ); // set & return the real part
        inline double SetImag( const double the_imag ); // set & return the imag part
        inline double Real( void ) const;     // return the real part
        inline double Imag( void ) const;     // return the imaginary part
        inline double Norm( void ) const;     // return magnitude squared
        inline double Mag( void ) const;      // return magnitude
        inline double Phase( void ) const;    // returns phase angle in radians

        inline double Degree( void ) const;   // returns phase angle in degree
        inline complex operator -( void );
        inline complex &operator =( const complex &z );
        inline complex &operator =( const double &z );
                                 // assignment operator avoids bitwise copy
        inline complex &operator+=( const double b );
        inline complex &operator+=( const complex &z );
        inline complex &operator-=( const double b );
        inline complex &operator-=( const complex &z );
        inline complex &operator*=( const double b );
        inline complex &operator/=( const double b );
        complex &operator*=( const complex &z );
        complex &operator/=( const complex &z );
        friend inline short operator == ( const complex &a, const complex &b);
        friend inline short operator ==( const complex &a, const double b );
        friend inline short operator ==( const double a, const complex &b );
        friend inline short operator !=( const complex &a, const complex &b ) ;
        friend inline short operator !=( const complex &a, const double b );
        friend inline short operator !=( const double a, const complex &b );
        friend inline complex operator +( const complex &a, const complex &b );

        friend inline complex operator +( const complex &a, const double b );
        friend inline complex operator +( const double a, const complex &b );
        friend inline complex operator -( const complex &a, const complex &b );

        friend inline complex operator -( const complex &a, const double b );
        friend inline complex operator -( const double a, const complex &b );
        friend complex operator *( const complex &a, const complex &b );
        friend inline complex operator *( const complex &a, const double b );
        friend inline complex operator *( const double a, const complex &b );
        friend complex operator /( const complex &a, const complex &b );
        friend complex operator /( const complex &a, const double b );
        friend complex operator /( const double a, const complex &b );
        friend inline double abs( const complex &a );  // absolute value ( == mag )
        friend inline double fabs( const complex &a ); // absolute value ( == mag )
        friend inline complex Conj( const complex &a ); // complex conjugate
        inline complex operator ~( void ) { return Conj(*this); }
        //add bb to aa
        friend inline void ComplexAddOn( complex &aa, const complex &bb ); 
        friend void ComplexMult( complex &aa, const complex &bb, const complex &cc ); // aa = bb * cc
        friend inline complex cos( const complex &x );
        friend inline complex sin( const complex &x );
        friend inline complex cosh( const complex &x );
        friend inline complex sinh( const complex &x );
        friend inline complex exp( const complex &x );
        friend complex log( const complex &x );
        friend complex log10( const complex &x );
        friend complex pow( const complex &x, double d );
        friend complex pow( const complex &x, const complex &p );
        friend complex sqrt( const complex &x );
        friend inline double dB( const complex &x, double multiplier);
};


inline complex::complex( )
{
    re = im = 0.0;
}


inline complex::complex( const complex &c )
{
    re = c.re;
    im = c.im;
}


inline complex &complex::operator =( const complex &x )
{
    re = x.re;
    im = x.im;
    return( *this );
}


inline complex &complex::operator =( const double &x )
{
    re = x;
    im = 0.0;
    return( *this );
}


inline double complex::SetReal( const double the_real )
{
    return( re =the_real );
}



inline double complex::SetImag( const double the_imag )
{
    return( im = the_imag );
}



inline double complex::Real( void ) const
{
    return( re );
}


inline double complex::Imag( void ) const
{
    return( im );
}



inline double complex::Norm( void ) const
{
    return( ( re * re ) + ( im * im ) );
}



inline double complex::Mag( void ) const
{
    return( hypot( re, im ) );
}



inline double complex::Phase( void ) const
{
    return( atan2( im, re ) ); // atan of ( im/re ) in ( -PI, PI )
}



inline double complex::Degree( void ) const
{
    return( atan2( im, re )*180/PI ); // Phase() * 180 / PI
}



inline complex complex::operator -( void )
{
    return( complex( -re, -im ) );
}



inline complex &complex::operator +=( const complex &b )
{
    re += b.re;
    im += b.im;
    return( *this );
}



inline complex &complex::operator +=( const double b )
{
    re += b;
    return( *this );
}



inline complex &complex::operator -=( const complex &b )
{
    re -= b.re;
    im -= b.im;
    return( *this );
}



inline complex &complex::operator -=( const double b )
{
    re -= b;
    return( *this );
}




inline complex &complex::operator *=( const double b )
{
    re *= b;
    im *= b;
    return( *this );
}




inline complex &complex::operator /=( const double b )
{
    re /= b;
    im /= b;
    return( *this );
}



inline double abs( const complex &a )
{
    return( a.Mag( ) );
}



inline double fabs( const complex &a )
{
    return( a.Mag( ) );
}



inline complex Conj( const complex &a )
{
    return( complex( a.re, -a.im ) );
}



inline short operator ==( const complex &a, const complex &b )
{
    return( ( a.re == b.re ) && ( a.im == b.im ) );
}



inline short operator ==( const complex &a, const double b )
{
    return( ( a.im == 0.0 ) && ( a.re == b ) );
}



inline short operator ==( const double a, const complex &b )
{
    return( ( b.im == 0.0 ) && ( a == b.re ) );
}



inline short operator !=( const complex &a, const complex &b )
{
    return( ( a.re != b.re ) || ( a.im != b.im ) );
}



inline short operator !=( const complex &a, const double b )
{
    return( ( a.im != 0.0 ) || ( a.re != b ) );
}



inline short operator !=( const double a, const complex &b )
{
    return( ( b.im != 0.0 ) || ( a != b.re ) );
}



inline complex operator +( const complex &a, const complex &b )
{
    return( complex( a.re + b.re, a.im + b.im ) );
}


    
inline complex operator +( const complex &a, const double b )
{
    return( complex( a.re + b, a.im ) );
}



inline complex operator +( const double a, const complex &b )
{
    return( complex( a + b.re, b.im ) );
}



inline complex operator -( const complex &a, const complex &b )
{
    return( complex( a.re - b.re, a.im - b.im ) );
}



inline complex operator -( const complex &a, const double b )
{
    return( complex( a.re - b, a.im ) );
}



inline complex operator -( const double a, const complex &b )
{
    return( complex( a - b.re, - b.im ) );
}



inline complex operator *( const double a, const complex &b )
{
    return( complex( a * b.re, a * b.im ) );
}



inline complex operator *( const complex &a, const double b )
{
    return( complex( a.re * b, a.im * b ) );
}



inline void ComplexAddOn( complex &aa, const complex &bb )
{
    aa.re += bb.re;
    aa.im += bb.im;
}



inline complex exp( const complex &z )
{
    register double r = exp( z.re );
    return( complex( r * cos( z.im ), r * sin( z.im ) ) );
}



inline complex cosh( const complex &z )
{
    return( complex( cos( z.im ) * cosh( z.re ), sin( z.im ) * sinh( z.re )));
}



inline complex sinh( const complex &z )
{
    return( complex( cos( z.im ) * sinh( z.re ), sin( z.im ) * cosh( z.re )));
}



inline complex cos( const complex &z )
{
    return( complex( cos( z.re ) * cosh( z.im ), -sin( z.re ) * sinh( z.im)));
}



inline complex sin( const complex &z )
{
    return( complex( sin( z.re ) * cosh( z.im ), cos( z.re ) * sinh( z.im )));
}


inline double dB( const complex &z, double multiplier )
{	
	multiplier =1.0;
    return (multiplier * 10.0*log10(z.Norm()));
}

#endif // COMPLEX_H

