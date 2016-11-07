#ifndef _DEFS_H
#define _DEFS_H

#define PI 3.141592653589793238462643383279

// Cartesian coordinates (x,y,z)
typedef struct{
	double x;
	double y;
	double z;
}cordC;
// Spherical coordinates (r, phi, theta)
typedef struct{
	double r;
	double ph;
	double th;
}cordS;



#endif
