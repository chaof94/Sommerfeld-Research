syms ba bb ua ub kx ky xp yp zp h 
A = [ 1 -1 ; -ba*ub/(bb*ua) -1];
c = 1i*ua/2*exp(-1j*(kx*xp + ky*yp))*exp(-1i*ba*(zp+h))/ba;
B = [1 ; ba*ub/(bb*ua)];
X = A\B;
simplify(X)
clc
pretty(X)
