clear all;clc;syms beta_a beta_b beta_c k_x k_y x_p y_p z_p h_1 h_2
z_p =0
A = [ 1,         -1,           -exp(-1i*beta_b*(h_1 + h_2)),         0; 
    beta_a/1,    beta_b/1,       -beta_b/1*exp(-1i*beta_b*(h_1 + h_2)),      0; ...
      0,    -exp(-1i*beta_b*(h_1 + h_2)),        -1,                 1;...
      0,  -beta_b/1*exp(-1i*beta_b*(h_1 + h_2)),   beta_b/1,          beta_c/1];
c = -1i*1/2*exp(-1j*(k_x*x_p + k_y*y_p))/beta_b;
B = [exp(-1i*beta_b*(h_1 + z_p));...
     exp(-1i*beta_b*(h_1 + z_p))*(beta_b/1);...
     exp(-1i*beta_b*(h_2 - z_p));...
     exp(-1i*beta_b*(h_2 - z_p))*(beta_b/1)];
X = A\B;
% simplify(X)
% clc
% pretty(X)
Ba = latex(X(1))
Bb1 = latex(X(2))
Bb2 = latex(X(3))
Bbc = latex(X(4))
