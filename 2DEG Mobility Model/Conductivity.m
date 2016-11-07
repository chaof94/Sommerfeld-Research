%% Calculate 2DEG condunctivity

%% \sigma = N e \mu
% Sigma = Electrical Conductivity
% e = Electron Charge
% mu = Electron Mobility

% From [1] 
% [1]  Computational model of 2DEG mobility in AlGaN/GaN heterostructures
% Karine Abgaryan*, Ilya Mutigullin, and Dmitry Reviznikov
% Phys. Status Solidi C 12, No. 4–5, 460–465 (2015) / DOI 10.1002/pssc.201400200

%% From Allen's paper
c = 3e10; %cm/s
ns = 1e12 *1e4 ; %cm^-2
eps_0 = 8.854e-12;
eps_s = 12*eps_0;
eps_ox = 4*eps_0;
a = 3.52e-6;
d = .14e-6;
q = 2*pi/a;
m_e = 9.31e-31;
m_star = .2*m_e;
e = 1.602e-19;
N = (0:.01:4)*ns;
omega_p = sqrt(N*e^2/m_star * q/(eps_ox*(eps_ox + eps_s*coth(q*d))));
yyaxis left
plot(N,omega_p)
ylabel('\omega_p [rad/s]')
yyaxis right
k = omega_p/c;
plot(N,k)
ylabel('k [/cm]')
