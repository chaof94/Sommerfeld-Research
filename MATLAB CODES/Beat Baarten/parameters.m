% Parameters definition
%
% 2D Scattering
% by Beat Hangartner
%
% History
% 26.5.02 created
% 29.5.02 modified Filname format
%
% Dimensions follow the SI-System
global I dx dy gamma beta kx ky k0 width len geom_size R_i_eff V_i
% global filename
% imaginary unit
I = sqrt(-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excitation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wavelength
if exist('lamda')==0
    global lamda
    lamda=1;
end;
%Electrical field strength
E0=1;
%Wave propagation direction
kx_=1;
ky_=0;
kz=0; %(2D)
kx=kx_/sqrt(kx_^2+ky_^2)*2*pi/lamda;
ky=ky_/sqrt(kx_^2+ky_^2)*2*pi/lamda;
k0=2*pi/lamda;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Background dielectric constant
epsilon_b=1;
%Grid spaceing
dx = 1;
dy = dx;
%material file
%filename='020-015-block1';
filename='020-010-block1';
%filename='002-002-block1';
%filename='010-010-block1';
%filename='080-040-block1';
%filename='080-040-grid';
%Area
width=str2num(filename(1:3));
len=str2num(filename(5:7));
geom_size = len*width/dx/dy;
% mesh Volume (Area)
V_i = dx*dy;
% effective radius
R_i_eff = sqrt(V_i/pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=1-kz^2/k0^2;
gamma=R_i_eff/k0*besselh(1,1,R_i_eff*k0)+I*2/(pi*k0^2);