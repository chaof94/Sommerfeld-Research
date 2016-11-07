% Main routine
%
% 2D Scattering project
% by Beat Hangartner
tic
parameters
g_eom
excitation
solver_mat2
toc
figure(1)
contourf(S.*conj(S))
title(['Scattering pattern lamda=',num2str(lamda),' Energy=',num2str(Energy)])
ylabel(['length = ',num2str(len)])
xlabel(['width = ',num2str(width)])
% colorbar
axis image
mov(1,10,1)