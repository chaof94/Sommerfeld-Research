% Excitation array
%
% 2D Scattering
% by Beat Hangartner
%
% History
% 26.5.02 created
%
global source
% compute plane wave source
source=(E0 * exp(I*kx*g_eometry(3,:) + I*ky*g_eometry(4,:))).';