%% Calculate 2DEG condunctivity
clear all;close all
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
N1 = 1e17 *1e6 ; %cm^-3
eps_0 = 8.854e-12;
eps_s = 12*eps_0;
eps_ox = 4*eps_0;
a = 3.52e-6;
d = .14e-6;
q = 2*pi/a;
m_e = 9.31e-31;
m_star = .067*m_e;
e = 1.602e-19;
e2 = (2.3068e-28);
N = (0:.01:4)*N1;
omega_p = sqrt(N*e2/m_star * q/(eps_ox*(eps_ox + eps_s*coth(q*d))));
close all
yyaxis left
plot(N,omega_p)
ylabel('\omega_p [rad/s]')
yyaxis right
k = omega_p/c;
plot(N,k)
ylabel('k [/cm]')

figure(2)
plot(k,omega_p)

figure(3)
w = (2*pi)*linspace(1e11,1e14,1e3);

% Data from LLoyd-Hughes Review of THz conductivity
m_e = 9.31e-31;
m_star = .067*m_e;
Ns = 1e17*1e6;
tau = 100e-15;%1.3889e-13
e2 = (2.3068e-28);
q = 1.60218e-19; 
f = 1/3;
ep0_0 = 12.9;
ep0_inf = 11;
w_p = sqrt(N1*q^2/(m_star*ep0_0*eps_0));
w_0 = sqrt(f)*w_p;
% w_0 = 1e12;

sigma = N1*q^2*tau/m_star./(1 - 1i*w*tau); % Drude
sigma = N1*q^2*tau/m_star./(1 - (1i*tau*(w - w_0^2./w)));

close all;clf;figure(1)
N = 2;
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% *7.75e-5.
% [ax, Plot1, Plot2] = plotyy(w/(2*pi)*1e-12, (real(sigma)), w/(2*pi)*1e-12, (imag(sigma)));
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')


% hold on
% ylabel('$\mathrm{Real}$','interpreter','latex','FontSize',11)
% yyaxis right
% plot(w/(2*pi)*1e-12,(real(sigma)),w/(2*pi)*1e-12,(imag(sigma)));
% symlog()

plot(w/(2*pi)*1e-12,(real(sigma)), 'linewidth',1.4);
hold on
plot(w/(2*pi)*1e-12,(imag(sigma)), 'linewidth',1.4);
% ylabel('$\mathrm{Imaginary}$','interpreter','latex','FontSize',11)
set(gcf,'Color','white');
set(gca,'FontName','times new roman','FontSize',11,'YScale', 'lin','XScale', 'log') % Set axes fonts to Times New Roman
% % xlim([6e12 10e12])
% xlim([6e-1 1e1])
% xlim([0 100])
%
ylabel('$\mathrm{\sigma_s}(\mathrm{S/m})$','interpreter','latex','FontSize',11)
xlabel('$f (\mathrm{THz})$','interpreter','latex')
legend({'$\mathrm{\sigma_s,_{real}}$', '$\mathrm{\sigma_s,_{imag}}$'},...
    'location','southeast','interpreter','latex','FontSize',11);
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% grid on
box on    
matlab2tikz('filename',sprintf('conductivity_gaas.tex'));