clear;close all

%Change res depending on the resolution of your data (it should be matched with C++ code in line 149)
res = 50;

data_ex=load('Ex.dat');
data_ey=load('Ey.dat');
data_ez=load('Ez.dat');
% x=-5:10/(res-1):5;
x = data_ex(:,1);
% z=-3:5.0/(res-1):2.0;
z = data_ex(:,2);
[xg,zg] = meshgrid(x,z);

ex_realpart = data_ex(:,3);
ex_imagpart = data_ex(:,4);
y_ex = ex_realpart+1i*ex_imagpart;
% y_ex = reshape(y_ex,[res,res]);

ey_realpart = data_ey(:,3);
ey_imagpart = data_ey(:,4);
y_ey = ey_realpart+1i*ey_imagpart;
% y_ey = reshape(y_ey,[res,res]);

ez_realpart = data_ez(:,3);
ez_imagpart = data_ez(:,4);
y_ez = ez_realpart+1i*ez_imagpart;
% y_ez = reshape(y_ez,[res,res]);

%For better figure interpolate data by res2 resolution
res2 = 400;
xg2 = -5:10/(res2-1):5;
zg2 = -3:5/(res2-1):2;
[xg3, zg3] = meshgrid(xg2, zg2);
% y_ex_i = interp2(xg, zg, y_ex, xg3, zg3);
% y_ey_i = interp2(xg, zg, y_ey, xg3, zg3);
% y_ez_i = interp2(xg, zg, y_ez, xg3, zg3);

% triangulate and plot
tri = delaunay(x, z);
%% E_x
figure(1)
box on
set(gcf,'color','white');
set(groot,'defaulttextinterpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
subplot(211)
trisurf(tri, x, z, real(y_ex));
title('$\Re (E_x)$','interpreter','latex');
% optional, could help make the plot look nicer
shading interp
colormap(brewermap([],'RdYlBu')) %RdYlBu 'RdYlBu' Spectral
% colormap('viridis')
% colormap(cubehelix)
colorbar
axis equal
% axis xy
hold on

line([-5,5],[0,0],[1,1],'linewidth',1.2,'color','black')
line([-5,5],[-1,-1],[1,1],'linewidth',1.2,'color','black')
view([0 90])



subplot(212)
trisurf(tri, x, z, imag(y_ex));
title('$\Im (E_x)$','interpreter','latex');
% optional, could help make the plot look nicer
shading interp
colormap(brewermap([],'RdYlBu')) %RdYlBu 'RdYlBu' Spectral
% colormap('viridis')
colorbar
axis equal
% axis xy
hold on

line([-5,5],[0,0],[1,1],'linewidth',1.2,'color','black')
line([-5,5],[-1,-1],[1,1],'linewidth',1.2,'color','black')
view([0 90])


%% Ey
figure(2)
box on
set(gcf,'color','white');
set(groot,'defaulttextinterpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
subplot(211)
trisurf(tri, x, z, real(y_ey));
title('$\Re (E_y)$','interpreter','latex');
% optional, could help make the plot look nicer
shading interp
colormap(brewermap([],'RdYlBu')) %RdYlBu 'RdYlBu' RdYlBu
% colormap('viridis')
colorbar
axis equal
% axis xy
hold on

line([-5,5],[0,0],[1,1],'linewidth',1.2,'color','black')
line([-5,5],[-1,-1],[1,1],'linewidth',1.2,'color','black')
view([0 90])



subplot(212)
trisurf(tri, x, z, imag(y_ey));
title('$\Im (E_y)$','interpreter','latex');
% optional, could help make the plot look nicer
shading interp
colormap(brewermap([],'RdYlBu')) %RdYlBu 'RdYlBu' RdYlBu
% colormap('viridis')
colorbar
axis equal
% axis xy
hold on

line([-5,5],[0,0],[1,1],'linewidth',1.2,'color','black')
line([-5,5],[-1,-1],[1,1],'linewidth',1.2,'color','black')
view([0 90])



%% Ez
figure(3)
box on
set(gcf,'color','white');
set(groot,'defaulttextinterpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
subplot(211)
trisurf(tri, x, z, real(y_ez));
title('$\Re (E_z)$','interpreter','latex');
% optional, could help make the plot look nicer
shading interp
colormap(brewermap([],'RdYlBu')) %RdYlBu 'RdYlBu' RdYlBu
% colormap('viridis')
axis equal
% axis xy
hold on

line([-5,5],[0,0],[1,1],'linewidth',1.2,'color','black')
line([-5,5],[-1,-1],[1,1],'linewidth',1.2,'color','black')
view([0 90])
colorbar

subplot(212)
trisurf(tri, x, z, imag(y_ez));
title('$\Im (E_z)$','interpreter','latex');
% optional, could help make the plot look nicer
shading interp
% colormap(brewermap([],'RdYlBu')) %RdYlBu 'RdYlBu' RdYlBu
% colormap('viridis')
axis equal
% axis xy
hold on

line([-5,5],[0,0],[1,1],'linewidth',1.2,'color','black')
line([-5,5],[-1,-1],[1,1],'linewidth',1.2,'color','black')
view([0 90])
colorbar
% myaa('publish');
cleanfigure();
matlab2tikz('filename',sprintf('Ez_z_dip.tex'),'showInfo', false)


%% This is the contour plot
% figure(4)
% box on
% set(gcf,'color','white');
% set(groot,'defaulttextinterpreter','latex');
% set(gca,'TickLabelInterpreter', 'latex');
% set(gca,...
%     'box','on',...
%     'FontName','times new roman',...
%     'FontSize',12);
% % BS = textread([char(fname(1)) '.txt']);
% X = data_ez(:,1);
% Y = data_ez(:,2);
% Xrange = linspace(min(X),max(X),length(X));
% Yrange = linspace(min(Y),max(Y),length(Y));
% 
% [XI,YI] = meshgrid(Xrange,Yrange);
% 
% z = imag(y_ex);
% Z = griddata(X,Y,z,XI,YI,'v4');
% g = contourf(XI,YI,Z,'linecolor','black');
% axis equal;
% shading interp;
% colormap(brewermap([],'RdYlBu')) %RdYlBu 'RdYlBu' RdYlBu
% colorbar;
% set(0,'defaulttextinterpreter','latex')
% set(gca, 'FontName', 'Times New Roman','FontSize',12);
% hold on
% line([-5,5],[0,0],[1,1],'linewidth',1.2,'color','black')
% line([-5,5],[-1,-1],[1,1],'linewidth',1.2,'color','black')
% xlabel('x');


