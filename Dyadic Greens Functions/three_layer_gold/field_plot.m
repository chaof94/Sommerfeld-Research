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
trisurf(tri, x, z, real(y_ey));
% optional, could help make the plot look nicer
shading interp
colormap(brewermap([],'RdYlBu')) %RdYlGn 'RdYlBu' Spectral
colorbar
axis equal
% axis xy
hold on
line([-5,5],[0,0],[1,1],'linewidth',1.2,'color','black')
line([-5,5],[-1,-1],[1,1],'linewidth',1.2,'color','black')
view([0 90])
% Use Brewermap color schemes
% colormap(brewermap([],'Spectral')) %RdYlGn 'RdYlBu' Spectral
%figure
% figure;
% subplot(3,2,1)
% imagesc(xg2, zg2, real(y_ex_i));
% hold on
% line([-5,5],[0,0])
% line([-5,5],[-1,-1])
% title('(a) Re(E_x)');
% grid
% axis equal
% colormap(brewermap([],'Spectral')) %RdYlGn 'RdYlBu' Spectral
% colorbar
% axis tight
% axis xy
% 
% 
% 
% subplot(3,2,3)
% imagesc(xg2,zg2,real(y_ey_i));
% hold on
% line([-5,5],[0,0])
% line([-5,5],[-1,-1])
% title('(b) Re(E_y)');
% grid
% axis equal
% colorbar
% axis tight
% axis xy
% 
% subplot(3,2,5)
% imagesc(xg2,zg2,real(y_ez_i));
% hold on
% line([-5,5],[0,0])
% line([-5,5],[-1,-1])
% title('(c) Re(E_z)');
% grid
% axis equal
% colormap(brewermap([],'Spectral')) %RdYlGn 'RdYlBu' Spectral
% colorbar
% axis tight
% axis xy
% 
% 
% subplot(3,2,2)
% imagesc(xg2, zg2, imag(y_ex_i));
% hold on
% line([-5,5],[0,0])
% line([-5,5],[-1,-1])
% title('(d) Im(E_x)');
% grid
% axis equal
% colormap(brewermap([],'RdYlGn')) %RdYlGn 'RdYlBu' Spectral
% colorbar
% axis tight
% axis xy
% 
% subplot(3,2,4)
% imagesc(xg2,zg2,imag(y_ey_i));
% hold on
% line([-5,5],[0,0])
% line([-5,5],[-1,-1])
% title('(e) Im(E_y)');
% grid
% axis equal
% colorbar
% axis tight
% axis xy
% 
% subplot(3,2,6)
% imagesc(xg2,zg2,imag(y_ez_i));
% hold on
% line([-5,5],[0,0])
% line([-5,5],[-1,-1])
% title('(f) Im(E_z)');
% grid
% axis equal
% colormap(brewermap([],'RdYlBu')) %RdYlGn 'RdYlBu' Spectral
% colorbar
% axis tight
% axis xy






