% Geometry array definition
%
% 2D Scattering
% by Beat Hangartner
%
% History
% 26.5.02 created
%
%
%%%%%%%%%%%%%%%%%%
% Geometry layout
%%%%%%%%%%%%%%%%%%%%%%
% Origin
% (0,0)
% |-----------------------------------> x-axis (width)
% | cell 1 | cell 2 |
% | x=0.5*dx | x=1.5*dx |
% | y=0.5*dy | y=0.5*dy |
% | | |
% | | |
% | | |
% |------------------------------------
% | cell 3 | cell 4 |
% | x=0.5*dx | x=1.5*dx |
% | y=1.5*dy | y=1.5*dy |
% | | |
% | | |
% | | |
% |-----------------------------------
% |
% \/
% y-axis (length)
global g_eometry
global bg_indices
global mat_indices
% global filename
%instantiate geometry cells in a linear array
%structure is as follows:
% 1. row index
% 2. row (complex) epsilon
% 3. row x-position
% 4. row y-position
g_eometry =zeros(4,geom_size);
g_eometry(1,:)=1:geom_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open material file
% and retrieve scatterer's layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fid = fopen(['./',filename,'/',filename,'.m']);
fid = fopen('020-010-block1.m');
file = fread(fid)';
fclose(fid);
%get rid of carriage return and line feed (PC only)
file=file(file~=10);
file=file(file~=13);
% impose vlaues coded in the file over the geometry
%('0' is ASCII character 48)
% only integers from 0 to 9 supported
% at the moment for epsilon -;)
g_eometry(2,:) = file-48;
% generate regular grid position information
for i=1:len/dy
    g_eometry(3,(i-1)*width/dx+1:i*width/dx)=0.5*dx:dx:width;
end
for i=1:width/dx
    g_eometry(4,[i-1+[1:width/dx:geom_size]])=0.5*dy:dy:len;
end
%extract background/material inidces; to be computed separatly
bg_indices = find(g_eometry(2,:)==0);
mat_indices = find(g_eometry(2,:)~=0);