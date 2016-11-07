% Movie
%
% 2D Scattering
% by Beat Hangartner
%
function E2=mov(start,stop,step)
E2=[];
write_jpg=1;
parameters
g_eom
range=start:step:(stop-step);
tic
for lamda=range
    parameters
    excitation
    solver_mat2
    if write_jpg==1
        figure(2)
        contourf(S.*conj(S))
        title(['Scattering pattern lamda=',num2str(lamda),...
            ' Energy=',num2str(Energy)])
        ylabel(['length = ',num2str(len)])
        xlabel(['width = ',num2str(width)])
        axis image
        pause(.01)
%         colorbar
%         saveas(gcf,['./',filename,'/' ,int2str(lamda*100), '.eps']);
%         clf reset
    end
    E2=[E2,Energy];
end
toc
clf reset
figure(3)
plot(range,E2);
% title('Average energy inside scatterer ', filename)
xlabel('wavelength')
ylabel('Energy')
% saveas(gcf,['./',filename,'/' ,'Energy-',num2str(start),'-',num2str(stop),...
%     '-',num2str(step), '.eps']);