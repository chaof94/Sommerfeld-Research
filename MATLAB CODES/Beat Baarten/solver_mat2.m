% Solver
%
% 2D Scattering
% by Beat Hangartner
%
% History
% 26.5.02 created
% 29.5.02 modified Algorithm corrections (epsilons)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute elements with contrast other than zero %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Energy
%number of material elements
mat_size = size(mat_indices,2);
bg_size = size(bg_indices,2);
rho=zeros(mat_size);
k=1;
for i=mat_indices
    rho(k,:)=sqrt((g_eometry(3,i) - g_eometry(3,mat_indices)).^2 +...
        (g_eometry(4,i) - g_eometry(4,mat_indices)).^2);
    k=k+1;
end
% put them into the hankelfunction
G_mat = I/4*V_i*k0^2*besselh(0,1,k0*rho);
%free memory
clear rho;
% add the correct diagonal elements
G_mat=full(spdiags((I*pi/2*beta*gamma)*k0^2*ones(mat_size,1),0,G_mat));
% apply epsilon
k=1;
for i=mat_indices
    G_mat(:,k) = G_mat(:,k)* g_eometry(2,i);
    k=k+1;
end
A_mat = eye(mat_size) - G_mat;
%free memory
clear G_mat;
solution_mat = (A_mat\source(mat_indices)).';
%free memory
clear A_mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute background solutions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute rho's between background and material
rho=zeros(bg_size,mat_size);
k=1;
for i=bg_indices
    rho(k,:)=sqrt((g_eometry(3,i) - g_eometry(3,mat_indices)).^2 + ...
        (g_eometry(4,i) - g_eometry(4,mat_indices)).^2);
    k=k+1;
end
% put them into the hankelfunction
G_bg = I/4* V_i*k0^2*besselh(0,1,k0*rho);
%free memory
clear rho;
solution_bg = source(bg_indices).' + ...
    (G_bg * ( solution_mat .* g_eometry(2, mat_indices)).' ) .';
%free memory
clear G_bg;
solution = zeros(1,geom_size);
solution(mat_indices)=solution_mat;
solution(bg_indices)=solution_bg;
% concatenate solution to matrix S
S=solution(1:width/dx);
for i=2:len/dy
    S = [solution((i-1)*width/dx+1:i*width/dx);S];
end
Energy = sum(solution(mat_indices).*conj(solution(mat_indices)))/mat_size;