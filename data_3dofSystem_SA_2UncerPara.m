% April 10th 2024
%
% Simulation of an uncertain N=3-dof system (publi by S. Addikhari MSSP 164 (2022) 108260-section_6)
% However, in that case, there are 2 uncertain parameters: 
%           - the stiffness matrix is multiply by a first uncertain parameter 
%           - the mass matrix is multiply by a second uncertain parameter 
%
% Each random parameter follows a random distribution: the pdf is either
% uniform or normal 
%
% data
%

% sizes of the uncertain system

n_ddl=3;

n_uncer_k=1;
n_uncer_m=1;

n_uncer=n_uncer_k+n_uncer_m;

% uncertainties

epsi_m=0.15*ones(n_uncer_m,1);  
epsi_k=0.15*ones(n_uncer_k,1);  

% mechanical parameters

m_mean=1*ones(n_uncer_m,1); 

k_mean=1*ones(n_uncer_k,1); 

xi_det=0.01;  % constant modal damping ratio: not mentioned in the paper


%  elementary matrices

k1=1;k2=1;k3=1;k4=1;k5=1;k6=3.5;
EM(1).I_k=zeros(n_ddl);
EM(1).I_k(1,1)=k1+k4+k6;
EM(1).I_k(2,2)=k2+k4+k5;
EM(1).I_k(3,3)=k5+k3+k6;
EM(1).I_k(1,2)=-k4;EM(1).I_k(2,1)=-k4;
EM(1).I_k(2,3)=-k5;EM(1).I_k(3,2)=-k5;
EM(1).I_k(3,1)=-k6;EM(1).I_k(1,3)=-k6;

m1=1;m2=1;m3=1;
EM(1).I_m=eye(n_ddl);

%====================================
% Deterministic Matrices of the N-dof system

M_det=zeros(n_ddl);
for i_EM=1:n_uncer_m
    M_det_EM(i_EM).M=m_mean(i_EM)*EM(i_EM).I_m;
    M_det=M_det+M_det_EM(i_EM).M;
    M_u(i_EM).M=epsi_m(i_EM)*M_det_EM(i_EM).M;
    K_u(i_EM).K=zeros(n_ddl);
end

K_det=zeros(n_ddl);
for i_EM=1:n_uncer_k
    K_det_EM(n_uncer_m+i_EM).K=k_mean(i_EM)*EM(i_EM).I_k;
    K_det=K_det+K_det_EM(n_uncer_m+i_EM).K;    
    K_u(n_uncer_m+i_EM).K=epsi_k(i_EM)*K_det_EM(n_uncer_m+i_EM).K ;   
    M_u(n_uncer_m+i_EM).M=zeros(n_ddl);
end

[Vec_det,Val_det]=eig(K_det,M_det);
Om2_det=diag(Val_det);
Om_det=sqrt(Om2_det);
fr_det=Om_det/2/pi;
Psi_det=zeros(n_ddl);

if norme_Psi==1
    Psi_det=Vec_det;
    M_mod_det=diag(eye(n_ddl));
elseif norme_Psi==2
    for i_mod=1:n_ddl
        Psi_det(:,i_mod)=Vec_det(:,i_mod)/Vec_det(1,i_mod);
    end
    M_mod_det=diag(Psi_det'*M_det*Psi_det);
end


max_w=max(Om_det)*1.5;
n_fr=100;
dw=max_w/n_fr;
w=0:dw:max_w;
fr=w/2/pi;
