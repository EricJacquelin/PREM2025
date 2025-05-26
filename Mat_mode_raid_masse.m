% 16th February 2024
% 
% Matrices related to uncertain matrice and modes:
% 
% * PKP(:,:,i+1)=Psi_det'*K_i*Psi_det
% * PMP(:,:,i+1)=Psi_det'*M_i*Psi_det
% 
% where Psi_det is the deterministic modal matrix 
%
% Vectors extracted from uncertain matrices
%
%  * B(l,i+1,k)=-PKP(l,k,i+1)   % B(l,k).B(i+1)=-PhiKPhi(i+1).K(l,k)
%
% Inputs: 
%		Psi_det, M_det, K_det
%		M_u, K_u
%		M3_sp
%		M4_sp
% Outputs:
%		B, MM3, MM4

% uncertain mass and stiffness modal matrices PMP and PKP

PMP(:,:,1)=Psi_det'*M_det*Psi_det;      % i=1 => mean mass matrix
PKP(:,:,1)=Psi_det'*K_det*Psi_det;      % i=1 => mean stiffness matrix
for i=1:n_uncer  
    PMP(:,:,i+1)=Psi_det'*M_u(i).M*Psi_det;
    PKP(:,:,i+1)=Psi_det'*K_u(i).K*Psi_det;
end

% Vector B related to PKP

for k=1:n_ddl
    for l=1:n_ddl
        B_int=zeros(PC_nb,1);
        B_int(1)=-PKP(l,k,1);
        for i=1:n_uncer  
            B_int(i+1)=-PKP(l,k,i+1)/Coef1_PC;
        end
        B((l-1)*PC_nb+1:l*PC_nb,k)=B_int;
    end
end

% MM3 and KM3: uncertain mass and stiffness modal matrix connected to 3th-order moment matrix

clear KM3_int MM3_int KM3 MM3

for l=1:n_ddl
    for n=1:n_ddl  % for n=1:n_ddl-1
        KM3_int=zeros(PC_nb);  % KM3_int=KM3^{ln}
        MM3_int=zeros(PC_nb);  % MM3_int=MM3^{ln}       
        for i=1:n_uncer+1   
            KM3_int=KM3_int+ PKP(l,n,i) * full(M3_sp(i).mom3);
            MM3_int=MM3_int+ PMP(l,n,i) * full(M3_sp(i).mom3);
        end
        KM3((l-1)*PC_nb+1:l*PC_nb,(n-1)*PC_nb+1:n*PC_nb)=(KM3_int);
        MM3((l-1)*PC_nb+1:l*PC_nb,(n-1)*PC_nb+1:n*PC_nb)=(MM3_int);
    end
end

% MM4: uncertain mass modal matrix connected to 4th-order moment matrix

clear  MM4_int  MM4

for l=1:n_ddl
    for m=1:PC_nb
        for n=1:n_ddl%-1
            MM4_int=zeros(PC_nb);
            for i=1:n_uncer+1
                MM4_int=MM4_int+ PMP(l,n,i) * full(M4_sp(i,m).mom4);
            end
            MM4_sp(l,m).MM4(:,(n-1)*PC_nb+1:n*PC_nb)=sparse(MM4_int);
        end
    end
end



