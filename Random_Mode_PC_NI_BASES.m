% October 17th, 2023
%
% Random modes derived with PC with a Non-intrusive method
% The PC are normalized <Phi_i^2>=1
%
% The PCE is also expanded on a specific basis.
% In this program, the possible bases are  
%   - the canonical basis [(1 0 0)T  (0 1 0)T  (0 0 1)T}
%   - the deterministic modal basis 
%   - the POD basis extracted from npt_id simulations
% 
% The advantage is a simplified expression thanks to a kronecker product
% The drawback is that the identification time is longer as the 
% coefficient vector that must be identified is bigger.
% 
% n_uncer: number of uncertain parameters
%
% called by Random_Mode_MAIN.m
%

%=========
% Choice of the basis
%______________________________________
% canonical basis
if choice_basis==1
    
    Psi_BASIS=eye(n_ddl);
    CAN_BASIS=Psi_BASIS;
    % set of the canonical coefficients: ie projection of each sample of the i_mode-th
    % mode, on the i_mod-th canonical basis vector.
    n_a=n_ddl;   % no truncature is done
    for i_mod=1:n_ddl
        for i_a=1:n_a
            BASIS(i_mod).a(:,i_a)=Psi_u_PC(i_mod).Psi_u(:,I_id)'*Psi_BASIS(:,i_a);    
        end
        BASIS(i_mod).Psi_BASIS=Psi_BASIS;
    end
%______________________________________    
% MODAL modes
elseif choice_basis==2
    % deterministic modal basis, ie solution of the eigenproblem
    % eig(K_det,Mdet)
    n_a=n_ddl;   % no truncature is done
    for i_mod=1:n_ddl
        BASIS(i_mod).Psi_BASIS=Psi_det;
        for i_a=1:n_a
            BASIS(i_mod).a(:,i_a)=Psi_u_PC(i_mod).Psi_u(:,I_id)'*BASIS(i_mod).Psi_BASIS(:,i_a);   
        end       
    end
%______________________________________    
% POD modes
elseif choice_basis==3
    for i_mod=1:n_ddl
        [BASIS(i_mod).a,BASIS(i_mod).Psi_BASIS,BASIS(i_mod).lambda,err,sig_rec,energie_tot]=POD_temp_svd(Psi_u_PC(i_mod).Psi_u(:,I_id));
    end
    
end

% 
% ===============================
% Dimensionless Eigenfrequency PCE

% PCE identification

disp('Dimensionless Eigenfrequencies PCE identification')
% tic
for numero_om_PCE=1:n_ddl
    numero_om_PCE;
    %===========================================================
    % samples to identify the PC model
   
    yy=Om_u_PC_ad(numero_om_PCE,I_id)';
    
    %===========================================================
    % samples to validate the PC model
    
    yy_val=Om_u_PC_ad(numero_om_PCE,I_val)';

    
    %===========================================================

    % sparse coefficients calculation with SVB and ARD
    clear ak_sparse Ind_sparse Ak_sparse Bk_sparse Vk_sparse L

    [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,yy,para_SVB_ini,e_ARD);
    
    PCE_Om(numero_om_PCE).ak=ak_sparse;
    PCE_Om(numero_om_PCE).Ind_sparse=Ind_sparse;
    PCE_Om(numero_om_PCE).deg_PC_sparse=deg_PC(Ind_sparse,:);
    
    % response obtained for the parameters whih help to identify  the PCE:
    % results should be very good!
    y_id=Phi(:,Ind_sparse)*ak_sparse;
    Om_PCE_id(numero_om_PCE,:)=y_id;
    
    ecart_id_om(numero_om_PCE,1)=norm(yy-y_id)/norm(yy)*100;
    
    % response obtained for the validtion parameters

    clear dep_PC_sparse Phi_val y_val_PCE
    deg_PC_sparse=PCE_Om(numero_om_PCE).deg_PC_sparse;
    Phi_val=Phi_def_norm(xx_val,don_sample,deg_PC_sparse);
    y_val_PCE=Phi_val*PCE_Om(numero_om_PCE).ak;
    
    Om_PCE_val(numero_om_PCE,:)=y_val_PCE';
    
    ecart_val_om(numero_om_PCE,1)=norm(yy_val-y_val_PCE)/norm(yy_val)*100;
    
end


disp('Dimensionless Squared Eigenfrequencies PCE identification')

for numero_om_PCE=1:n_ddl
     numero_om_PCE;
     
     clear yy yy_val
    %===========================================================
    % samples to identify the PC model
   
    yy=Om2_u_PC_ad(numero_om_PCE,I_id)';
    
    %===========================================================
    % samples to validate the PC model
    
    yy_val=Om2_u_PC_ad(numero_om_PCE,I_val)';

    %===========================================================

    % sparse coefficients calculation with SVB and ARD
    clear ak_sparse Ind_sparse Ak_sparse Bk_sparse Vk_sparse L
    [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,yy,para_SVB_ini,e_ARD);
    
    PCE_Om2(numero_om_PCE).ak=ak_sparse;
    PCE_Om2(numero_om_PCE).Ind_sparse=Ind_sparse;
    PCE_Om2(numero_om_PCE).deg_PC_sparse=deg_PC(Ind_sparse,:);
    
    % response obtained for the parameters whih help to identify  the PCE:
    % results should be very good!
    y_id=Phi(:,PCE_Om2(numero_om_PCE).Ind_sparse)*PCE_Om2(numero_om_PCE).ak;
    Om2_PCE_id(numero_om_PCE,:)=y_id';
    
    ecart_id_om2(numero_om_PCE,1)=norm(yy-y_id)/norm(yy)*100;
    
    % response obtained for the validtion parameters

    clear dep_PC_sparse Phi_val y_val_PCE
    Phi_val=Phi_def_norm(xx_val,don_sample,PCE_Om2(numero_om_PCE).deg_PC_sparse);
    y_val_PCE=Phi_val*PCE_Om2(numero_om_PCE).ak;
    
    Om2_PCE_val(numero_om_PCE,:)=y_val_PCE';
    
    ecart_val_om2(numero_om_PCE,1)=norm(yy_val-y_val_PCE)/norm(yy_val)*100;
       
end

% ===============================
% modal mass PCE

% PCE identification

disp('modal mass PCE identification')

for numero_om_PCE=1:n_ddl
    numero_om_PCE;
    %===========================================================
    % samples to identify the PC model
   
    clear yy yy_val
    yy=M_mod_u_PC(numero_om_PCE,I_id)';
    
    %===========================================================
    % samples to validate the PC model
    
    yy_val=M_mod_u_PC(numero_om_PCE,I_val)';
    
    %===========================================================

    % sparse coefficients calculation with SVB and ARD
    clear ak_sparse Ind_sparse Ak_sparse Bk_sparse Vk_sparse L
    [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,yy,para_SVB_ini,e_ARD);
    
    PCE_M_mod(numero_om_PCE).ak=ak_sparse;
    PCE_M_mod(numero_om_PCE).Ind_sparse=Ind_sparse;
    PCE_M_mod(numero_om_PCE).deg_PC_sparse=deg_PC(Ind_sparse,:);
    
    % response obtained for the parameters whih help to identify  the PCE:
    % results should be very good!
    y_id=Phi(:,PCE_M_mod(numero_om_PCE).Ind_sparse)*PCE_M_mod(numero_om_PCE).ak;
    M_mod_PCE_id(numero_om_PCE,:)=y_id';
    
    ecart_id_M_mod(numero_om_PCE,1)=norm(yy-y_id)/norm(yy)*100;
    
    % response obtained for the validtion parameters

    clear dep_PC_sparse Phi_val y_val_PCE
    Phi_val=Phi_def_norm(xx_val,don_sample,PCE_M_mod(numero_om_PCE).deg_PC_sparse);
    y_val_PCE=Phi_val*PCE_M_mod(numero_om_PCE).ak;
    
    M_mod_PCE_val(numero_om_PCE,:)=y_val_PCE';
    
    ecart_val_M_mod(numero_om_PCE,1)=norm(yy_val-y_val_PCE)/norm(yy_val)*100;

end
% toc


% ===============================
% Eigenvector PCE

disp('Eigenvector PCE identification')
clear PCE.Psi
clear y_val_PCE

for numero_VP_PCE=1:n_ddl
    numero_VP_PCE;
    %===========================================================
    % samples to identify/validate the PC model
    
    clear yy yy_val 
   
    yy=reshape(Psi_u_PC(numero_VP_PCE).Psi_u(:,I_id),n_ddl*npt_id,1);
    yy_val=reshape(Psi_u_PC(numero_VP_PCE).Psi_u(:,I_val),n_ddl*npt_val,1);
        
    %===========================================================
    % PCE identification
    
    % sparse coefficients calculation with SVB and ARD

    Phi_Psi=kron(Phi,BASIS(numero_VP_PCE).Psi_BASIS);
    [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi_Psi,yy,para_SVB_ini,e_ARD);

    PCE_Psi(numero_VP_PCE).ak_sparse=ak_sparse;
    PCE_Psi(numero_VP_PCE).Ind_sparse=Ind_sparse;
        
    for compo_Psi=1:n_ddl
        clear Ind_temp
        Ind_temp=intersect(find(Ind_sparse>(compo_Psi-1)*PC_nb),find(Ind_sparse<=compo_Psi*PC_nb));
        PCE_Psi(numero_VP_PCE).deg_PC_sparse(compo_Psi).deg=deg_PC(Ind_sparse(Ind_temp)-(compo_Psi-1)*PC_nb,:);
    end
    
        % response obtained for the parameters which help to identify  the PCE:
        % results should be very good!
    y_id=Phi_Psi(:,PCE_Psi(numero_VP_PCE).Ind_sparse)*PCE_Psi(numero_VP_PCE).ak_sparse;

    Psi_u_PC(numero_VP_PCE).Psi_u_direct_id_Vec=reshape(Psi_u_PC(numero_VP_PCE).Psi_u(:,I_id),npt_id*n_ddl,1);  
    Psi_u_PC(numero_VP_PCE).Psi_u_PCE_id=y_id;
    ecart_id_Phi(numero_VP_PCE,1)=norm(Psi_u_PC(numero_VP_PCE).Psi_u_direct_id_Vec-Psi_u_PC(numero_VP_PCE).Psi_u_PCE_id)/norm(Psi_u_PC(numero_VP_PCE).Psi_u_direct_id_Vec)*100;     
           
    % response obtained for the validation parameters
        
    Phi_val=Phi_def_norm(xx_val,don_sample,deg_PC);
    Phi_Psi_val=kron(Phi_val,BASIS(numero_VP_PCE).Psi_BASIS);
    y_val_PCE=Phi_Psi_val(:,PCE_Psi(numero_VP_PCE).Ind_sparse)*PCE_Psi(numero_VP_PCE).ak_sparse;
    Psi_u_PC(numero_VP_PCE).Psi_u_PCE_val=y_val_PCE;
         
    Psi_u_PC(numero_VP_PCE).Psi_u_val_Vec=reshape(Psi_u_PC(numero_VP_PCE).Psi_u(:,I_val),npt_val*n_ddl,1);    
    ecart_val_Phi(numero_VP_PCE,1)=norm(Psi_u_PC(numero_VP_PCE).Psi_u_val_Vec-Psi_u_PC(numero_VP_PCE).Psi_u_PCE_val)/norm(Psi_u_PC(numero_VP_PCE).Psi_u_val_Vec)*100;
     
end

