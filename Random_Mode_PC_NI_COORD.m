% July 20th 2023
%
% Random modes derived with PC with a Non-intrusive method
% The PC are normalized <Phi_i^2>=1
%  
% n_uncer: number of uncertain parameters
%
% called by Random_Mode_MAIN.m
%

%=========
% Choice of the basis and calculation of the eigenvector elements on the selected basis

%______________________________________
% canonical basis
if choice_basis==1
    
    Psi_BASIS=eye(n_ddl);  % canonical basis
    
    % set of the canonical coefficients: ie projection of each sample of the i_mod-th
    % mode, on the i_mod-th canonical basis vector: as a consequence, we
    % get the initial coordinate
    n_a=n_ddl;   
    for i_mod=1:n_ddl  % number of the mode
        BASIS(i_mod).Psi_BASIS=Psi_BASIS;
        for i_a=1:n_a  % number of the coordinate
            BASIS(i_mod).a(:,i_a)=Psi_u_PC(i_mod).Psi_u(:,I_id)'*BASIS(i_mod).Psi_BASIS(:,i_a);    
        end
    end
    
%______________________________________    
% MODAL modes
elseif choice_basis==2   
    n_a=n_ddl;   % no truncature is done
    for i_mod=1:n_ddl
        BASIS(i_mod).Psi_BASIS=Psi_det;
        for i_a=1:n_a
          BASIS(i_mod).a(:,i_a)=Psi_u_PC(i_mod).Psi_u(:,I_id)'*M_det*BASIS(i_mod).Psi_BASIS(:,i_a)/M_mod_det(i_a);   
        end       
    end
    
%______________________________________
% POD modes
elseif choice_basis==3
    for i_mod=1:n_ddl
        % POD basis is derived from the identification sample results     
        [BASIS(i_mod).a,BASIS(i_mod).Psi_BASIS,BASIS(i_mod).lambda,err,sig_rec,energie_tot]=POD_temp_svd(Psi_u_PC(i_mod).Psi_u(:,I_id));
    end
    
end

% ===============================
% Dimensionless Eigenfrequency PCE

% PCE identification

disp('Eigenfrequencies PCE identification - Omega_n')
 
for i_mode=1:n_ddl
     i_mode;
     
    %===========================================================
    % samples to identify the PC model
   
    yy=Om_u_PC_ad(i_mode,I_id)';
    
    %===========================================================
    % samples to validate the PC model
    
    yy_val=Om_u_PC_ad(i_mode,I_val)';

    
    %===========================================================

    % sparse coefficients calculation with SVB and ARD
    clear ak_sparse Ind_sparse Ak_sparse Bk_sparse Vk_sparse L
    [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,yy,para_SVB_ini,e_ARD);
    clear ak_sparse
    ak_sparse=Phi(:,Ind_sparse)\yy;

    
    PCE_Om(i_mode).ak=ak_sparse;
    PCE_Om(i_mode).Ind_sparse=Ind_sparse;
    PCE_Om(i_mode).deg_PC_sparse=deg_PC(Ind_sparse,:);

    y_id=Phi(:,PCE_Om(i_mode).Ind_sparse)*PCE_Om(i_mode).ak;
    Om_PCE_id(i_mode,:)=y_id';
    
    ecart_id_om(i_mode,1)=norm(yy-y_id)/norm(yy)*100;

    
    % response obtained for the validtion parameters


    clear dep_PC_sparse Phi_val y_val_PCE

    Phi_val=Phi_def_norm(xx_val,don_sample,PCE_Om(i_mode).deg_PC_sparse);
    y_val_PCE=Phi_val*PCE_Om(i_mode).ak;

    
    Om_PCE_val(i_mode,:)=y_val_PCE';
    
    ecart_val_om(i_mode,1)=norm(yy_val-y_val_PCE)/norm(yy_val)*100;
 
end
% toc

disp('Squared Eigenfrequencies PCE identification: Omega^2_n')
% tic
for i_mode=1:n_ddl
     i_mode;
     
     clear yy yy_val
    %===========================================================
    % samples to identify the PC model
   
    yy=Om2_u_PC_ad(i_mode,I_id)';
    
    %===========================================================
    % samples to validate the PC model
    
    yy_val=Om2_u_PC_ad(i_mode,I_val)';

    %===========================================================

    % sparse coefficients calculation with SVB and ARD
    clear ak_sparse Ind_sparse Ak_sparse Bk_sparse Vk_sparse L
    [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,yy,para_SVB_ini,e_ARD);
    clear ak_sparse
    ak_sparse=Phi(:,Ind_sparse)\yy;
    
    PCE_Om2(i_mode).ak=ak_sparse;
    PCE_Om2(i_mode).Ind_sparse=Ind_sparse;
    PCE_Om2(i_mode).deg_PC_sparse=deg_PC(Ind_sparse,:);
    
    % response obtained for the parameters whih help to identify  the PCE:
    % results should be very good!
    y_id=Phi(:,PCE_Om2(i_mode).Ind_sparse)*PCE_Om2(i_mode).ak;
    Om2_PCE_id(i_mode,:)=y_id';
    
    ecart_id_om2(i_mode,1)=norm(yy-y_id)/norm(yy)*100;
    
    % response obtained for the validtion parameters

    clear dep_PC_sparse Phi_val y_val_PCE
    Phi_val=Phi_def_norm(xx_val,don_sample,PCE_Om2(i_mode).deg_PC_sparse);
    y_val_PCE=Phi_val*PCE_Om2(i_mode).ak;
    
    Om2_PCE_val(i_mode,:)=y_val_PCE';
    
    ecart_val_om2(i_mode,1)=norm(yy_val-y_val_PCE)/norm(yy_val)*100;
       
end

% ===============================
% modal mass PCE

% PCE identification

disp('modal mass PCE identification')
 
for i_mode=1:n_ddl
     i_mode;
     
    %===========================================================
    % samples to identify the PC model
   
    clear yy yy_val
    yy=M_mod_u_PC(i_mode,I_id)';
    
    %===========================================================
    % samples to validate the PC model
    
    yy_val=M_mod_u_PC(i_mode,I_val)';
    
    %===========================================================

    % sparse coefficients calculation with SVB and ARD
    clear ak_sparse Ind_sparse Ak_sparse Bk_sparse Vk_sparse L
    [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,yy,para_SVB_ini,e_ARD);
    clear ak_sparse
    ak_sparse=Phi(:,Ind_sparse)\yy;    
    
    PCE_M_mod(i_mode).ak=ak_sparse;
    PCE_M_mod(i_mode).Ind_sparse=Ind_sparse;
    PCE_M_mod(i_mode).deg_PC_sparse=deg_PC(Ind_sparse,:);
    
    % response obtained for the parameters whih help to identify  the PCE:
    % results should be very good!
    y_id=Phi(:,PCE_M_mod(i_mode).Ind_sparse)*PCE_M_mod(i_mode).ak;
    M_mod_PCE_id(i_mode,:)=y_id';
    
    ecart_id_M_mod(i_mode,1)=norm(yy-y_id)/norm(yy)*100;
    
    % response obtained for the validtion parameters

    clear dep_PC_sparse Phi_val y_val_PCE
    Phi_val=Phi_def_norm(xx_val,don_sample,PCE_M_mod(i_mode).deg_PC_sparse);
    y_val_PCE=Phi_val*PCE_M_mod(i_mode).ak;
    
    M_mod_PCE_val(i_mode,:)=y_val_PCE';
    
    ecart_val_M_mod(i_mode,1)=norm(yy_val-y_val_PCE)/norm(yy_val)*100;

end

% ===============================
% Eigenvector PCE

disp('Eigenvector PCE identification')
clear  PCE.Psi
clear yy yy_val y_id
 
n_a=n_ddl;
for i_mode=1:n_ddl
    
    Psi_u_PC(i_mode).Psi_u_id=zeros(n_ddl,npt_id);
    Psi_u_PC(i_mode).Psi_u_val=zeros(n_ddl,npt_val);
    
    %===========================================================
    % samples to identify/validate the PC model
    for i_Psi_element=1:n_a
        
             
        yy= BASIS(i_mode).a(:,i_Psi_element);     
        
        %===========================================================
        % PCE identification
        
        % sparse coefficients calculation with SVB and ARD
        clear ak_sparse Ind_sparse Ak_sparse Bk_sparse Vk_sparse L
        [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,BASIS(i_mode).a(I_id,i_Psi_element),para_SVB_ini,e_ARD);
        clear ak_sparse
        ak_sparse=Phi(:,Ind_sparse)\yy;
        
        PCE_Psi(i_mode).ak(i_Psi_element).ak=ak_sparse;
        PCE_Psi(i_mode).Ind_sparse(i_Psi_element).Ind=Ind_sparse;
        PCE_Psi(i_mode).deg_PC_sparse(i_Psi_element).deg=deg_PC(Ind_sparse,:);
        
        % response obtained for the parameters whih help to identify  the PCE:
        % results should be very good!
        y_id=Phi(:,PCE_Psi(i_mode).Ind_sparse(i_Psi_element).Ind)* PCE_Psi(i_mode).ak(i_Psi_element).ak;
        
        Psi_u_PC(i_mode).Psi_u_id=Psi_u_PC(i_mode).Psi_u_id+kron(y_id',BASIS(i_mode).Psi_BASIS(:,i_Psi_element));
               
        % response obtained for the validtion parameters
        
        clear dep_PC_sparse Phi_val y_val_PCE
        deg_PC_sparse=PCE_Psi(i_mode).deg_PC_sparse(i_Psi_element).deg;
        Phi_val=Phi_def_norm(xx_val,don_sample,deg_PC_sparse);
        y_val_PCE=Phi_val*PCE_Psi(i_mode).ak(i_Psi_element).ak;
                         
        Psi_u_PC(i_mode).Psi_u_val=Psi_u_PC(i_mode).Psi_u_val+kron(y_val_PCE',BASIS(i_mode).Psi_BASIS(:,i_Psi_element));
        
    end
    
    Psi_u_PC(i_mode).Psi_u_direct_id_Vec=reshape(Psi_u_PC(i_mode).Psi_u(:,I_id),npt_id*n_ddl,1);    
    Psi_u_PC(i_mode).Psi_u_PCE_id=reshape(Psi_u_PC(i_mode).Psi_u_id,npt_id*n_ddl,1);
    ecart_id_Phi(i_mode,1)=norm(Psi_u_PC(i_mode).Psi_u_direct_id_Vec-Psi_u_PC(i_mode).Psi_u_PCE_id)/norm(Psi_u_PC(i_mode).Psi_u_direct_id_Vec)*100;
                              
    Psi_u_PC(i_mode).Psi_u_val_Vec=reshape(Psi_u_PC(i_mode).Psi_u(:,I_val),npt_val*n_ddl,1);    
    Psi_u_PC(i_mode).Psi_u_PCE_val=reshape(Psi_u_PC(i_mode).Psi_u_val,npt_val*n_ddl,1);
    ecart_val_Phi(i_mode,1)=norm(Psi_u_PC(i_mode).Psi_u_val_Vec-Psi_u_PC(i_mode).Psi_u_PCE_val)/norm(Psi_u_PC(i_mode).Psi_u_val_Vec)*100;
    

end
