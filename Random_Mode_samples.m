% 4th August 2023
%
% Simulation of an uncertain N-dof system; each rigidity, each mass is random: 
% Each random parameter follows a statistical distribution: the pdf is either
% uniform or normal 
%
% To identify the PCE through a NI method, samples are required.
% They are calculted in this subroutine

%

%===========================================================
% random mode samples

for i_u=1:npt_tot
       
    K_u_PC=zeros(n_ddl);
    M_u_PC=zeros(n_ddl);
    for i_EM=1:n_uncer
        K_u_PC =K_u_PC +z_u_PC(i_u,i_EM)*K_u(i_EM).K;		
    	M_u_PC=M_u_PC+z_u_PC(i_u,i_EM)*M_u(i_EM).M;
    end
    
    M_tot_PC=M_det+M_u_PC;
    K_tot_PC=K_det+K_u_PC;    
 
    M_tot_Vec(:,i_u)=diag(M_tot_PC);
    
    [Vec_u,Val_u]=eig(K_tot_PC,M_tot_PC);
    Om2_u_PC(:,i_u)=diag(Val_u);
    Om_u_PC(:,i_u)=sqrt(Om2_u_PC(:,i_u));

   
    if norme_Psi==1
        for i_mod=1:n_ddl
            Psi_u_PC(i_mod).Psi_u(:,i_u)=Vec_u(:,i_mod);
            M_mod_u_PC(i_mod,i_u)=Psi_u_PC(i_mod).Psi_u(:,i_u)'*M_tot_PC*Psi_u_PC(i_mod).Psi_u(:,i_u);
        end
    elseif norme_Psi==2
        for i_mod=1:n_ddl
            Psi_u_PC(i_mod).Psi_u(:,i_u)=Vec_u(:,i_mod)/Vec_u(1,i_mod);
            M_mod_u_PC(i_mod,i_u)=Psi_u_PC(i_mod).Psi_u(:,i_u)'*M_tot_PC*Psi_u_PC(i_mod).Psi_u(:,i_u);
        end
    end
    
        
end

for i_mod=1:n_ddl
    Om2_u_PC_ad(i_mod,:)=Om2_u_PC(i_mod,:)/Om2_det(i_mod);
    Om_u_PC_ad(i_mod,:)=Om_u_PC(i_mod,:)/Om_det(i_mod);
end
