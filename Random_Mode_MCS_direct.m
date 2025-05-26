% 22th August 2023
%
% Simulation of an uncertain N-dof system; each rigidity, each mass is random: 
% Each random parameter follows a statistical distribution: the pdf is either
% uniform or normal 
%
% Called by Random_Mode_Main_NI.m / Random_Mode_Main_INTRUSIVE.m

ij=sqrt(-1);

% modes aléatoires

t_ini_MCS=cputime;


for i_u=1:npt_MCS_direct
    
    
    Raid_u=zeros(n_ddl);
    Masse_u=zeros(n_ddl);
    for i_EM=1:n_uncer
        Raid_u =Raid_u +z_u_MCS(i_u,i_EM)*K_u(i_EM).K;		
    	Masse_u=Masse_u+z_u_MCS(i_u,i_EM)*M_u(i_EM).M;
    end
    
    M_tot=M_det+Masse_u;
    K_tot=K_det+Raid_u;
    
    [Vec_u,ValP_u]=eig(K_tot,M_tot);
    Om2_u(:,i_u)=diag(ValP_u);
    Om_u(:,i_u)=sqrt(Om2_u(:,i_u));
    
    if norme_Psi==1
        for i_mod=1:n_ddl
            Psi_u(i_mod).Psi_u(:,i_u)=Vec_u(:,i_mod);
            M_mod_u(i_mod,i_u)=diag(Psi_u(i_mod).Psi_u(:,i_u)'*M_tot*Psi_u(i_mod).Psi_u(:,i_u));
        end
    elseif norme_Psi==2
        for i_mod=1:n_ddl
            Psi_u(i_mod).Psi_u(:,i_u)=Vec_u(:,i_mod)/Vec_u(1,i_mod);
            M_mod_u(i_mod,i_u)=diag(Psi_u(i_mod).Psi_u(:,i_u)'*M_tot*Psi_u(i_mod).Psi_u(:,i_u));
        end
    end
       

end

T_MCS=cputime-t_ini_MCS;

% Résultats

Om2_moy_direct=mean(Om2_u,2);
Om_moy_direct=mean(Om_u,2);
Om_std_direct=std(Om_u,1,2);

compare_Om_det_moy=[Om_det Om_moy_direct];


for i_mod=1:n_ddl
        Psi_u(i_mod).Psi_moy=mean(Psi_u(i_mod).Psi_u,2);
end



