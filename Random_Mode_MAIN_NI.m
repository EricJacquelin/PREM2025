% 22th August 2023
%
% Simulation of an uncertain N-dof system; each rigidity, each mass is random: 
% Each random parameter follows a statistical distribution: the pdf is either
% uniform or normal 
%
% The random mode are derived from two MCS procedure
%
% - random_Mode_MCS_direct.m: direct MCS 
% - random_mode_PC_NI_BASES.m or random_mode_PC_NI_COORD.m: MCS from NI-PC-direct
%
% For the PC_NI methods, 3 parameters must be fixed:
%   1) method: = 1 if each coordinate of the eigenvector is individually identified
%                => random_mode_PC_NI_COORD.m;                
%                it seems the better choice (at least, it is quicker!)              
%              = 2 if each eigenvector is globally identified 
%                => random_mode_PC_NI_BASES.m; 
%   2) choice_basis: = 1: canonical basis; =2: modal basis  ; =3: POD basis
%   3) norme_Psi:    = 1: normalization wrt the mass matrix => M_mod=eyes(n_ddl) ; 
%                    = 2: the 1st element of each eigenvector=1 
%
% The data are called in data_NdofSystem.m
% 

%

clear

% subdirectory

path(path,'.\PCE');
path(path,'.\SVB');
path(path,'.\POD');

disp('  ')
disp('  ')
disp('-----New simulation------')

ij=sqrt(-1);
imprime = 0; % to save my plots into files
courbes = 1     % = 0 no curve is plotted; = 1 curves are plotted;

%=====================================
% Random mode determination method, bassis, normalization

% direct PCE AND identification of each coefficient given in the canonical basis OR the deterministic modal basis OR the POD basis method = 1
% direct PCE AND global identification through the canonical basis OR deterministic modal OR POD basis: method = 2
method=1;

choice_basis=1; % = 1: canonical basis ; = 2: modal basis  ; = 3: POD basis

norme_Psi=2;    % = 1: normalization with respect to the mass matrix => M_mod=eyes(n_ddl) ; = 2: the 1st element of each eigenvector=1 

%=====================================
% dynamical system data file: uncomment/comment the studied dynamical system data file

% data_3dofSystem_SA_2UncerPara     % 3 dof dynamical system, 2 uncertain parameters
data_3dofSystem_SA_9UncerPara     % 3 dof dynamical system, 9 uncertain parameters

%=====================================
% statistics

% MCS
npt_MCS_direct=10000;  % number of points used for the final MCS validation tests

% pdf 
if n_uncer==2
    pdf_law='uniform'
elseif n_uncer==9
    pdf_law='normal'
end

% choice of the PC family
if strcmp(pdf_law,'uniform')==1
    poly='Legendre';
elseif strcmp(pdf_law,'normal')==1
    poly='Hermite';
end

%=====================================
% MCS samples

clear z_u_MCS 
z_u_MCS=LHS_matlab(npt_MCS_direct,n_uncer,pdf_law);  

%=====================================
% Random modes: DIRECT MCS 
%
% Inputs: z_u_MCS(i_u,i_EM)           with i_u=1:npt_MCS_direct and i_EM=1:n_uncer
%
% Outputs: Psi_u(i_mod).Psi_u(:,i_u) 
%          Om_u(:,i_u)                with i_u=1:npt_MCS_direct and i_mod=1:3

t_MCS_ini=cputime;
Random_Mode_MCS_direct
t_MCS=cputime-t_MCS_ini;

%=====================================
% Random modes: PCE metamodel

t_ini=cputime;

% data for NI-PCE included the identification/validation samples using Random_Mode_samples.m 

data_RM_PC_NI 

% Identification of the PCE coefficients

if method==1
    Random_Mode_PC_NI_COORD
elseif method==2
    Random_Mode_PC_NI_BASES
end

% Evaluation of the MCS validation data

% => Eigenfrequencies
Om_PCE=zeros(npt_MCS_direct,n_ddl);
for numero_om=1:n_ddl
    clear deg_PC_sparse Phi_simul 
    deg_PC_sparse=PCE_Om(numero_om).deg_PC_sparse;
    Phi_simul=Phi_def_norm(z_u_MCS,don_sample,deg_PC_sparse);
    Om_PCE(:,numero_om)=Om_det(numero_om)*Phi_simul*PCE_Om(numero_om).ak;
end
Om_PCE=Om_PCE';

% => Squared eigenfrequencies
Om2_PCE=zeros(npt_MCS_direct,n_ddl);
for numero_om=1:n_ddl
    clear deg_PC_sparse Phi_simul 

    degre_om_PC_sparse.nb(numero_om)=size(PCE_Om2(numero_om).deg_PC_sparse,1);
    degre_om_PC_sparse.max(numero_om)=max(sum(PCE_Om2(numero_om).deg_PC_sparse,2));
    
    deg_PC_sparse=PCE_Om2(numero_om).deg_PC_sparse;
    Phi_simul=Phi_def_norm(z_u_MCS,don_sample,deg_PC_sparse);
    Om2_PCE(:,numero_om)=Om2_det(numero_om)*Phi_simul*PCE_Om2(numero_om).ak;
end
Om2_PCE=Om2_PCE';

% => Eigenvectors
n_a=n_ddl;
for numero_VP=1:n_ddl
    
    if method==1
        clear yy_simul 
        Psi_u_PC(numero_VP).Psi_u=zeros(n_ddl,npt_MCS_direct);
        for i_a=1:n_a
            clear deg_PC_sparse Phi_simul
            degre_Psi_PC_sparse.nb(i_a,numero_VP)=size(PCE_Psi(numero_VP).deg_PC_sparse(i_a).deg,1);
            degre_Psi_PC_sparse.max(i_a,numero_VP)=max(sum(PCE_Psi(numero_VP).deg_PC_sparse(i_a).deg,2));
            
            deg_PC_sparse=PCE_Psi(numero_VP).deg_PC_sparse(i_a).deg;
            Phi_simul=Phi_def_norm(z_u_MCS,don_sample,deg_PC_sparse);
            yy_simul=Phi_simul*PCE_Psi(numero_VP).ak(i_a).ak;
            Psi_u_PC(numero_VP).Psi_u=Psi_u_PC(numero_VP).Psi_u+kron(yy_simul',BASIS(numero_VP).Psi_BASIS(:,i_a));
        end
     elseif method==2
        clear deg_PC_sparse Phi_simul  yy_simul
        Phi_simul=Phi_def_norm(z_u_MCS,don_sample,deg_PC);
        Phi_Psi_simul=kron(Phi_simul,BASIS(numero_VP).Psi_BASIS);
        yy_simul=Phi_Psi_simul(:,PCE_Psi(numero_VP).Ind_sparse)*PCE_Psi(numero_VP).ak_sparse;
        Psi_u_PC(numero_VP).Psi_u=reshape(yy_simul,n_ddl,npt_MCS_direct);
    end
end

% => Modal Masses

for i_u=1:npt_MCS_direct        
    Masse_u=zeros(n_ddl);
    for i_EM=1:n_uncer
    	Masse_u=Masse_u+z_u_MCS(i_u,i_EM)*M_u(i_EM).M;
    end    
    M_tot=M_det+Masse_u;
    for k=1:n_ddl
        M_mod_PCE(k,i_u)=Psi_u_PC(k).Psi_u(:,i_u)'*M_tot*Psi_u_PC(k).Psi_u(:,i_u);
    end
end

t_simul=cputime-t_ini;

%% ========================================
% Comparisons between direct MCS and PCE metamodel MCS

disp('----- Comparison ------')


% Om2: mean
for i=1:n_ddl
    Om2_moy_PCE_a(i,1)=Om2_det(i)*PCE_Om2(i).ak(1);  % 
end

Om2_moy_PCE=mean(Om2_PCE,2);
Om2_std_PCE=std(Om2_PCE,1,2);

if n_uncer==2   &&  strcmp(pdf_law,'uniform')==1
    Om2_moy_theo=Om2_det/2/epsi_m*log((1+epsi_m)/(1-epsi_m)); % works if 1
    %               uncertain stiffness & 1 uncertain mass: both uncertainties are
    %               proportional to the mean mass and stiffness matrices
    compare_Om2_moy=[Om2_det Om2_moy_theo Om2_moy_direct Om2_moy_PCE Om2_moy_PCE_a];
else
    compare_Om2_moy=[Om2_det Om2_moy_direct Om2_moy_PCE Om2_moy_PCE_a];
end

% ========================================
% Comparison with respect to MCS

% Om2: MCS
for i=1:n_ddl
    ecart_om2_MCS(i)=norm(Om2_u(i,:)-Om2_PCE(i,:))/norm(Om2_u(i,:))*100;
end

% eigenvectors 
for k=1:n_ddl
    
    for i=1:n_ddl  % all the eigenvectors have a unit first element
        Psi_u(k).Psi_u_norm1(i,:)=Psi_u(k).Psi_u(i,:)./Psi_u(k).Psi_u(1,:);
        Psi_u_PC(k).Psi_u_norm1(i,:)=Psi_u_PC(k).Psi_u(i,:)./Psi_u_PC(k).Psi_u(1,:);
    end

    Psi_alea_moy(:,k)=mean(Psi_u(k).Psi_u_norm1,2);
    Psi_alea_std(:,k)=std(Psi_u(k).Psi_u_norm1,1,2);
    
    Psi_alea_PC_moy(:,k)=mean(Psi_u_PC(k).Psi_u_norm1,2)';
    Psi_alea_PC_std(:,k)=std(Psi_u_PC(k).Psi_u_norm1,1,2)';

    ecart_Psi_fro(k)=norm(Psi_u_PC(k).Psi_u_norm1-Psi_u(k).Psi_u_norm1,'fro')/norm(Psi_u(k).Psi_u_norm1,'fro')*100;
end


% ========================================
% FRF

% damping
if exist('d_det','var')==0
    if exist('xi_det','var')==0
        disp('No damping is given.')
        d_det=0.*Om_det;
    else
        d_det=2*xi_det.*Om_det;
    end
end

% FRF

% Please select the dofs related to the FRF
i_dof1=1;
j_dof2=1;

% Construction of the FRFs
FRF_MCS_direct=zeros(npt_MCS_direct,length(w));
FRF_MCS_PC=zeros(npt_MCS_direct,length(w));
for i_u= 1 : npt_MCS_direct
    clear Vec_MCS_direct Vec_MCS_PC 
    for i_mod=1: n_ddl
        Vec_MCS_direct(:,i_mod)=Psi_u(i_mod).Psi_u(:,i_u) ;
        Vec_MCS_PC(:,i_mod)=Psi_u_PC(i_mod).Psi_u(:,i_u) ; 
    end
    if exist('D','var')==1
        d_MCS_direct=diag(Vec_MCS_direct'*D*Vec_MCS_direct)./M_mod_u(:,i_u);  % to be verified: never tested
        d_MCS_PC=diag(Vec_MCS_PC'*D*Vec_MCS_PC)./M_mod_PCE(:,i_u);            % to be verified: never tested
    elseif exist('xi_det','var')==1
        d_MCS_direct=2*xi_det.*Om_u(:,i_u);
        d_MCS_PC=2*xi_det.*Om_PCE(:,i_u);
    else
        damp=zeros(n_ddl,1);
        d_MCS_direct=zeros(n_ddl,1);
        d_MCS_PC=zeros(n_ddl,1);       
    end
        
    FRF_MCS_direct(i_u,:)=FRF_modale_damping(i_dof1,j_dof2,Vec_MCS_direct,Om_u(:,i_u),w,d_MCS_direct,M_mod_u(:,i_u));
    FRF_MCS_PC(i_u,:)=FRF_modale_damping(i_dof1,j_dof2,Vec_MCS_PC,Om_PCE(:,i_u),w,d_MCS_PC,M_mod_PCE(:,i_u));
end
FRF_mean_abs_MCS_direct=abs(mean(FRF_MCS_direct,1));
FRF_mean_abs_MCS_PC=abs(mean(FRF_MCS_PC,1));
FRF_std_abs_MCS_direct=abs(std(FRF_MCS_direct,1,1));
FRF_std_abs_MCS_PC=abs(std(FRF_MCS_PC,1,1));


for i_u=1:npt_MCS_direct
    ecart_FRF_MCS_PCE(i_u)=norm(FRF_MCS_direct(i_u,:)-FRF_MCS_PC(i_u,:))/norm(FRF_MCS_direct(i_u,:))*100;
    ecart_ABS_FRF_MCS_PCE(i_u)=norm(abs(FRF_MCS_direct(i_u,:))-abs(FRF_MCS_PC(i_u,:)))/norm(abs(FRF_MCS_direct(i_u,:)))*100;
    ecart_ARG_FRF_MCS_PCE(i_u)=norm(angle(FRF_MCS_direct(i_u,:))-angle(FRF_MCS_PC(i_u,:)))/norm(angle(FRF_MCS_direct(i_u,:)))*100;    
end
[Max_ecart_FRF,Indice_Max_ecart_FRF]=max(ecart_FRF_MCS_PCE);
[Max_ecart_ABS_FRF,Indice_Max_ecart_ABS_FRF]=max(ecart_ABS_FRF_MCS_PCE);
[Max_ecart_ARG_FRF,Indice_Max_ecart_ARG_FRF]=max(ecart_ARG_FRF_MCS_PCE);
[Max_ecart_FRF Indice_Max_ecart_FRF/1000 Max_ecart_ABS_FRF Indice_Max_ecart_ABS_FRF/1000 Max_ecart_ARG_FRF Indice_Max_ecart_ARG_FRF/1000];


t_simul;

%% ========================================
% Errors

ecart_om2_MCS 
ecart_Psi_fro

Err_FRF_mean=norm(FRF_mean_abs_MCS_direct-FRF_mean_abs_MCS_PC)/norm(FRF_mean_abs_MCS_direct)*100
Err_FRF_std=norm(FRF_std_abs_MCS_direct-FRF_std_abs_MCS_PC)/norm(FRF_std_abs_MCS_direct)*100

Err_FRF_MCS=norm(FRF_MCS_direct-FRF_MCS_PC,'fro')/norm(FRF_MCS_direct,'fro')*100


%% ========================================
% Figures 

if courbes==1

    val_n_ddl=int2str(n_ddl);
    val_n_uncer=int2str(n_uncer);
    val_i_dof1=int2str(i_dof1);
    val_j_dof2=int2str(j_dof2);
    val_deg=int2str(PC_order);
    val_npt_tot=int2str(npt_tot);
    val_npt_MCS=int2str(npt_MCS_direct);
    val_epsi_MK=int2str(100*epsi_m(1));    % we suppose delta_M=delta_K, which is true for the examples of the paper
    
    disp('----- Figures  ------')
    
    fname='arial';
    fsize=16;
    fname_title='arial';
    fsize_title=14;
    fname_legend='arial';
    fsize_legend=12;
    fname_axis='arial';
    fsize_axis=12;
    fig_x=400;fig_y=420;fig_wi_x=350;fig_wi_y=250;
    
    % ================
    % Figures FRF - mean
    
    h_frf_moy=figure(10*method+1);
    semilogy(w,FRF_mean_abs_MCS_direct,'k','linewidth',1)
    hold on
    semilogy(w,FRF_mean_abs_MCS_PC,'r--','linewidth',2)
    hold off
    xlabel('\omega (rad/s)','FontName',fname,'fontsize',fsize);
    ylabel('|FRF|_{mean}','FontName',fname,'fontsize',fsize);
    leg=legend('direct mean FRF','PCE mean FRF');
    set(leg,'Location', 'Best','FontName',fname_legend,'fontsize',fsize_legend)
    
    grid
    title('Non-intrusive PCE','FontName',fname_title,'fontsize',fsize_title);
    set(h_frf_moy,'position',[[fig_x fig_y-0*fig_wi_y fig_wi_x fig_wi_y]])
    
    ha=gca;set(ha,'linewidth',1.5,'FontName',fname_axis,'FontSize',fsize_axis,'Box','on');
    
    if imprime==1
        nom_fic_fig=strcat('figure\fig\NI-FRF-',pdf_law,'-moy_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        savefig(nom_fic_fig)
        nom_fic_eps=strcat('figure\eps\NI-FRF-',pdf_law,'-moy_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_eps,'-depsc')
        nom_fic_jpg=strcat('figure\jpg\NI-FRF-',pdf_law,'-moy_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_jpg,'-djpeg')
        nom_fic_pdf=strcat('figure\pdf\NI-FRF-',pdf_law,'-moy_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_pdf,'-dpdf')
    end
    
    % ================
    % Figures FRF - std
    
    h_frf_std=figure(method*10+2);
    semilogy(w,FRF_std_abs_MCS_direct,'k','linewidth',1)
    hold on
    semilogy(w,FRF_std_abs_MCS_PC,'r--','linewidth',2)
    hold off
    xlabel('\omega (rad/s)','FontName',fname,'fontsize',fsize);
    ylabel('|FRF|_{std}','FontName',fname,'fontsize',fsize);
    leg=legend('direct std FRF','PCE std FRF');
    set(leg,'Location', 'Best','FontName',fname_legend,'fontsize',fsize_legend)
    
    grid
    title('Non-intrusive PCE','FontName',fname_title,'fontsize',fsize_title);
    set(h_frf_std,'position',[[fig_x fig_y-1.4*fig_wi_y fig_wi_x fig_wi_y]])
    
    ha=gca;set(ha,'linewidth',1.5,'FontName',fname_axis,'FontSize',fsize_axis,'Box','on');
    
    if imprime==1
        nom_fic_fig=strcat('figure\fig\NI-FRF-',pdf_law,'-std_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        savefig(nom_fic_fig)
        nom_fic_eps=strcat('figure\eps\NI-FRF-',pdf_law,'-std_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_eps,'-depsc')
        nom_fic_jpg=strcat('figure\jpg\NI-FRF-',pdf_law,'-std_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_jpg,'-djpeg')
        nom_fic_pdf=strcat('figure\pdf\NI-FRF-',pdf_law,'-std_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_pdf,'-dpdf')
    end
    
    % ================
    % Figures FRF - worst case
    
    I_rand=[Indice_Max_ecart_ABS_FRF];
    val_I_rand=int2str(I_rand);
    
    h_frf=figure(method*10+7);
    subplot(2,1,1)
    semilogy(w, abs(FRF_MCS_direct(I_rand,:)),'k','linewidth',1)    
    hold on
    semilogy(w,abs(FRF_MCS_PC(I_rand,:)),'r--','linewidth',2)
    hold off
    xlabel('\omega (rad/s)','FontName',fname,'fontsize',fsize);
    ylabel('|FRF|','FontName',fname,'fontsize',fsize);
    axis([0 max_w 1e-2 1e2]);
    grid
    title('Non-intrusive PCE','FontName',fname_title,'fontsize',fsize_title);
    
    subplot(2,1,2)
    plot(w, angle(FRF_MCS_direct(I_rand,:)),'k','linewidth',1)    %FRF_MCS_PC(i_u,:)
    hold on
    plot(w,angle(FRF_MCS_PC(I_rand,:)),'r--','linewidth',2)
    hold off
    xlabel('\omega (rad/s)','FontName',fname,'fontsize',fsize);
    ylabel('angle(FRF)','FontName',fname,'fontsize',fsize);
    grid
    set(h_frf,'position',[[fig_x+0.8*fig_wi_x fig_y-0.7*fig_wi_y 1.5*fig_wi_x 1.5*fig_wi_y]])
    
    ha=gca;set(ha,'linewidth',1.5,'FontName',fname_axis,'FontSize',fsize_axis,'Box','on');
    
    if imprime==1
        nom_fic_fig=strcat('figure\fig\NI-FRF-',pdf_law,'-ech_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_I_rand-',val_I_rand,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        savefig(nom_fic_fig)
        nom_fic_eps=strcat('figure\eps\NI-FRF-',pdf_law,'-ech_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_I_rand-',val_I_rand,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_eps,'-depsc')
        nom_fic_jpg=strcat('figure\jpg\NI-FRF-',pdf_law,'-ech_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_I_rand-',val_I_rand,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_jpg,'-djpeg')
        nom_fic_pdf=strcat('figure\pdf\NI-FRF-',pdf_law,'-ech_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_i_',val_i_dof1,'_j_',val_j_dof2,'_I_rand-',val_I_rand,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_pdf,'-dpdf')
    end
    
    %===================
    % PDF frequencies
    
    n_max_om=min(n_ddl,10);
    h_pdf=figure(10*(method)+3);
    clf
    % hold on
    for numero_om=1:n_max_om;
        clear sig_direct sig_PCE pdf_sig_direct pdf_sig_PCE
        
        n_pt_hist=250;                        % number of points in each bin
        n_hist=fix(npt_MCS_direct/n_pt_hist); % histogram number of points
        
        sig_direct=Om_u(numero_om,:);
        sig_PCE=Om_PCE(numero_om,:);
        clear pdf_sig
        pdf_sig_direct=ksdensity(sig_direct);
        pdf_sig_PCE=ksdensity(sig_PCE);
        
        h1=plot(linspace(min(abs(sig_direct)),max(abs(sig_direct)),length(pdf_sig_direct)),pdf_sig_direct,'k','linewidth',2);
        hold on
        h2=plot(linspace(min(abs(sig_PCE)),max(abs(sig_PCE)),length(pdf_sig_PCE)),pdf_sig_PCE,'r--','linewidth',2);
    end
    hold off
    xlabel('frequency (rad/s)','FontName',fname,'FontSize',fsize)
    ylabel('PDF ','FontName',fname,'FontSize',fsize)
    leg=legend('direct MCS','PCE MCS');
    set(leg,'Location', 'Best','FontName',fname_legend,'fontsize',fsize_legend)
    grid
    title('Non-intrusive PCE','FontName',fname_title,'fontsize',fsize_title);
    set(h_pdf,'position',[[fig_x-1.1*fig_wi_x fig_y-0*fig_wi_y fig_wi_x fig_wi_y]])
    
    ha=gca;set(ha,'linewidth',1.5,'FontName',fname_axis,'FontSize',fsize_axis,'Box','on');
    
    if imprime==1
        nom_fic_fig=strcat('figure\fig\NI-PDF-',pdf_law,'-freq_pr_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        savefig(nom_fic_fig)
        nom_fic_eps=strcat('figure\eps\NI-PDF-',pdf_law,'-freq_pr_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_eps,'-depsc')
        nom_fic_jpg=strcat('figure\jpg\NI-PDF-',pdf_law,'-freq_pr_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_jpg,'-djpeg')
        nom_fic_pdf=strcat('figure\pdf\NI-PDF-',pdf_law,'-freq_pr_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_id-',val_npt_tot,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
        print(nom_fic_pdf,'-dpdf')
    end
    
end