% 22th August 2023
%
% Simulation of an uncertain N-dof system; each rigidity, each mass is random: 
% Each random parameter follows a statistical distribution: the pdf is either
% uniform or normal 
%
% The random mode are derived from two MCS procedure
%
% - random_Mode_MCS_direct.m: direct MCS 
% - PCE metamodel+MCS
%
% The data are called in data_NdofSystem.m
% 

%

clear

disp('  ')
disp('----- New simulation ------')

ij=sqrt(-1);
imprime =0; % to save my plots into files

%% =========================================================%
% Path of the PCE toolbox

if exist('P')==1
    path(P);
else
    P=path;
end

path(path,'.\PCE');

%=====================================
% data

norme_Psi=1;    % =1: normalization with respect to the mass matrix => M_mod=eyes(n_ddl) ; =2: the 1st element of each eigenvector=1 

% === data file: to be changed when a new system is studied === %

data_3dofSystem_SA_2UncerPara     % 3 dof 2 uncertain parameters
% data_3dofSystem_SA_9UncerPara     % 3 dof 9 uncertain parameters: PC_order=3, n_para_maxi_term_PCE=2 => OK

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
% random modes: direct MCS
% Entrées : z_u_MCS(i_u,i_EM) avec i_u=1:npt_MCS_direct et i_EM=1:n_uncer
% Sorties : Psi_u(i_mod).Psi_u(:,i_u)
%           Om_u(:,i_u)                pour i_u de 1 à npt_MCS_direct

t_MCS_ini=cputime;
Random_Mode_MCS_direct
t_MCS=cputime-t_MCS_ini;

%=====================================
% random modes: PCE metamodel

% data for intrusive PCE

disp('----- PC construction ------')
t_ini=cputime;
data_RM_PC_INT

% Evaluation of the additional matrices

if n_ddl==3 && n_uncer==9  && strcmp(pdf_law,'normal')==1 && PC_order==3
    disp('----- Moment matrices & PCE matrices  ------')
    load Mat_3_4_sp_136_normal  % => M3_sp M4_sp
    disp('----- PCE matrices  ------')
    Mat_mode_raid_masse
else
    disp('----- 3tuple moment matrix  ------')
    Matrices_moments_triple_norm
    
    disp('----- Moment matrices  ------')
    Matrices_moments_norm
    
    % Uncertain mass and stiffness "modal matrices"
    disp('----- PCE matrices  ------')
    Mat_mode_raid_masse

end

% PCE coefficient

disp('----- PCE coefficient calculation  ------')

for k=1:n_ddl
    
    clear MM3_k MM4_k KM3_k
    KM3_k=KM3;
    KM3_k(:,(k-1)*PC_nb+1:k*PC_nb,:,:)='';
    MM3_k=MM3(:,(k-1)*PC_nb+1:k*PC_nb,:,:);
    if exist('MM4_sp','var')==0
        for l=1:n_ddl
            for m=1:PC_nb
                for n=1:n_ddl%-1
                    MM4_sp(l,m).MM4(:,(n-1)*PC_nb+1:n*PC_nb)=MM4(:,(n-1)*PC_nb+1:n*PC_nb,l,m);
                end
            end
        end
    end
    MM4_k=MM4_sp;
    for l=1:n_ddl
        for m =1:PC_nb
            MM4_k(l,m).MM4(:,(k-1)*PC_nb+1:k*PC_nb)='';
        end
    end
     
    clear A
    A=[KM3_k    -Om2_det(k)*MM3_k];
    
    Y0=zeros(n_ddl*PC_nb,1);
    [Y,norm_res]=New_Raph(A,MM4_k,B(:,k),Om2_det(k),Y0,PC_nb);
    
    PCE_Psi(k).lambda=reshape(Y(1:(n_ddl-1)*PC_nb),PC_nb,n_ddl-1); % PC_nb x n_ddl matrix
    PCE_Om2(k).ak=Y((n_ddl-1)*PC_nb+1:n_ddl*PC_nb);    

end

% PCE MCS

t_postprocess_ini=cputime;
Phi_simul=Phi_def_norm(z_u_MCS,don_sample,deg_PC);

Om2_PCE=zeros(npt_MCS_direct,n_ddl);
for k=1:n_ddl
    Om2_PCE_coef(:,k)=Phi_simul*PCE_Om2(k).ak;
    Om2_PCE(:,k)=Om2_det(k)*Om2_PCE_coef(:,k);
    
    lambda_k=Phi_simul*PCE_Psi(k).lambda;
    liste_n=1:1:n_ddl;
    liste_n(k)='';
    Psi_det_n=Psi_det(:,liste_n);
    Psi_u_PC(k).Psi_u= repmat(Psi_det(:,k),1,npt_MCS_direct);  %=Psi_det(:,k);  %
    for i=1:n_ddl
        Psi_u_PC(k).Psi_u(i,:)=Psi_u_PC(k).Psi_u(i,:)+Psi_det_n(i,:)*lambda_k';
    end
end
Om2_PCE=Om2_PCE';
Om_PCE=sqrt(Om2_PCE);

for i_u=1:npt_MCS_direct        
    Masse_u=zeros(n_ddl);
    for i_EM=1:n_uncer
    	Masse_u=Masse_u+z_u_MCS(i_u,i_EM)*M_u(i_EM).M;
    end    
    M_tot=M_det+Masse_u;
    for k=1:n_ddl
        M_mod_PCE(k,i_u)=Psi_u_PC(k).Psi_u(:,i_u)'*M_tot*Psi_u_PC(k).Psi_u(:,i_u);
    end
    if n_ddl==3 % to check orthogonality
        Mat_mod_PCE_1(1,i_u)=Psi_u_PC(2).Psi_u(:,i_u)'*M_tot*Psi_u_PC(3).Psi_u(:,i_u);
        Mat_mod_PCE_2(1,i_u)=Psi_u_PC(3).Psi_u(:,i_u)'*M_tot*Psi_u_PC(1).Psi_u(:,i_u);
        Mat_mod_PCE_3(1,i_u)=Psi_u_PC(1).Psi_u(:,i_u)'*M_tot*Psi_u_PC(2).Psi_u(:,i_u);    
    end
end
t_simul=cputime-t_ini;
    
%% ========================================
% Comparisons between direct MCS and PCE metamodel MCS

disp('----- Comparison ------')

% Om2: mean, std

for i=1:n_ddl
    Om2_moy_PCE_a(i,1)=Om2_det(i)*PCE_Om2(i).ak(1);
end

Om2_moy_PCE=mean(Om2_PCE,2);
Om2_std_PCE=std(Om2_PCE,1,2);

if n_uncer==2   &&  strcmp(pdf_law,'uniform')==1
    Om2_moy_theo=Om2_det/2/epsi_m*log((1+epsi_m)/(1-epsi_m)); % works if 1
    %               uncertain stiffness & 1 uncertain mass: both uncertainties are
    %               proportional to the mean mass and stiffness matrices
    compare_Om2_moy=[Om2_det Om2_moy_theo Om2_moy_direct Om2_moy_PCE Om2_moy_PCE_a]
else
    compare_Om2_moy=[Om2_det Om2_moy_direct Om2_moy_PCE Om2_moy_PCE_a]
end
        
% ========================================
% Comparison with respect to MCS
 
for i=1:n_ddl
    ecart_om2_MCS(i)=norm(Om2_u(i,:)-Om2_PCE(i,:),'fro')/norm(Om2_u(i,:),'fro')*100;
end
ecart_om2_MCS   


% eigenvectors 

for k=1:n_ddl
    
    for i=1:n_ddl % all the eigenvectors have a unit first element
        Psi_u(k).Psi_u_norm1(i,:)=Psi_u(k).Psi_u(i,:)./Psi_u(k).Psi_u(1,:);
        Psi_u_PC(k).Psi_u_norm1(i,:)=Psi_u_PC(k).Psi_u(i,:)./Psi_u_PC(k).Psi_u(1,:);
    end

    Psi_alea_moy(:,k)=mean(Psi_u(k).Psi_u_norm1,2);
    Psi_alea_std(:,k)=std(Psi_u(k).Psi_u_norm1,1,2);
    
    Psi_alea_PC_moy(:,k)=mean(Psi_u_PC(k).Psi_u_norm1,2)';
    Psi_alea_PC_std(:,k)=std(Psi_u_PC(k).Psi_u_norm1,1,2)';
    
    ecart_Psi_fro(k)=norm(Psi_u_PC(k).Psi_u_norm1-Psi_u(k).Psi_u_norm1,'fro')/norm(Psi_u(k).Psi_u_norm1,'fro')*100;
   
end
ecart_Psi_fro

%% ========================================
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

FRF_det=FRF_modale_damping(i_dof1,j_dof2,Psi_det,Om_det,w,d_det,M_mod_det);
FRF_abs_det=abs(FRF_det);

FRF_MCS_direct=zeros(npt_MCS_direct,length(w));
FRF_MCS_PC=zeros(npt_MCS_direct,length(w));
for i_u= 1 : npt_MCS_direct
    clear Vec_MCS_direct Vec_MCS_PC 
    for i_mod=1: n_ddl
        Vec_MCS_direct(:,i_mod)=Psi_u(i_mod).Psi_u(:,i_u) ;
        Vec_MCS_PC(:,i_mod)=Psi_u_PC(i_mod).Psi_u(:,i_u) ;
    end
    if exist('D','var')==1
        d_MCS_direct=diag(Vec_MCS_direct'*D*Vec_MCS_direct);  % to be verified: never tested
        d_MCS_PC=diag(Vec_MCS_PC'*D*Vec_MCS_PC);              % to be verified: never tested
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

Err_FRF_mean=norm(FRF_mean_abs_MCS_direct-FRF_mean_abs_MCS_PC)/norm(FRF_mean_abs_MCS_direct)*100
Err_FRF_std=norm(FRF_std_abs_MCS_direct-FRF_std_abs_MCS_PC)/norm(FRF_std_abs_MCS_direct)*100

Err_FRF_MCS=norm(FRF_MCS_direct-FRF_MCS_PC,'fro')/norm(FRF_MCS_direct,'fro')*100

for i_u=1:npt_MCS_direct
    ecart_FRF_MCS_PCE(i_u)=norm(FRF_MCS_direct(i_u,:)-FRF_MCS_PC(i_u,:))/norm(FRF_MCS_direct(i_u,:))*100;
    ecart_ABS_FRF_MCS_PCE(i_u)=norm(abs(FRF_MCS_direct(i_u,:))-abs(FRF_MCS_PC(i_u,:)))/norm(abs(FRF_MCS_direct(i_u,:)))*100;
    ecart_ARG_FRF_MCS_PCE(i_u)=norm(angle(FRF_MCS_direct(i_u,:))-angle(FRF_MCS_PC(i_u,:)))/norm(angle(FRF_MCS_direct(i_u,:)))*100;
end
[Max_ecart_FRF,Indice_Max_ecart_FRF]=max(ecart_FRF_MCS_PCE);
[Max_ecart_ABS_FRF,Indice_Max_ecart_ABS_FRF]=max(ecart_ABS_FRF_MCS_PCE);
[Max_ecart_ARG_FRF,Indice_Max_ecart_ARG_FRF]=max(ecart_ARG_FRF_MCS_PCE);
[Max_ecart_FRF Indice_Max_ecart_FRF/1000 Max_ecart_ABS_FRF Indice_Max_ecart_ABS_FRF/1000 Max_ecart_ARG_FRF Indice_Max_ecart_ARG_FRF/1000];

t_simul

%% ========================================
% Figures 

disp('----- Figures  ------')

val_n_ddl=int2str(n_ddl);
val_n_uncer=int2str(n_uncer);
val_i_dof1=int2str(i_dof1);
val_j_dof2=int2str(j_dof2);
val_deg=int2str(PC_order);
val_npt_MCS=int2str(npt_MCS_direct);
val_epsi_MK=int2str(100*epsi_m(1));    % we suppose delta_M=delta_K

fname='arial';
fsize=16;
fname_title='arial';
fsize_title=14;
fname_legend='arial';
fsize_legend=12;
fname_axis='arial';
fsize_axis=12;
fig_x=1000;fig_y=420;fig_wi_x=350;fig_wi_y=250;

%===================
% Figures FRF - mean

h_frf_moy=figure(1);
semilogy(w,FRF_mean_abs_MCS_direct,'k','linewidth',1)
hold on
semilogy(w,FRF_mean_abs_MCS_PC,'r--','linewidth',2)
hold off
xlabel('\omega (rad/s)','FontName',fname,'fontsize',fsize);
ylabel('|FRF|_{mean}','FontName',fname,'fontsize',fsize);
leg=legend('direct mean FRF','PCE mean FRF');
set(leg,'Location', 'Best','FontName',fname_legend,'fontsize',fsize_legend)

grid
title('Intrusive PCE','FontName',fname_title,'fontsize',fsize_title);
set(h_frf_moy,'position',[[fig_x fig_y-0*fig_wi_y fig_wi_x fig_wi_y]])

ha=gca;set(ha,'linewidth',1.5,'FontName',fname_axis,'FontSize',fsize_axis,'Box','on');

if imprime==1
    nom_fic_fig=strcat('figure\fig\INT-FRF-',pdf_law,'-moy_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    savefig(nom_fic_fig)
    nom_fic_eps=strcat('figure\eps\INT-FRF-',pdf_law,'-moy_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_eps,'-depsc')
    nom_fic_jpg=strcat('figure\jpg\INT-FRF-',pdf_law,'-moy_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_jpg,'-djpeg')
    nom_fic_pdf=strcat('figure\pdf\INT-FRF-',pdf_law,'-moy_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_pdf,'-dpdf')
end

%===================
% Figures FRF - std

h_frf_std=figure(2);
semilogy(w,FRF_std_abs_MCS_direct,'k','linewidth',1)
hold on
semilogy(w,FRF_std_abs_MCS_PC,'r--','linewidth',2)
hold off
legend('direct std FRF','PCE std FRF')
xlabel('\omega (rad/s)','FontName',fname,'fontsize',fsize);
ylabel('|FRF|_{std}','FontName',fname,'fontsize',fsize);
leg=legend('direct std FRF','PCE std FRF');
set(leg,'Location', 'Best','FontName',fname_legend,'fontsize',fsize_legend)

grid
title('Intrusive PCE','FontName',fname_title,'fontsize',fsize_title);
set(h_frf_std,'position',[[fig_x fig_y-1.4*fig_wi_y fig_wi_x fig_wi_y]])

ha=gca;set(ha,'linewidth',1.5,'FontName',fname_axis,'FontSize',fsize_axis,'Box','on');

if imprime==1
    nom_fic_fig=strcat('figure\fig\INT-FRF-',pdf_law,'-std_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    savefig(nom_fic_fig)
    nom_fic_eps=strcat('figure\eps\INT-FRF-',pdf_law,'-std_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_eps,'-depsc')
    nom_fic_jpg=strcat('figure\jpg\INT-FRF-',pdf_law,'-std_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_jpg,'-djpeg')
    nom_fic_pdf=strcat('figure\pdf\INT-FRF-',pdf_law,'-std_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_pdf,'-dpdf')
end

%===================
% Figures FRF - worst case

I_rand=[Indice_Max_ecart_ABS_FRF];
val_I_rand=int2str(I_rand);

h_frf=figure(7);
subplot(2,1,1)
semilogy(w, abs(FRF_MCS_direct(I_rand,:)),'k','linewidth',1)   
hold on
semilogy(w,abs(FRF_MCS_PC(I_rand,:)),'r--','linewidth',2)
hold off
xlabel('\omega (rad/s)','FontName',fname,'fontsize',fsize);
ylabel('|FRF|','FontName',fname,'fontsize',fsize);
axis([0 max_w 1e-2 1e2]);
grid
title('Intrusive PCE','FontName',fname_title,'fontsize',fsize_title);

subplot(2,1,2)
plot(w, angle(FRF_MCS_direct(I_rand,:)),'k','linewidth',1)    
hold on
plot(w,angle(FRF_MCS_PC(I_rand,:)),'r--','linewidth',2)
hold off
xlabel('\omega (rad/s)','FontName',fname,'fontsize',fsize);
ylabel('angle(FRF)','FontName',fname,'fontsize',fsize);
grid
set(h_frf,'position',[[fig_x+0.5*fig_wi_x fig_y-0.7*fig_wi_y 1.5*fig_wi_x 1.5*fig_wi_y]])

ha=gca;set(ha,'linewidth',1.5,'FontName',fname_axis,'FontSize',fsize_axis,'Box','on');

if imprime==1
    nom_fic_fig=strcat('figure\fig\INT-FRF-',pdf_law,'-ech_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_I_rand-',val_I_rand,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    savefig(nom_fic_fig)
    nom_fic_eps=strcat('figure\eps\INT-FRF-',pdf_law,'-ech_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_I_rand-',val_I_rand,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_eps,'-depsc')
    nom_fic_jpg=strcat('figure\jpg\INT-FRF-',pdf_law,'-ech_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_I_rand-',val_I_rand,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_jpg,'-djpeg')
    nom_fic_pdf=strcat('figure\pdf\INT-FRF-',pdf_law,'-ech_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_i_',val_i_dof1,'_j_',val_j_dof2,'_I_rand-',val_I_rand,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_pdf,'-dpdf')
end

%===================
% PDF frequencies

n_max_om=min(n_ddl,10);
h_pdf=figure(3);

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
title('Intrusive PCE','FontName',fname_title,'fontsize',fsize_title);

set(h_pdf,'position',[[fig_x-1.1*fig_wi_x fig_y-0*fig_wi_y fig_wi_x fig_wi_y]])

ha=gca;set(ha,'linewidth',1.5,'FontName',fname_axis,'FontSize',fsize_axis,'Box','on');

if imprime==1
    nom_fic_fig=strcat('figure\fig\INT-PDF-',pdf_law,'-freq_pr_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    savefig(nom_fic_fig)
    nom_fic_eps=strcat('figure\eps\INT-PDF-',pdf_law,'-freq_pr_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_eps,'-depsc')
    nom_fic_jpg=strcat('figure\jpg\INT-PDF-',pdf_law,'-freq_pr_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_jpg,'-djpeg')
    nom_fic_pdf=strcat('figure\pdf\INT-PDF-',pdf_law,'-freq_pr_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_pdf,'-dpdf')
end

%===================
% PDF FRF(freq_i)

i_fr_pr=3;
val_i_fr_pr=num2str(i_fr_pr);
i=min(find(fr>=fr_det(i_fr_pr)));   % freq_i is close to the 1st eigenfrequency
freq_i=fr(i);
val_freq=num2str(round(freq_i*2*pi,2));

clear sig_direct sig_PCE pdf_sig_direct pdf_sig_PCE

n_pt_hist=200;                        % number of points in each bin
n_hist=fix(npt_MCS_direct/n_pt_hist); % histogram number of points

sig_direct=abs(FRF_MCS_direct(:,i));
sig_PCE=abs(FRF_MCS_PC(:,i));
clear pdf_sig
pdf_sig_direct=ksdensity(sig_direct);
pdf_sig_PCE=ksdensity(sig_PCE);


h_pdf_frf=figure(4);
h1=plot(linspace(min(abs(sig_direct)),max(abs(sig_direct)),length(pdf_sig_direct)),pdf_sig_direct,'k','linewidth',1);
hold on
h2=plot(linspace(min(abs(sig_PCE)),max(abs(sig_PCE)),length(pdf_sig_PCE)),pdf_sig_PCE,'r--','linewidth',2);
hold off
nom_x=strcat('FRF(',val_freq,' rad/s) (m)');
nom_y=strcat('PDF FRF(',val_freq,' rad/s)');
xlabel(nom_x,'FontName',fname,'FontSize',fsize)
ylabel(nom_y,'FontName',fname,'FontSize',fsize)
leg=legend('direct MCS','PCE MCS');
set(leg,'Location', 'Best','FontName',fname_legend,'fontsize',fsize_legend)

title('Intrusive PCE','FontName',fname_title,'fontsize',fsize_title);
set(h_pdf_frf,'position',[[fig_x-1.1*fig_wi_x fig_y-1.4*fig_wi_y fig_wi_x fig_wi_y]])
hold off

ha=gca;set(ha,'linewidth',1.5,'FontName',fname_axis,'FontSize',fsize_axis,'Box','on');

if imprime==1
    nom_fic_fig=strcat('figure\fig\INT-PDF_FRF-',pdf_law,'-fr_pr_',val_i_fr_pr,'_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    savefig(nom_fic_fig)
    nom_fic_eps=strcat('figure\eps\INT-PDF_FRF-',pdf_law,'-fr_pr_',val_i_fr_pr,'_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_eps,'-depsc')
    nom_fic_jpg=strcat('figure\jpg\INT-PDF_FRF-',pdf_law,'-fr_pr_',val_i_fr_pr,'_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_jpg,'-djpeg')    
    nom_fic_pdf=strcat('figure\pdf\INT-PDF_FRF-',pdf_law,'-fr_pr_',val_i_fr_pr,'_',val_n_ddl,'_',val_n_uncer,'_deg-',val_deg,'_n_MCS-',val_npt_MCS,'_epsi-MK',val_epsi_MK);
    print(nom_fic_pdf,'-dpdf')
end

path=P;