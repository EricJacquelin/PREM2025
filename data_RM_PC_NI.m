% July 20th 2023
%
% Random modes derived with PC with a Non-intrusive method
% The PC are normalized <Phi_i^2>=1
%  
% n_uncer: number of uncertain parameters
% 
% Called by
% - random_mode_PC_NI_direct.m: MCS from NI-PC-direct
% - random_mode_PC_NI_modal.m: MCS from NI-PC-modal: PCE is also expanded on a modal expansion
% - random_mode_PC_NI_POD.m: MCS from NI-PC-POD: PCE is also expanded on a POD expansion
%
% The samples of the radom modes are calculated in 
% Random_Mode_samples.m
%




%=========================================================%
% data

% for the 3-dof system, npt_tot=100 and PC_order=5 give good results

npt_id=50  % total number of simulations, which are divided into 2sets: one for identification, and the other one for validation
PC_order=3 % maximum PC degree
n_para_maxi_term_PCE=PC_order;%2;          % maximum number of parameters that are multiplied in each term of the expansion 
n_para_maxi_term_PCE=min([n_para_maxi_term_PCE PC_order n_uncer]);

%=========================================================%
% sampling for identifying and validating the PCE

if npt_id>100
    npt_val=floor(npt_id/10);
elseif npt_id>19
    npt_val=10;
elseif npt_id>9
    npt_val=floor(npt_tot/2);
else
    error('Increase the number of samples')
end
npt_tot=npt_id+npt_val;

% samples drawn with the matlab LHS programs
z_u_PC(1:npt_id,:)=LHS_matlab(npt_id,n_uncer,pdf_law);  
z_u_PC(1+npt_id:npt_tot,:)=LHS_matlab(npt_val,n_uncer,pdf_law);  


%=========================================================%
% PC construction

clear deg_PC
deg_PC=degree_PC(n_uncer,PC_order,n_para_maxi_term_PCE);  
PC_nb=size(deg_PC,1);
nb(PC_order)=PC_nb;
if strcmp(pdf_law,'uniform')
    don_sample(1).N_para=n_uncer;
    don_sample(2).N_para=0;
elseif strcmp(pdf_law,'normal')
    don_sample(1).N_para=0;
    don_sample(2).N_para=n_uncer;
end
don_sample(1).law='uniform';
don_sample(2).law='normal';

  
%===========================================================
% random mode samples

Random_Mode_samples  % outputs: Om2_u_PC(:,i_u)   Om_u_PC(:,i_u)  M_mod_u_PC(:,i_u)    Psi_u_PC(i_mod).Psi_u(:,i_u)

%===========================================================
% sampling to identify the PC model

I_id=1:npt_id;
xx=z_u_PC(I_id,:);  
I_val=1+npt_id:npt_tot;
xx_val=z_u_PC(I_val,:);


%===========================================================
% PC : matrix Phi

Phi=Phi_def_norm(xx,don_sample,deg_PC);

%===========================================================
% initialisation SVB-ARD

A0=1e-2;
B0=1e-4;
C0=1e-2;
D0=1e-4;
e_ARD=5;          % influence the sparcity of the response ("low" e_ARD produces a (almost) non sparse PCE; "high" e_ARD produces a (very) sparse PCE
para_SVB_ini=[A0,B0,C0,D0];


