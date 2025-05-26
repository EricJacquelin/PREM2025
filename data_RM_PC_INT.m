%% February 13th 2024
%
% Random modes derived with PC with an intrusive method
% The PC are normalized <Phi_i^2>=1
%  
% n_uncer: number of uncertain parameters
% 
% Called by
% - Random_Mode_MAIN_INTRUSIVE.m
%

%% =========================================================%
% data

% for the 3-dof system, npt_tot=100 and PC_order=5 give good results

PC_order=3% maximum PC degree
n_para_maxi_term_PCE=2; % maximum number of parameters that are multiplied in each term of the expansion 
n_para_maxi_term_PCE=min([n_para_maxi_term_PCE PC_order n_uncer]);

%% =========================================================%
% PC construction

clear deg_PC
% deg_PC=degree_PC_hyp(n_uncer,PC_order,1);  % there are several algo to find the PC set
deg_PC=degree_PC(n_uncer,PC_order,n_para_maxi_term_PCE);  % there are several algo to find the PC set.
PC_nb=size(deg_PC,1)
nb(PC_order)=PC_nb;
if strcmp(pdf_law,'uniform')
    don_sample(1).N_para=n_uncer;
    don_sample(2).N_para=0;
else
    don_sample(1).N_para=0;
    don_sample(2).N_para=n_uncer;
end
don_sample(1).law='uniform';
don_sample(2).law='normal';

%% Coefficient of \xi in \psi_1(\xi), which is either a Hermite or a
% Legendre normalized polynomial of degree 1.

if strcmp(pdf_law,'normal')==1
    Coef1_PC=1;
elseif strcmp(pdf_law,'uniform')==1
    Coef1_PC=sqrt(3);
end

