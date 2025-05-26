% 8 November 2018
%
% Sparse Variational Bayesian with Automatic Relevance Determination
%
% Baysian inference coupled with ARD: ONE term is rejected at each loop s and
% the variational lower bound is then evaluated. The model s with the
% maximum L is then selected
% 
% [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,y,para_SVB_ini,e_ARD)
% 
% Inputs 
%   Phi: N x Np matrix: evalutation of the Np PC for the N samples points of the experimental design
%   y: Nx1 vector: output of the model for the N samples points of the experimental design  
%   para_SVB_ini(1)= A0 parameter of a 1st Gamma distribution
%   para_SVB_ini(2)= B0 parameter of a 1st Gamma distribution
%   para_SVB_ini(3)= C0 parameter of a 2nd Gamma distribution
%   para_SVB_ini(4)= D0 parameter of a 2nd Gamma distribution
%   e_ARD: parameter which give the threshold of the ARD
% 
% Outputs
%   ak_sparse: identified parameters of the model= mean of the conditional posterior normal distribution of the model parameters
%   Ind_sparse: list of nonzeros model parameters
%   Ak_sparse,Bk_sparse: parameters of the posterior Gamma distribution of the model error variance
%   Vk_sparse: (factor of the) covariance matrix of the conditional posterior normal distribution of the model parameters
%   L: Variational Lower Bound


function [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ_V1(Phi,y,para_SVB_ini)

[~,Np_ini]=size(Phi);
cpt_max=max(200,Np_ini*2);

Np=Np_ini;
s=0;
Non_Inds=[];
Ind=1:1:Np;
Inds=1:1:Np;
Np=length(Inds);
while Np>1 && (s<cpt_max)
    
    s=s+1    ;
   
    clear invVk Vk ak Dk Lq ARD
    ak=zeros(Np_ini,1);
    Vk=zeros(Np_ini);
    moy_alp=zeros(Np_ini,1);
    ARD=Inf*ones(Np_ini,1);

    
    % Variational Bayesian step
    
    Phi_int=Phi(:,Inds);
    [Ak,Bk,~,~,Lq,aak,VVk,mmoy_alp]=SVB_EJ(Phi_int,y,para_SVB_ini);
    
    L(s)=Lq(end);
    Vk(Inds,Inds)=VVk;
    ak(Inds)=aak;
    moy_alp(Inds)=mmoy_alp;
    
    % Sparse step
    
    ARD(Inds)=1./moy_alp(Inds);
    [~,Non_Ind]=min(ARD);
    Non_Inds=union(Non_Inds,Non_Ind);
    Inds=setdiff(Ind,Non_Inds);
    resu(s).Inds=Inds;
    resu(s).ak=ak(Inds);
    resu(s).Vk=Vk(Inds,Inds);
    resu(s).Bk=Bk;
    resu(s).Ak=Ak;
    
    if Np==length(Inds) && abs(L(s))~=Inf
        length(Inds)
        s=cpt_max;
    end
     Np=length(Inds);
%     [s Np]
    
end
        
[~,II]=max(L);
ak_sparse=resu(II).ak;
Ind_sparse=resu(II).Inds;
Vk_sparse=resu(II).Vk;
Ak_sparse=resu(II).Ak;
Bk_sparse=resu(II).Bk;



    
    


