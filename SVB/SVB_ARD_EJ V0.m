% 8 November 2018
%
% Sparse Variational Bayesian with Automatic Relevance Determination
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
%
%  ak_sparse is a Np sparse-vector: the non-zeros values are collected in
%  Ind_sparse


function [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,y,para_SVB_ini,e_ARD)

[N,Np_ini]=size(Phi);
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
    
    Phis=Phi(:,Inds);
    [Ak,Bk,~,~,Lq,aak,VVk,mmoy_alp]=SVB_EJ(Phis,y,para_SVB_ini);
    
    L(s)=Lq(end);
    Vk(Inds,Inds)=VVk;
    ak(Inds)=aak;
    moy_alp(Inds)=mmoy_alp;
    
    % Sparse step
    
%     ARD(Inds,1)=1./moy_alp;
    ARD(Inds)=1./moy_alp(Inds);    
	logARD=log(ARD(Inds,1));
    lnT_ARD=min(logARD)+(max(logARD)-min(logARD))/e_ARD;
	Non_Ind0=find(log(ARD)<=lnT_ARD);
	Non_Ind=intersect(Inds,Non_Ind0);
    Non_Inds=union(Non_Inds,Non_Ind);
    Inds=setdiff(Ind,Non_Inds);
    resu(s).Inds=Inds;
    resu(s).ak=ak(Inds);
    resu(s).Vk=Vk(Inds,Inds);
    resu(s).Bk=Bk;
    resu(s).Ak=Ak;
    
    if Np==length(Inds) && abs(L(s))~=Inf
        s
		s=cpt_max;
    end
    Np=length(Inds);
end
        
[~,II]=max(L);
ak_sparse=resu(II).ak;
Ind_sparse=resu(II).Inds;
Vk_sparse=resu(II).Vk;
Ak_sparse=resu(II).Ak;
Bk_sparse=resu(II).Bk;



    
    


