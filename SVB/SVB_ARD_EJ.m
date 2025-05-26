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


function [ak_sparse,Ind_sparse,Ak_sparse,Bk_sparse,Vk_sparse,L]=SVB_ARD_EJ(Phi,y,para_SVB_ini,e_ARD)

[N,Np]=size(Phi);

s=0;

Ind=1:1:Np;
L_fin=-inf;
while Np>1 
    
    s=s+1    ;

    clear invVk Vk ak Dk Lq ARD
        
    % Variational Bayesian step
    
    [Ak,Bk,~,~,Lq,ak,Vk,moy_alp]=SVB_EJ(Phi,y,para_SVB_ini);
    
    L(s)=Lq(end);
    
    % Sparse step
    
%     ARD(Inds,1)=1./moy_alp;
    ARD=1./moy_alp;    
	logARD=log(ARD);
    lnT_ARD=min(logARD)+(max(logARD)-min(logARD))/e_ARD;
	Non_Ind=find(log(ARD)<=lnT_ARD);
    
    
    if L(s)>L_fin
        L_fin=L(s);
        ak_sparse=ak;
        Ind_sparse=Ind;
        Ak_sparse=Ak;
        Bk_sparse=Bk;
        Vk_sparse=Vk;
    end
    
        Ind(Non_Ind)=[];
    Phi(:,Non_Ind)=[];
    Np=size(Phi,2);

   
end
        




    
    


