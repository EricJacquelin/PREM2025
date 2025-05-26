% 8 November 2018
%
% Sparse Variational Bayesian 
% 
% Baysian inference coupled with ARD: terms are rejected if their ARD are
% lower than a threshold T_ARD 
% 
% [Ak,Bk,Ck,Dk,L,ak,Vk,moy_alp]=SVB_EJ(Phi,y,para_SVB_ini)
% 
% Inputs 
%   Phi: N x Np matrix: evalutation of the Np PC for the N samples points of the experimental design
%   y: Nx1 vector: output of the model for the N samples points of the experimental design  
%   para_SVB_ini(1)= A0 parameter of a 1st Gamma distribution
%   para_SVB_ini(2)= B0 parameter of a 1st Gamma distribution
%   para_SVB_ini(3)= C0 parameter of a 2nd Gamma distribution
%   para_SVB_ini(4)= D0 parameter of a 2nd Gamma distribution
% 
% Outputs
%   Ak,Bk: parameters of the posterior Gamma distribution of the model error variance
%   Ck,Dk: parameters of the posterior Gamma distribution of the model parameter variance
%   ak: identified parameters of the model= mean of the conditional posterior normal distribution of the model parameters
%   Vk: (factor of the) covariance matrix of the conditional posterior normal distribution of the model parameters
%   L: Variational Lower Bound
%   moy_alp: mean of the model parameter variance 

function [Ak,Bk,Ck,Dk,L,ak,Vk,moy_alp]=SVB_EJ(Phi,y,para_SVB_ini)

% Some initial quantities

[N,Np]=size(Phi);

A0=para_SVB_ini(1);
B0=para_SVB_ini(2);
C0=para_SVB_ini(3);
D0=para_SVB_ini(4);

PhiVkPhi=zeros(N,1);
PhiPhi=Phi'*Phi;
PhiY=Phi'*y;
n2y=norm(y)^2;

cpt_max=500;
T_L=1e-5;
cpt=0;

% == start of algo 1

k=0;

% Initialize

Ak=A0+N/2;
Bk=B0;
Ck=C0+1/2;
Dk=D0*ones(Np,1);

% moy_alp=Ck./Dk;
moy_alp=C0./Dk;  % by Jan Drugowitsch
moy_lam=diag(moy_alp);

% Variational Bayesian step

Lq(1)=-1e10;
Lq(2)=Lq(1)/(1+T_L);

while (abs(Lq(k+2)-Lq(k+1))>T_L*abs(Lq(k+2))) && (cpt<cpt_max)

    k=k+1;
    cpt=cpt+1;
    invVk=PhiPhi+moy_lam;
    Vk=inv(invVk);
    ak=Vk*PhiY;
    Bk=B0+1/2*(n2y-ak'*invVk*ak);
    moy_zetaa2=ak.^2*Ak/Bk+diag(Vk);
    Dk=D0+1/2*moy_zetaa2;
    moy_alp=Ck./Dk;
    moy_lam=diag(moy_alp);
    for i=1:N
        PhiVkPhi(i,1)=Phi(i,:)*Vk*Phi(i,:)';
    end
%     chol(Vk);
%     [chol_Vk,chol_p] = chol(Vk);
%     if    chol_p~=0
%         display('cho(Vk) => Vk is not positive definite')
%     end
% 	logdetVk=2 * sum(log(diag(chol_Vk)), 1);
    logdetVk=2 * sum(log(diag(chol(Vk))), 1);
    Lq(k+2)=-N/2*log(2*pi)+1/2*logdetVk-B0*Ak/Bk...
        -1/2*( Ak/Bk*norm(y-Phi*ak)^2 +sum(PhiVkPhi))...
        +gammaln(Ak)-Ak*log(Bk)+Ak-gammaln(A0)+A0*log(B0)...
        -Ck*sum(log(Dk))+Np*(1/2-gammaln(C0)+C0*log(D0)+gammaln(Ck));

%     if k==1
%        Lq(k+1)=Lq(k+2)-2*T_L;
%    end
end
L=Lq(3:end);

if cpt==cpt_max
    warning('The required precision of Lq is not reached')
end
if cpt==0
    warning('No Bayesian inference occured')
end

