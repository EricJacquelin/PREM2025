% 16 gennaio 2019
% 
% Analisi delle componenti principali dei signali temporali (Proper orthogonal decomposition)
% => sig n_t x k matrice: k signali che sono discretizzati in tempo con n_t
%   punti
% => fac: factor for selecting n terms in POD: the selecting terms are such
%     their eigenvalues lambda_i is greater than max(lambda_i) x fac
%     default value fac=10^(-5)
% => crit: criterion for selecting the method of selection of n terms
%      "factor": the selected terms: lambda_i >= max(lambda_i) x fac
%      "energy": the selected terms: sum_1^n lambda_i >= fac x sum_all lambda_i
%      default value crit=factor; 
% 
% TO BE COMPLETED: I have not programmed the energy criterion
% 
% 
%   lambda: n piu grandi autovalori di la matrice sig'*sig
%   V: autovettori di la matrice sig'*sig che correspondono a lambda
%   a: proiezione di sig su V => V*a è un approssimazione di sig
%   sig_rec: =V*a è un approssimazione di sig
%   err: errore tra sig e sig_rec
%   ener_tot=energia nei signali
%

function [a,V,lambda,err,sig_rec,ener_tot]=POD_temp(sig,fac,crit)

%===========================================================

if exist('crit')==0
    crit='factor';
    if exist('fac')==0 
        fac=10^(-5);
    end
end
if exist('fac')==0 && strcmp(crit,'energy')==1
    error('You should enter a factor for this criterion')
end
if exist('fac')==0 
    fac=10^(-5);
end
% fac

% Eigenvalues of the correlation matrix
R=sig*sig';
[n_t,k]=size(sig);

[zz,xx]=eig(R);   
[xx_sort,Ind]=sort(abs(diag(xx)),'descend');
n=length(find(abs(xx_sort)/abs(max(xx_sort))>fac));

zz=zz(:,Ind);
V=zeros(n_t,n);
for i=1:n
    V(:,i)=zz(:,i)/norm(zz(:,i));
end

ener_tot=sum(xx_sort);
lambda=xx_sort(1:n);

a=zeros(k,n);
for j=1:n
    a(:,j)=sig'*V(:,j);
end
sig_rec=V*a';

err=norm(sig-sig_rec)/norm(sig)*100;
