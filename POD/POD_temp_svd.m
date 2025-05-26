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
% 
% TO BE COMPLETED: I have not programmed the energy criterion
% 
%   lambda: n piu grandi autovalori di la matrice sig'*sig
%   V: autovettori di la matrice sig'*sig che correspondono a lambda
%   a: proiezione di sig su V => V*a è un approssimazione di sig
%   sig_rec: =V*a è un approssimazione di sig
%   err: errore tra sig e sig_rec
%   energie= somma di tutti i autovalori
%
% Si usa la svd della matrice sig
%

function [a,V,lambda,err,sig_rec,energie]=POD_temp_svd(sig,fac,crit)

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

%===========================================================

% svd of the signal matrix
[n_t,k]=size(sig);

[zz,xx,ww] = svd(sig,0);

[xx_sort,Ind]=sort(abs(diag(xx)),'descend');
n=length(find(abs(xx_sort)/abs(max(xx_sort))>fac));


zz=zz(:,Ind);
V=zeros(n_t,n);
for i=1:n
    V(:,i)=zz(:,i)/norm(zz(:,i));
end

lambda=xx_sort(1:n);
energie=sum(xx_sort.^2);

a=zeros(k,n);
for j=1:n
    a(:,j)=sig'*V(:,j);
end
sig_rec=V*a';

err=norm(sig-sig_rec)/norm(sig)*100;
