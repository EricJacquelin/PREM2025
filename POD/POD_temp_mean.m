% 16 gennaio 2019
% 
% Analisi delle componenti principali dei signali temporali (Proper orthogonal decomposition)
% => sig n_t x k matrice: k signali che sono discretizzati in tempo con n_t
%   punti
% 
%   lambda: n piu grandi autovalori di la matrice sig'*sig
%   V: autovettori di la matrice sig'*sig che correspondono a lambda
%   a: proiezione di sig su V => V*a è un approssimazione di sig
%   sig_rec: =V*a è un approssimazione di sig
%   err: errore tra sig e sig_rec
%   

function [a,V,lambda,sig_mean,err,sig_rec]=POD_temp(sig_ini)

%===========================================================

[n_t,k]=size(sig_ini);

sig_M=max(sig_ini,[],2);
sig_m=min(sig_ini,[],2);

sig_mean=0.5*(sig_M+sig_m);
sig=sig_ini-repmat(sig_mean,1,k);

% matrice de correlazione
R=sig*sig';

[zz,xx]=eig(R);   
[xx_sort,Ind]=sort(abs(diag(xx)),'descend');
n=length(find(abs(xx_sort)/abs(max(xx_sort))>10^(-5)));

zz=zz(:,Ind);
V=zeros(n_t,n);
for i=1:n
    V(:,i)=zz(:,i)/norm(zz(:,i));
end

lambda=xx_sort(1:n);

a=zeros(k,n);
for j=1:n
    a(:,j)=sig'*V(:,j);
end
sig_rec=sig_mean+V*a';

err=norm(sig-sig_rec)/norm(sig)*100;
