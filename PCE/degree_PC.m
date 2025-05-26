% 09/10/2018
%
% function Coef=degree_PC(n_para,deg,n_lie)
%
% n_para: nb of parameters
% deg: truncature degree
% n_lie: maximal number of factors that must appeared in each terms of the PCE
%
% Coef: n_term x n_para matrix: Coef(i,j) give the degree of 1D-PC related
%       to the j-th variable, which arise in the i-th term of the PCE
%

function Coef=degree_PC(n_para,deg,n_lie)

ind_ini=(0:1:deg)';

% Calcul des exposants de chaque variable sans contrainte, sauf celle du
% degré maxi de chaque variable

ind=ind_ini;
for i=1:n_para-1
    clear ind_int
    ind_int=kron(ind_ini,ones(size(ind,1),1));   
    ind=[ind_int repmat(ind,(deg+1),1)];
    clear non_zeros
    non_zeros=find(sum(ind,2)>deg);
    ind(non_zeros,:)='';
end
n_cas=size(ind,1);

% Calcul des exposants de chaque variable du modèle avec les autres contraintes

count=0;
for i=1:n_cas
    d=ind(i,:);
    non_zeros=length(find(d));
    if and(sum(d)<=deg,non_zeros<n_lie+1)
        count=count+1;
        Coef(count,:)=d;
    end
end

%% =========================================================%
% deg_PC is reorganized so that, if p<q, degree(Psi(p))<=degree(Psi(q))

[sum_deg,I_deg]=sort(sum(Coef,2));
Coef=Coef(I_deg,:);

%% =========================================================%
% deg_PC is reorganized so that, for i=1:n_para, Psi(i+1)=psi_1(xi_i)
% (reminder: Psi(1)=1)

I_inter=n_para:-1:1;
Coef_int=Coef(:,I_inter);
Coef=Coef_int;
