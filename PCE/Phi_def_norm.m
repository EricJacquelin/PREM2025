% 9 October 2018
%
% Polynomial chaos evaluation
%
% Phi=Phi_def_norm(S,don_sample,deg_PC)
%
% S: npt x n_para matrix: S(:,i) : npt sample vector of the i-th RV
%
% don_sample(i): information on the samples that have been drawn
%                don_sample(i).N_para: number of parameters according to
%                                 don_sample(i).law
%                don_sample(1).law => U[-1,1] 
%                don_sample(2).law => N(0,1)
%
% deg_PC: n_term x n_para matrix  deg_PC(i,j): degree of the j-th PC, for i-th term
%
% Phi: npt x n_term matrix:
%
% This work is for any number of RV
%

function Phi=Phi_def_norm(S,don_sample,deg_PC)

npt=size(S,1);
[n_term,n_para]=size(deg_PC);

for i_term=1:n_term
    clear Phi_int d non_zeros non_zeros_1 non_zeros_2 ind_1 ind_2
    d=deg_PC(i_term,:);
    non_zeros=find(d);
    non_zeros_1=find(non_zeros<don_sample(1).N_para+1);
    non_zeros_2=find(non_zeros>don_sample(1).N_para);
    ind_1=non_zeros(non_zeros_1);
    ind_2=non_zeros(non_zeros_2);
    
    compt=0;
    if min(size(ind_1))>0 % il existe des para tirés selon une loi uniforme
        for i=ind_1
            compt=compt+1;
            Phi_int(:,compt)=Psi_Legendre_n_rec(S(:,i),d(i));
        end
    end
    if min(size(ind_2))>0 % il existe des para tirés selon une loi normale
        for i=ind_2
            compt=compt+1;
            Phi_int(:,compt)=Psi_Hermite_n_rec(S(:,i),d(i));
        end
    end
    if min(size(non_zeros))==0
        compt=compt+1;
        Phi_int(:,compt)=ones(npt,1);
    end
    Phi(:,i_term)=prod(Phi_int,2);
end
