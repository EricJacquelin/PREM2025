% March 14, 2024
%
% Jacobian of the nonlinear function Fonction_nl.m, function of the unknowns 
% Y \in R^{(P+1)(N+1) x 1}
% P+1 <=> PC_nb
% N <=> n_ddl
% Considering the notations in the paper Y <=> [Y^t a^t]^t
%

function dfnl=Jacobian_nl(Y,MM4_k,PC_nb)
 
n_ddl=round(size(Y,1)/PC_nb);

O1=zeros((n_ddl-1)*PC_nb);
O2=zeros(PC_nb);

for l=1:n_ddl
    clear dfnl_int
    for m=1:PC_nb
        dfnl_int(m,:)=Y'*[O1               MM4_k(l,m).MM4'
                           MM4_k(l,m).MM4  O2           ];
    end 
    dfnl((l-1)*PC_nb+1:l*PC_nb,:)=dfnl_int;
end

