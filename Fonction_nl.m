% March 12 2024
%
% Nonlinear (quadratic) fonction of the unknowns Y \in R^{(P+1)(N+1) x 1}
% P+1 <=> PC_nb
% N <=> n_ddl
% Considering the notations in the paper Y <=> [Y^t a^t]^t
%

function fnl=Fonction_nl(Y,MM4_k,PC_nb)

n_ddl=round(size(Y,1)/PC_nb);

fnl=zeros(n_ddl*PC_nb,1);

O1=zeros((n_ddl-1)*PC_nb);
O2=zeros(PC_nb);

for l=1:n_ddl
    clear fnl_int
    for m=1:PC_nb
        fnl_int(m,1)=1/2*Y'*[O1                MM4_k(l,m).MM4'
                             MM4_k(l,m).MM4    O2           ]*Y;
    end 
    fnl((l-1)*PC_nb+1:l*PC_nb)=fnl_int;
end

