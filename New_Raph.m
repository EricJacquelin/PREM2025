% December 15 2015
%
% Newton_Raphson to solve a nonlinear (quadratic) equation of 
% the unknowns Y \in R^{(P+1)(N+1) x 1}
% P+1 <=> PC_nb
% N <=> n_ddl
% Considering the notations in the paper Y <=> [Y^t a^t]^t
%

function [Y,norm_res]=New_Raph(A,MM4_k,bb,om2_det_k,Y0,PC_nb)

n_ddl=round(size(Y0,1)/PC_nb);

epsilon=1e-10;
Compt_max=100;
Compt=0;
Y=Y0;

fnl=Fonction_nl(Y,MM4_k,PC_nb);
dfnl=Jacobian_nl(Y,MM4_k,PC_nb);

res=A*Y-om2_det_k*fnl-bb;
norm_res(Compt+1)=norm(res);

while norm(res)>epsilon && Compt<Compt_max
    Compt=Compt+1;
    S=A-om2_det_k*dfnl;
    dY=-S\res;
    Y=Y+dY;

    fnl=Fonction_nl(Y,MM4_k,PC_nb);
    dfnl=Jacobian_nl(Y,MM4_k,PC_nb);

    res=A*Y-om2_det_k*fnl-bb;
    norm_res(Compt+1)=norm(res);
end
Compt;
if Compt>=Compt_max
    Compt
end