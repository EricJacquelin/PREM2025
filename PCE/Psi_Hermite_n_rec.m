% 13th June 2015
%
% Recursion definition of normalized Hermite polynomial
% h_i(x)=H_i(x)/sqrt(i!)
%

function Psi_H=Psi_Hermite_n_rec(ksi,order)

% Polynomials 

npt=length(ksi);

Psi=zeros(npt,53);

Psi(:,1)=ones(npt,1);
Psi(:,2)=ksi;
if order>1
    for i=3:order+1
        Psi(:,i)=(ksi.*Psi(:,i-1)-sqrt(i-2)*Psi(:,i-2))/sqrt(i-1);
    end
end

Psi_H=Psi(:,order+1);