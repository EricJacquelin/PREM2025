% 13th June 2015
%
% Recursion definition of normalized Legendre polynomial
% h_i(x)=H_i(x) * sqrt((2i+1)/2)
%

function Psi_L=Psi_Legendre_n_rec(ksi,order)

% Polynomials 

npt=length(ksi);

% Psi=zeros(npt,53);

Psi(:,1)=ones(npt,1);
Psi(:,2)=ksi*sqrt(3);
if order>1
    for i=3:order+1
        Psi(:,i)=sqrt((2*i-3)*(2*i-1))/(i-1)*ksi.*Psi(:,i-1)-(i-2)/(i-1)*sqrt((2*i-1)/(2*i-5))*Psi(:,i-2);
    end
end

Psi_L=Psi(:,order+1);