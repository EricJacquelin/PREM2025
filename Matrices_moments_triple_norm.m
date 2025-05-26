% 9th February 2024
%
%
% Triple Moment matrices are evaluated:
%
% PC_mom3(i).mom3(j,k) = <Psi(i) Psi(j) Psi(k)>
%
% and Psi(i)= psi(i_1) x ... psi(i_n_uncer)=Psi_[i_1 i_2 ... i_n_uncer]
% [i_1 i_2 ... i_n_uncer]=deg_PC(i) calculated in data_RM_PC_INT.m
%
% Therefore
% PC_mom3(i).mom3(j,k) = <psi(i_1) psi(j_1) psi(k_1)> x ... x  <psi(i_n_uncer) psi(j_n_uncer) psi(k_n_uncer)> 
%                      = \Prod_{r=1}^n_uncer   <psi(i_r) psi(j_r) psi(k_r)> 
%                      = \Prod_{r=1}^n_uncer   mom3_sp(i_r+1).mom3(j_r+1,k_r+1)
% 


clear mom3_sp
for i=0:PC_order+1   % we goes to index PC_order+1, as mom3_sp((PC_order+1)+1).mom3 is required to calculate M4(i_u,PC_nb).mom4
    clear mom3
    for j=0:PC_order
        for k=0:PC_order            
                s=(i+j+k)/2;
                if max([i j k])> s || s~=fix(s)
                    mom3(j+1,k+1)=0;
                else
                    if strcmp(pdf_law,'normal')==1
                        mom3(j+1,k+1)=sqrt(factorial(i)*factorial(j)*factorial(k))/...
                            (factorial(s-i)*factorial(s-j)*factorial(s-k));
                    elseif strcmp(pdf_law,'uniform')==1
                        As=factorial(2*s)/2^s/(factorial(s))^2;
                        Asm=factorial(2*(s-i))/2^(s-i)/(factorial((s-i)))^2;
                        Asn=factorial(2*(s-j))/2^(s-j)/(factorial((s-j)))^2;
                        Asp=factorial(2*(s-k))/2^(s-k)/(factorial((s-k)))^2;
                        mom3(j+1,k+1)=sqrt((2*i+1)*(2*j+1)*(2*k+1))/(2*s+1)*Asm*Asn*Asp/As;
                    end
                end            
        end
    end
    mom3_sp(i+1).mom3=sparse(mom3);
end

for i=1:PC_nb
    clear mom3
    for j=1:PC_nb
        for k=1:PC_nb
            for r=1:n_uncer
                moment_3_ana(r)=mom3_sp(deg_PC(i,r)+1).mom3(deg_PC(j,r)+1,deg_PC(k,r)+1);
                if moment_3_ana(r)==0                   
                    break
                end
            end
            mom3(j,k)=prod(moment_3_ana,2) ;
        end
    end
    PC_mom3_sp(i).mom3=sparse(mom3);
end


            
