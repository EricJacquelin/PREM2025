%% 9th February 2024
%
% Moment matrices of a product of (n-1) PCs Psi(\Xi) times \Xi_i are evaluated
% for n=3 and n=4
%
% if i_u>0
% M3_sp(i_u+1,i).mom3(j,k) = <xi_i_u Psi(j) Psi(k)>
% M4_sp(i_u+1,i).mom4(j,k) = <xi_i_u Psi(j) Psi(k) Psi(l)>
% if i_u=0
% M3_sp(0+1,i).mom3(j,k) = <Psi(j) Psi(k)>
% M4_sp(0+1,j).mom4(j,k) = <Psi(j) Psi(k) Psi(l)>
%
% and Psi(j)= psi(j_1) x ... psi(j_n_uncer)=Psi_[j_1 j_2 ... j_n_uncer]
% [j_1 j_2 ... j_n_uncer]=deg_PC(j) calculated in data_RM_PC_INT.m
% psi(j_r) is function of \xi_r
%
% Therefore
%
% M3_sp(i_u,i).mom3(j,k) = <\xi_i psi(j_i) psi(k_i)> x \Prod_{r=1,r<>i}^n_uncer   <psi(j_r) psi(k_r)> 
% From the recurrence relation: \xi_i psi(j_i) = A_j_i psi(j_i + 1) - B_j_i psi(j_i) + C_j_i psi(j_i - 1)
% <\xi_i psi(j_i) psi(k_i) >=A_j_i <psi(j_i + 1) psi(k_i)> -
%                            B_j_i <psi(j_i)     psi(k_i)> + 
%                            C_j_i <psi(j_i - 1) psi(k_i)>
% Furthermore the orthonormality gives 
%    <psi(j_i) psi(k_i)> = delta_{j_i,k_i} (Kronecker's symbol)
%
% M4_sp(i_u,i).mom4(j,k) = <\xi_i psi(j_i) psi(k_i) psi(l_i)> x \Prod_{r=1,r<>i}^n_uncer   <psi(j_r) psi(k_r) psi(l_r)> 
% From the recurrence relation: \xi_i psi(j_i) = A_j_i psi(j_i + 1) - B_j_i psi(j_i) + C_j_i psi(j_i - 1)
% <\xi_i psi(j_i) psi(k_i) psi(l_i)>=A_j_i <psi(j_i + 1) psi(k_i) psi(l_i)> -
%                                    B_j_i <psi(j_i)     psi(k_i) psi(l_i)> + 
%                                    C_j_i <psi(j_i - 1) psi(k_i) psi(l_i)>

%% Coefficients of the recurrence equation

for i=0:PC_order
    if strcmp(pdf_law,'normal')==1
        A_PC(i+1)= sqrt(i+1);
        B_PC(i+1)=0;
        C_PC(i+1)=sqrt(i);
    elseif strcmp(pdf_law,'uniform')==1
        A_PC(i+1)= (i+1)/sqrt((2*i+1)*(2*i+3));
        B_PC(i+1)=0;
        C_PC(i+1)=i/sqrt((2*i-1)*(2*i+1));
    end
end

%% Coefficient of \xi in \psi_1(\xi), which is either a Hermite or a
% Legendre normalized polynomial of degree 1.

if strcmp(pdf_law,'normal')==1
    Coef1_PC=1;
elseif strcmp(pdf_law,'uniform')==1
    Coef1_PC=sqrt(3);
end


%% M3 moment matrix

clear M3_sp 
M3_sp(1).mom3=sparse(eye(PC_nb));    % i_u=0     
for i_u=1:n_uncer                    % i_u>0  
    ii_u=find(deg_PC(i_u+1,:)==1);   % i_u+1 corresponds to the row number in deg_PC 
                                     % ii_u corresponds to the uncertain parameter number related to the (i_u+1)-th deg_PC row
    clear  mom2
    for i=1:PC_nb
        clear i_ind
        i_ind=deg_PC(i,:);
        i_ind(i_u)='';
        for j=1:PC_nb
            mom2(i,j)=0;
            clear j_ind
            j_ind=deg_PC(j,:);
            j_ind(i_u)='';
            
            dif_ind_norm=norm(i_ind-j_ind);            
            if dif_ind_norm~=0   
                continue
            end
            % factor related to i_u - we work with indices, then with mom3_sp(i).mom3
            i_i_u=deg_PC(i,i_u);
            j_i_u=deg_PC(j,i_u);
            if i_i_u+1~=j_i_u && i_i_u-1~=j_i_u
                continue
            end     
            if i_i_u+1==j_i_u
                mom2(i,j)=mom2(i,j)+A_PC(i_i_u+1);
            end
            if i_i_u-1==j_i_u
                mom2(i,j)=mom2(i,j)+C_PC(i_i_u+1);
            end            

        end
        
    end
    M3_sp(ii_u+1).mom3=sparse(mom2); 
end


%% M4 moment matrix 


for i=1:PC_nb                        % i_u=0     
    M4_sp(1,i).mom4=PC_mom3_sp(i).mom3;
end

for i_u=1:n_uncer                    % i_u>0    
    ii_u=find(deg_PC(i_u+1,:)==1);   % i_u+1 corresponds to the row number in deg_PC
                                     % ii_u corresponds to the uncertain parameter number related to the (i_u+1)-th deg_PC row
    Ind_list=1:n_uncer;
    Ind_list(i_u)='';
    
    for i=1:PC_nb
        for j=1:PC_nb
		j;
            for k=1:PC_nb
			k;
                
                % factor related to i_u - we work with indices, then with mom3_sp(i).mom3
                i_i_u=deg_PC(i,ii_u);
                j_i_u=deg_PC(j,ii_u);
                k_i_u=deg_PC(k,ii_u);
                
                if i_i_u>0
                    moment_4_ana_i_u=A_PC(i_i_u+1)*full(mom3_sp(i_i_u+2).mom3(j_i_u+1,k_i_u+1))+...
                        C_PC(i_i_u+1)*full(mom3_sp(i_i_u+0).mom3(j_i_u+1,k_i_u+1));
                else
                    moment_4_ana_i_u=1/Coef1_PC*full(mom3_sp(2).mom3(j_i_u+1,k_i_u+1));
                end
                if moment_4_ana_i_u==0
                    mom4(j,k)=0;
                    continue
                end
                
                % product of factors not related to i_u - we work with numbers, then with mom3_sp(deg_PC(i,r)+1).mom3
                for r=Ind_list
                    moment_3_ana(r)=mom3_sp(deg_PC(i,r)+1).mom3(deg_PC(j,r)+1,deg_PC(k,r)+1);
                    if moment_3_ana(r)==0
                        moment_4_ana_0=0;
                        break
                    end
                end
                moment_4_ana_0=prod(moment_3_ana(Ind_list),2) ;
                
                
                mom4(j,k)=moment_4_ana_0*moment_4_ana_i_u;
            end
        end
        M4_sp(ii_u+1,i).mom4=sparse(mom4); 
    end
end


