%  September 26th  2023
%
% Determining the FRF between 2 dof by a modal expansion
%
% function S=FRF_modale(i_dof1,j_dof2, Vec, om, w, d ,M_mod )
% 
% i_dof1,j_dof2: dofs with respect to which the FRF is calculated 
% Vec: modal vector matrix (n_ddl x n_ddl matrix)
% om: eigenfrequency vector (n_ddl x 1 vector)
% M_mod: modal mass vector; if modal mass vetor M_mod is not given, it is supposed that the
% d = (Vec' x D x Vec)/M_mod     where D is the damping matrix
%

function FRF=FRF_modale_damping(i_dof1,j_dof2, Vec, om, w, d ,M_mod )

ij=sqrt(-1);

n_ddl=size(Vec,1);
if exist('M_mod','var')==0
    M_mod=ones(n_ddl,1);
end


if exist('d','var')==1
    if length(d)==1
        damp=d.*ones(n_ddl,1);
    else
        damp=d;
    end
else
    damp=zeros(n_ddl,1);
end

om2=om.^2;

FRF=zeros(size(w));
for i=1:n_ddl
    FRF=FRF+Vec(i_dof1,i)*Vec(j_dof2,i)./(-w.^2+om2(i)+ij*damp(i)*w)/M_mod(i); 
end


