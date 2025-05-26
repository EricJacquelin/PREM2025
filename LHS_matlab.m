%  2nd May 2025
%
% Latin Hypercube Sampling 
%
% function S=LHS_matlab(npt,n_var,law,para);
%
% The Matlab's lhs routines are used to generate lhs samples for either a
% uniform or a normal distribution
%
% If  law=='uniform',   the S distribution is uniform
%           para(1)=-1  para(2)=1
% If  law=='normal',    the S distribution is normal
%           mu=zeros(1,n_var)  sigma=diag(ones(n_var,1)
%

function S=LHS_matlab(npt,n_var,law)


rng shuffle
if strcmp(law,'uniform')==1
    S=lhsdesign(npt,n_var)*2  - 1;
elseif strcmp(law,'normal')==1
    S=lhsnorm(zeros(1,n_var),diag(ones(n_var,1)),npt) ;
end


