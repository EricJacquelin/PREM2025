function [w_final, ind_final, L_final] = sparse_vb(X, y, r)
% -------------------------------------------------------------------------
% Finction for generating sparse matrix by variational Bayesian inference
% based on ARD (Automatic relevance determination)
% Input:
% X = Basis matrix (Number of column represents number of terms in the
% polynomial
% y = Response vector (Same number of rows as X with 1 column)
% r = Tuning parameter for ARD (If r is greater then less terms are
% discarded in each iteration and vice versa)
% Output:
% w = Coefficient vector with final number of terms in the polynomial
% ind = Indices of finally selected number of terms in the polynomial
% -------------------------------------------------------------------------

[N, D] = size(X);
ind = 1:1:D;

k = 1;
siz(k) = D;
L_final = -inf;
while D>1
    [w, V, invV, logdetV, an, bn, E_a, L] = vb_linear_fit_ard(X, y);
    ard = inv(diag(E_a));
    ard1 = diag(ard);
    ard_ln = log(ard1);
    t_ard = min(ard_ln)+((max(ard_ln)-min(ard_ln))/r);
    
    ind1 = find(ard_ln<=t_ard);
    X(:,ind1) = [];
    w(ind1) = [];
    [N, D] = size(X);
    ind(ind1) = [];
    k = k+1;
    siz(k) = length(ind);
    if L > L_final
        w_final = w;
        ind_final = ind;
        L_final = L;
    end
%     y_pred = X*w;
%     mse = sum((y-y_pred).^2)/N;
%     if mse>1e-3 || isequal(siz(k),siz(k-1))
%         return
%     end
end