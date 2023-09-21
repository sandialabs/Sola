function [x12,relres] = krylov_sqrt(A,G,b,maxiter,tol,varargin)
%
% [x12,relres] = krylov_sqrt(A,G,b,maxiter,tol,varargin)
%
% This function computes A^{1/2}b using Lanczos approach
% described in Algorithm 2.2 of
%   "Quantifying Uncertainties in Bayesian Linear Inverse Problems using
%   Krylov Subspace Methods" - Saibaba, Chung, and Petroske, 2018
%
% Inputs:
%   A (n x n) - Sparse matrix or funMat type
%   G (n x n) - Sparse matrix or funMat type. Preconditioner, such that G'*G approx A^{-1}
%   b (n x 1) - right hand side
%   maxiter   - maximum number of Lanczos iterations
%         tol - tolerance for stopping
%    varargin - test {'True', 'False'} Optional parameter to verify accuracy of Lanczos relationships

if nargin > 5
    test = varargin{1};
else
    test = 'False';
end

n = size(b,1);
nrmb = norm(b);

%Initialize Lanczos quantities
V = zeros(n,maxiter);
T = zeros(maxiter+1,maxiter+1);

%First step
vj = b/nrmb;
vjm1 = b*0;
beta = 0;
ykp = 0;
relres = zeros(maxiter,1);

for j = 1:maxiter
    V(:,j) = vj;
    wj = G*(A(G'*vj));
    alpha = wj'*vj;
    wj = wj - alpha*vj -beta*vjm1;
    beta = norm(wj);
    
    %Set vectors for new iterations
    vjm1 = vj;
    vj =    wj/beta;
    
    %Reorthogonalize vj (CGS2) % Change to something more sophisticated
    vj = vj - V(:,1:j)*(V(:,1:j)'*vj);  vj = vj/norm(vj);
    vj = vj - V(:,1:j)*(V(:,1:j)'*vj);  vj = vj/norm(vj);
    
    %Set the tridiagonal matrix
    T(j,j) = alpha;
    T(j+1,j) = beta; T(j,j+1) = beta;
    
    %Compute partial Lanczos solution
    Tk = T(1:j,1:j);    Vk = V(:,1:j);
    e1 = zeros(j,1);  e1(1) = nrmb;
    T12 = sqrtm(Tk);
    yk = T12*e1;
    
    %   relres(j) = norm(G*(A(G'*(Vk*(Tk\e1))))-b); % Lanczos residual
    %Check differences b/w successive iterations
    relres(j) = norm([ykp; 0] - yk)/norm(yk);
    if  relres(j) < tol
        relres = relres(1:j);
        break
    else
        ykp = yk;
    end
end
x12 = G \ (Vk*yk);


%Test the accuracy of Lanczos
if strcmp(test,'True')
    
    maxiter = size(relres,1);
    AG = 0*Vk;
    for i = 1:maxiter
        AG(:,i) = A(G'*Vk(:,i));
    end
    AVk = G*AG;
    
    figure, imagesc(log10(abs(Vk'*AVk -Tk))), colorbar
    norm(Vk'*AVk -Tk)
    norm(Vk'*Vk- eye(maxiter))
    ek = zeros(maxiter,1);  ek(end) = beta;
    norm(AVk - Vk*Tk - vj*ek')
end
end