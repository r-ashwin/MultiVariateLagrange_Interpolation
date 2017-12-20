function Fhat = multiVariateLagrange(mu,F,muhat,ordr)
% MULTIVARIATELAGRANGE()
% *interpolates a parametric vector using polynomials in the Lagrange form* 
% Arguments:
%   mu - The parameter snapshots (Mxp). Each row represents one value of
%   the parameter (1xp)
%   F -  The snapshots of the vector (matrix) to be interpolated (NxM).
%   Each column is the vector corresponding to one row of mu
%   muhat - The parameter value at which interpolation is sought (1xp)
%   ordr - The order of Lagrange polynomial for interpolation (integer
%   scalar)
% Output:
%   Fhat - The interpolated vector (matrix) (Nx1)

% Add dependecy path
% addpath('/permn');
n = ordr; % order of interpolation
m = length(muhat);    % number of variables
rho = nchoosek(n+m,n);
if rho > length(mu)
    disp('-------------------------------------------------------');
    disp('Not enough snapshots! Try lowering the "ordr" argument.');
    disp('returning a vector of zeros.....');
    disp('-------------------------------------------------------');
    Fhat = zeros(length(F), 1);
    return;
end

ind = knnsearch(mu,muhat,'k',rho);
F = F';

% Generate A matrix and a vector
p = permn([0:n], m);
p = p(sum(p,2)<=n, :);
A = []; a = [];
for j = 1:length(p)
    A = [A, dotProd(mu(ind,:), p(j,:))];
    a = [a, dotProd(muhat,     p(j,:))];
end
   
for i = 1:size(F,2)
    alpha = linsolve(A,F(ind,i));   
    Fhat(i,1) = a*alpha;
end
    

function c = dotProd(A, b)
c = bsxfun(@power, A, b);
c = prod(c, 2);
end

end