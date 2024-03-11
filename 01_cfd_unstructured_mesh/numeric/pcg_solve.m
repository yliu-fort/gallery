function [x, res_tol, final_res, iter] = pcg_solve(A, b, tol, m_max,K1,K2, x0)
% SOLVEPCG   Preconditioned Conjugate Gradients method for solving
%    a system of linear equations Au = f.
if nargin < 2
    mxsize = A;
    A = gallery('poisson', mxsize);
    b = rand(mxsize^2,1);
end

if nargin < 3, tol = 1e-06; end
if nargin < 4, m_max = 1000; end
if nargin < 7, x0 = 0*b; end
        
% initial parameters
x = x0;
r = b - A * x;
z = K2 \ (K1 \ r);
p = z;

% initial iteration
iter = 0;norm_f = norm(b);
res_tol = norm(r)/norm_f ;

% start computing
while( (res_tol(end) > tol) && (iter < m_max))
    
  q = A * p;
  q2 = p' * q;
  
  alpha = (r' * p) / q2;
  
  x = x + alpha * p;
  r = r - alpha * q;
  
  z = K2 \ (K1 \ r);
  
  beta = (z' * q) / q2;
  
  p = z - beta * p;
  
  iter = iter + 1;
  res_tol(iter) = norm(r)/norm_f;
end

final_res = res_tol(end);

end