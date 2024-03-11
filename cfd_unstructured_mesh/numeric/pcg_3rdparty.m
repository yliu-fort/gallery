function [x,res_tol,iter] = pcg_3rdparty(A, b, tol, m_max, C1, C2, x0)
% SOLVEPCG   Preconditioned Conjugate Gradients method for solving
%    a system of linear equations Au = f.

if isempty(C1), C1 = speye(size(A));end
if isempty(C2), C2 = speye(size(A));end

x = x0;
r = b - A * x;
p = C2 \ (C1 \ r);
norm_f = norm(b);

  
iter = 0;
while( (norm(r)/norm_f > tol) & (iter < m_max))
    
  q = A * p;

  alpha = (r' * p) / (p' * q);
  
  x = x + alpha * p;
  r = r - alpha * q;
  z = C2 \ (C1 \ r);
  
  beta = (z' * q) / (q' * p);
  
  p = z - beta * p;
  
  iter = iter + 1;
end

res_tol = norm(r)/norm_f;
