function [x, final_res, iter] = pcg_inverse(A, b, tol, m_max, K1, K2, x0)
% SOLVEPCG   Preconditioned Conjugate Gradients method for solving
%    a system of linear equations Au = f.

if nargin < 3, tol = 1e-06; end
if nargin < 4, m_max = 1000; end


% initial parameters
x = x0;
r = b - A * x;
p = K2 * (K1 * r);

% initial iteration
iter = 0;norm_f = norm(b);

% start computing
while( (norm(r)/norm_f > tol) && (iter < m_max))
    
  q = A * p; % heavy load
  
  q2 = p' * q;
  
  alpha = (r' * p) / q2;
  
  x = x + alpha * p;
  r = r - alpha * q;
  
  z = K2 * (K1 * r); % heavy load
  
  beta = (z' * q) / q2;
  
  p = z - beta * p;
  
  iter = iter + 1;

end

final_res = norm(r)/norm_f;

end