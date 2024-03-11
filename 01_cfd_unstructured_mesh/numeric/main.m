% Comparison of solvePCG with Matlab's pcg function for the poisson
% matrix, Jacobi- and Gauss-Seidel preconditioner.

% Author: Andreas Klimke
% Date: May 2003
% Version: 1.1

% Residual error tolerance tol is 1*10^(-5)
tol = 1e-5;

% Maximum number of iterations 
m_max = 1000;

% Compute for N grid points per axis direction
N = [10:5:50];

% index counter
index=1;

% Start loop
for n = N
  
  % Compute grid step width h
  h = 1/(n+1);
  
  % Initial start vector is the (n^2) x 1 zero vector 
  u_s = zeros(n^2,1);
  
  % Right-hand side is = 1, scaled according to grid resolution
  f = h^2 * ones(n^2,1);
  
  % Get the Poisson matrix
  A=gallery('poisson', n);

  D = diag(diag(A));
  C1 = tril(A);     % (since tril(A) = L+D )
  C2 = 0.25*(C1');  % (since inv(D) = 1/4*Identity for the Poisson
                    % matrix.) 
  C= C1*C2;
  E = speye(n^2);
  u = zeros(n^2,1);

  tic;
  % PCG with Jacobi
  [u, m_J(index)]=solvePCG(A, f, u_s, D, E, tol, m_max);
  time_J(index)=toc;
  tic;
  % PCG with symmetric Gauss-Seidel
  [u, m_GS(index)]=solvePCG(A, f, u_s, C1, C2, tol, ...
			    m_max);
  time_GS(index)=toc;
  tic;
  % PCG with symmetric Gauss-Seidel, no splitting of C
  [u, m_GS(index)]=solvePCG(A, f, u_s, C, E, tol, ...
			    m_max);
  time_GS2(index)=toc;
  tic;
  % For comparison: Matlab's PCG with Jacobi
  [u,flag,r,m_Matlab_J(index)] = pcg(A, f, tol, m_max, D, E);
  time_Matlab_J(index)=toc;
  tic;
  % For comparison: Matlab's PCG with Gauss-Seidel
  [u,flag,r,m_Matlab_GS(index)] = pcg(A, f, tol, m_max, C1, C2);
  time_Matlab_GS(index)=toc;
	
  index = index + 1;
end

subplot(2,1,1);
plot(N.^2, m_J, 'bo-', N.^2, m_GS, 'rs-', N.^2, m_Matlab_J, 'k.--', N.^2, ...
     m_Matlab_GS, 'kx--', 'LineWidth', 1);

xlabel('Number of unknowns');
ylabel('Number of iterations (e=1e-5)');
legend('solvePCG with Jacobi', ... 
       'solvePCG with symmetric Gauss-Seidel', ...
       'Matlab''s pcg with Jacobi', ...
       'Matlab''s pcg with symmetric Gauss-Seidel');

subplot(2,1,2);
plot(N.^2, time_J,'bo-', N.^2, time_GS, 'ms-', N.^2, time_GS2,'rx-', ...
		 N.^2, time_Matlab_GS, 'kd--', 'LineWidth', 1);

xlabel('Number of unknowns');
ylabel('Computation time (e=1e-5) [s]');
legend('solvePCG with Jacobi', ... 
       'solvePCG with symmetric Gauss-Seidel', ...
       'solvePCG mit s. GS, no splitting of preconditioner', ...
			 'Matlab''s pcg with symmetric GS');





