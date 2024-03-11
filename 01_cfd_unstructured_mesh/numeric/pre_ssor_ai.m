function [output1, output2] = pre_ssor_ai(A, w,preconditioner_flag)
% SSOR - AI preconditioner

% Preconditioner
switch preconditioner_flag
    case 'SSOR'
        % SSOR preconditioner
        D = diag(diag(A));L = triu(A, 1);D_s = (1/w)*D + L;
        K1 = 1/sqrt(2 - w)*D_s*sqrt(inv(1/w*D));
        K2 = K1';
        M = K1*K2;
    case 'GaussSeidel'
        % Gauss-Seidel preconditioner
        K1 = tril(A);
        K2 = inv(diag(diag(A)))*(K1');
    case 'SSOR-AI'
        % SSOR preconditioner
        D = diag(diag(A));L = tril(A, -1);
        iDw = inv((1/w)*D);
        K1 = sqrt(2 - w)*sqrt(iDw)*(spones(D) - L*iDw + (L*iDw)^2 );%- (L*iDw)^3);
        K2 = K1';
        M = K2*K1;
    otherwise
        K1 = speye(size(A));
        K2 = K1;
end

if nargout < 2, output1 = M;end
if nargout > 1, output1 = K1; output2 = K2;end

end