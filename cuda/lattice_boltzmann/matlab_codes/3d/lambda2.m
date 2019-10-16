function out = lambda2(a11,a12,a13,a21,a22,a23,a31,a32,a33)
J = [a11 a12 a13;a21 a22 a23;a31 a32 a33];
S = 0.5*(J + J.');
O = 0.5*(J - J.');
A = S.^2 + O.^2;

out = median(eig(A));
end