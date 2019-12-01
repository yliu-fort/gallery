function J = desaturated_rainbow(m)
%JET    Variant of HSV
%   JET(M), a variant of HSV(M), is an M-by-3 matrix containing
%   the default colormap used by CONTOUR, SURF and PCOLOR.
%   The colors begin with dark blue, range through shades of
%   blue, cyan, green, yellow and red, and end with dark red.
%   JET, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   See also HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 5.7.4.2 $  $Date: 2005/06/21 19:31:40 $

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end
n = ceil(m/2);

% s = m-[0 m/6 m/3 m/2 m*2/3 m*5/6 m];
% rs = [0 0 1/4 1 1 1 .5];
% gs = [0 0 1/4 1 1/4 0 0];
% bs = [.5 1 1 1 1/4 0 0];


% s = m-[0 m/8 m/4 m/2 m*3/4 m*7/8 m];
% rs = [0 0 3/4 1 1 1 .5];
% gs = [0 0 3/4 1 3/4 0 0];
% bs = [.5 1 1 1 3/4 0 0];

%f1 = .4;

%s = m-[0 m*.2 m*.35 m*.475 m*.5 m*.525 m*.65 m*.8   m];
%rs =  [0    0    f1   1     1    1    1    1  .5];
%gs =  [0    0    f1   1     1    1    f1   0   0];
%bs =  [.5   1    1    1     1    1    f1   0   0];

g = [
         0.0 ...
         0.27843137254900002 ...
         0.27843137254900002 ...
         0.85882352941200002;
         0.14299999999999999 ...
         0 ...
         0 ...
         0.36078431372500003;
         0.28499999999999998 ...
         0 ...
         1 ...
         1;
         0.42899999999999999 ...
         0 ...
         0.50196078431400004 ...
         0;
         0.57099999999999995 ...
         1 ...
         1 ...
         0;
         0.71399999999999997 ...
         1 ...
         0.38039215686299999 ...
         0;
         0.85699999999999998 ...
         0.419607843137 ...
         0 ...
         0;
         1 ...
         0.87843137254899994 ...
         0.30196078431399997 ...
0.30196078431399997];

s = m-m.*g(:,1);
rs =  g(:,2);
gs =  g(:,3);
bs =  g(:,4);

r = interp1(s,rs,1:m)';
g = interp1(s,gs,1:m)';
b = interp1(s,bs,1:m)';

% r = [ones(1,m-n) 1:-1/(n-1):0]';
% g = [0:1/((m-n)-1):1 1:-1/(n-1):0]';
% b = [0:1/((m-n)-1):1 ones(1,n) ]';
J = [r g b];

J = J(end:-1:1,:);
