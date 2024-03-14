function c = coldhot(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUETecplot(M), is an M-by-3 matrix that defines a colormap.
%   The colors match Tecplot's Diverging Redblue plot.
%   REDBLUETecplot, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblueTecplot)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
% Inspired from https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap
%   Fernando Zigunov, 2018
if nargin < 1, m = size(get(gcf,'colormap'),1); end

%Colormap key colors from Tecplot
cm = [	[0,1,1],
        [0,0,1],
        [0,0,0.50196078431400004],
        [1,0,0],
        [1,1,0]];
        
%pos=[0;0.45;0.5;0.55;1];
pos=[0;0.495;0.5;0.505;1];

%Interpolates given the colormap positions
posDesired=linspace(0,1,m);
c = interp1(pos,cm,posDesired.');