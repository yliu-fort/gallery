function nearest = Nearest( coordinates,IC,ncells )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[~,ind] = sort(sum((IC - ones(ncells,1)*coordinates).^2,2));
nearest = ind(1);

end

