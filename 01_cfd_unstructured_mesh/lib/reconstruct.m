function gradField = reconstruct(gamma, field, operator, o_offdiag, mesh)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if isempty(gamma), gamma = mesh.lf_weighted(~mesh.is_bface,:); end

%% Calculate weighted field
weighted_field = gamma.*(field(mesh.neighbour(~mesh.is_bface)) - field(mesh.owner(~mesh.is_bface)));

%% Construct corrector
weighted_field = repmat(weighted_field,2,1);
corrector = [accumarray(o_offdiag, weighted_field(:,1),[mesh.ncells,1]); ...
             accumarray(o_offdiag, weighted_field(:,2),[mesh.ncells,1])];

%% Reconstruct gradient
gradField = reshape(operator*corrector,[],2);

end

