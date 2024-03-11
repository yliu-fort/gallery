function gradField = gpuReconstruct(gamma, field, operator, o_offdiag, owner, neighbour, is_bface)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Calculate weighted field
weighted_field = gamma(~is_bface,:).*(field(neighbour(~is_bface)) - field(owner(~is_bface)));

%% Construct corrector
weighted_field = repmat(weighted_field,2,1);
corrector = [accumarray(o_offdiag, weighted_field(:,1)); ...
             accumarray(o_offdiag, weighted_field(:,2))];

%% Reconstruct gradient
gradField = reshape(operator*corrector,[],2);

end

