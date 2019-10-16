function output = ParallelCrossDiffusionTerm(field, operator1, operator2, mesh)
% J(T,f) = jf*lf
% based on interpolated phi_v and tf in the face regardless whether the
% face is located on boundary

% Interpolate cell value to node, operator1 -> nnodes * ncells
tmp_var = operator1*field;

% Calculate diffusion correction term, operator2 -> nfaces * 1
[nodelist,~,~] = find(mesh.link_face_to_node);
%output.x = tmp_var.x(nodelist(1:2:end)) - tmp_var.x(nodelist(2:2:end));
%output.y = tmp_var.y(nodelist(1:2:end)) - tmp_var.y(nodelist(2:2:end));
output = tmp_var(nodelist(1:2:end),:) - tmp_var(nodelist(2:2:end),:);

output = output.*operator2;

end

