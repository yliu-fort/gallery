function output = cell2nodeinterpOperator(mesh, gpu_flag)
%Return cell center to vertices field interpolation operator
%   Author: Yuxuan Liu
if nargin < 2, gpu_flag = false; end
assert(islogical(gpu_flag), 'function::diffusionOperator: Unrecongnized device_flag type.\n');

% calculating interp weights for each node
%transpose = mesh.link_cell_to_node';
%NZMAX = nnz(transpose);

% initialize field
%output = spalloc(mesh.nnodes, mesh.ncells,NZMAX);

%{
for ivt = 1:mesh.nnodes
    
    local_cell_list = find(transpose(:,ivt));
    operator = 1./(sqrt(...
        (mesh.IC(local_cell_list,1) - mesh.verts(ivt,1)).^2 ...
        +(mesh.IC(local_cell_list,2) - mesh.verts(ivt,2)).^2)...
        );
    

    output(ivt,local_cell_list) = operator./sum(operator);
    
end
%}

[i,j,~] = find(mesh.link_cell_to_node);

k = 1./sqrt(...
        (mesh.IC(j,1) - mesh.verts(i,1)).^2 ...
        +(mesh.IC(j,2) - mesh.verts(i,2)).^2);
    
kbar = accumarray(i,k);

output = sparse(i,j,k./kbar(i));

%% gpu computing support

if gpu_flag
    output = gpuArray(output);
    fprintf('cell2vert interp operator has been transfered to gpuArray successfully.\n');
end


end

