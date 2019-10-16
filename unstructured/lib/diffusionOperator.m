function [diffOperator,crossDiffOperator] = diffusionOperator(nu, mesh, gpu_flag)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3, gpu_flag = false; end
assert(islogical(gpu_flag), 'function::diffusionOperator: Unrecongnized device_flag type.\n');
%%
 tf = zeros(mesh.nfaces,2);
for ifc = 1:mesh.nfaces
    node_list = find(mesh.link_face_to_node(:,ifc));
    %c1 = verts(link_face_to_node(iface,1),:);
    %c2 = verts(link_face_to_node(iface,2),:);
    c1 = mesh.verts(node_list(1),:);
    c2 = mesh.verts(node_list(2),:);
    tf(ifc,:) = c2 - c1;
end

%tf = geo.verts(geo.link_face_to_node(:,2),:) - geo.verts(geo.link_face_to_node(:,1),:);

%%

diffOperator = nu*mesh.areaf./mesh.df; %Df(gfc,:) 05/28/18 removed '-'
crossDiffOperator = -nu*dot(tf,mesh.lf,2)./mesh.df; % remove ./geo.areaf.^2.*diffOperator after "dot(tf,lf,2)" 06/03/18

%% gpu computing support

if gpu_flag
    diffOperator = gpuArray(diffOperator);
    fprintf('diffusion operator has been transfered to gpuArray successfully.\n');
    crossDiffOperator = gpuArray(crossDiffOperator);
    fprintf('cross-diffusion operator has been transfered to gpuArray successfully.\n');
end

end