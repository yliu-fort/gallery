function leastSquareOperator = reconstructOperator(mesh, gpu_flag)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2, gpu_flag = false; end
assert(islogical(gpu_flag), 'function::diffusionOperator: Unrecongnized device_flag type.\n');

%% Construct least-square gradient reconstruction operator

        leastSquareOperator_tmp = cell(mesh.ncells,1);
        
        for icell = 1:mesh.ncells
            ifc = find(mesh.link_face_to_cell(:,icell));
            
            LS = mesh.lf(ifc,:);
            LS_weighted = mesh.lf_weighted(ifc,:); % weighted lf
            
            leastSquareOperator_tmp{icell} = inv(LS_weighted'*LS);
        end
        
        leastSquareOperator_tmp{1} = sparse(leastSquareOperator_tmp{1});
        leastSquareOperator_tmp2 = blkdiag(leastSquareOperator_tmp{:});
        
        % rearrange orders
        leastSquareOperator = [leastSquareOperator_tmp2(1:2:end,[1:2:end 2:2:end]);...
                               leastSquareOperator_tmp2(2:2:end,[1:2:end 2:2:end])];

%% gpu computing support

if gpu_flag
    leastSquareOperator = gpuArray(leastSquareOperator);
    fprintf('grad-reconstruct operator has been transfered to gpuArray successfully.\n');
end
     
end

