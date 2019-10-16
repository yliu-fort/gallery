function phi_v = Cell2Verts( phi_c,geo )
% Super slow.(improved)
% Medium slow.
% Need to improve performance.

%% Read data


%% Interpolation from cell to vertices

for ivt = 1:geo.nnodes
    %[local_face,~] = find(link_face_to_node == ivt); %removed @ 7:41pm 12/02
    [local_cell,~] =  find(geo.link_cell_to_node == ivt);
    %local_cell = nonzeros(unique(link_face_to_cell(local_face,:))); %removed @ 7:41pm 12/02
    d_cell_to_node = 1./(sqrt((geo.IC(local_cell,1) - geo.verts(ivt,1)).^2 +...
                          (geo.IC(local_cell,2) - geo.verts(ivt,2)).^2));
    %d_cell_to_node = d_cell_to_node.^2; % optional method
    dtot_cell_to_node = sum(d_cell_to_node);
    phi_v(ivt,1) = sum(phi_c(local_cell,1).*d_cell_to_node/dtot_cell_to_node); % !avoid using dot divide
end

end