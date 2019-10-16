function  nbcell = GeometricInfo( geo )
% Construct link from one cell to all neighbour cells.
% Compute normal distance from one cell to neighbour cells.

%% Read Meshdata
%{
[...
                            ncells, ~, ~, ~, ...
                            IC, FC, ~, ...
                            sn, snsign, ~, ~, ...
                            nface, ~, ~, ...
                            link_cell_to_face, ...
                            link_face_to_cell ...
                            ] =   GeometryProperties{:};
%}   

%% Read Meshdata
%ncells = GeometryProperties.ncells;

%link_cell_to_face = GeometryProperties.link_cell_to_face;
%link_face_to_cell = GeometryProperties.link_face_to_cell; 

%% Compute nbcell and df
 %                      
for icell = 1:geo.ncells
    for ifc = 1:geo.nface(icell)
        gfc = geo.link_cell_to_face(icell,ifc);
        nbcell(icell,ifc) = setdiff(geo.link_face_to_cell(gfc,:),icell);
    end
end

end

