function Fn = SurfNormFlux(U, Boundary, Dirchlet, geo )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%{
[...
                            ~, nfaces, ~, ~, ...
                            ~, ~, ~, ...
                            sn, ~, ~, areaf, ...
                            ~, ~, ~, ...
                            ~, ...
                            link_face_to_cell, ...
                            ~, ...
                            ~, ...
                            ~, ...
                            link_face_to_bface ...
                            ] =   GeometryProperties{:};
                    %}    
for ifc = 1:geo.nfaces
    gcell = nonzeros(geo.link_face_to_cell(ifc,:));
    bface = geo.link_face_to_bface(ifc);
    is_dirchlet = any(Dirchlet == bface);

    if bface ~= 0
        Fn(ifc,1) = U(gcell,:)*(geo.sn(ifc,:)')*geo.areaf(ifc); %neumann
        if is_dirchlet
            Fn(ifc,1) = Boundary(bface,:)*(geo.sn(ifc,:)')*geo.areaf(ifc); % test
            % for irregular bc
        end
    else
        Fn(ifc,1) = 0.5*sum(U(gcell,:))*(geo.sn(ifc,:)')*geo.areaf(ifc);   %interior
    end

end

end
