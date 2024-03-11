function phi = gpuInterpPhi(U, owner, neighbour, is_bface, Sf, link_bface_to_face, nfaces )
% High performance phi interpolation function

%% Initialize phi
phi = zeros(nfaces,1,'gpuArray');

%% built-in operation

phi(~is_bface) = 0.5*((U.field(owner(~is_bface),1) + U.field(neighbour(~is_bface),1)).*Sf(~is_bface,1) + ...
                      (U.field(owner(~is_bface),2) + U.field(neighbour(~is_bface),2)).*Sf(~is_bface,2));

%% MEX interp function

%phi(~mesh.is_bface) = interpPhi(U.field',mesh.owner(~mesh.is_bface),mesh.neighbour(~mesh.is_bface),mesh.Sf(~mesh.is_bface,:)');


%% Correct boundary conditions
% Compute correction
ifc = link_bface_to_face;
is_dirchlet = U.boundary.dirchlet;
is_neumann = U.boundary.neumann;

ucorr_x = U.boundary.field(:,1).*is_dirchlet + U.field(owner(ifc),1).*is_neumann;
ucorr_y = U.boundary.field(:,2).*is_dirchlet + U.field(owner(ifc),2).*is_neumann;
operator = ucorr_x.*Sf(ifc,1) + ucorr_y.*Sf(ifc,2);

% Correct boundary condition
phi(ifc) = operator;

end