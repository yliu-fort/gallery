function phi = ParallelSurfNormFlux(U, mesh )
% High performance phi interpolation function

%% Initialize phi
phi = zeros(mesh.nfaces,1);

%% built-in operation

phi(~mesh.is_bface) = 0.5*sum((U.field(mesh.owner(~mesh.is_bface),:) +...
                               U.field(mesh.neighbour(~mesh.is_bface),:))...
                               .*mesh.Sf(~mesh.is_bface,:),2);

%% MEX interp function

%phi(~mesh.is_bface) = interpPhi(U.field',mesh.owner(~mesh.is_bface),mesh.neighbour(~mesh.is_bface),mesh.Sf(~mesh.is_bface,:)');


%% Correct boundary conditions
% Compute correction
ifc = mesh.link_bface_to_face;
is_dirchlet = U.boundary.dirchlet;
is_neumann = U.boundary.neumann;

ucorr.x = U.boundary.field(:,1).*is_dirchlet + U.field(mesh.owner(ifc),1).*is_neumann;
ucorr.y = U.boundary.field(:,2).*is_dirchlet + U.field(mesh.owner(ifc),2).*is_neumann;
operator = ucorr.x.*mesh.Sf(ifc,1) + ucorr.y.*mesh.Sf(ifc,2);

% Correct boundary condition
phi(ifc) = operator;

end