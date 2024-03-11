function mesh = loadMesh(casename, path, debug)
%% read polyMesh

%% input management
if nargin < 1, casename = 'mesh'; end
if nargin < 2, path = 'polyMesh'; end
if nargin < 3, debug = 0; end
path = [casename '/' path '/'];

%% load polymesh
meshfields = {'boundary','faces','neighbour','owner','points'};

for i = meshfields
    tic
    mesh.(i{:}) = read_field_openfoam([path i{:}],debug);
    fprintf('mesh %s has been loaded.\n',i{:});
    toc
end

%% Build mesh cell index and information
i = [1:mesh.owner.length 1:mesh.neighbour.length];
j = [mesh.owner.field.value;mesh.neighbour.field.value];
mesh.cells.index = sparse(i,double(j),true);
mesh.cells.length = size(mesh.cells.index,2);

%% write mesh infos
mesh.info.nfaces = mesh.faces.length;
mesh.info.nnodes = mesh.points.length;
mesh.info.ncells = mesh.cells.length;
mesh.info.nInternalfaces = mesh.neighbour.length;
mesh.info.nBoundaryfaces = mesh.info.nfaces - mesh.info.nInternalfaces;

%% modify boundary patches
mesh.boundary.patches = mesh.boundary.field.patch;
mesh.boundary = rmfield(mesh.boundary,'field');
patch_names = fieldnames(mesh.boundary.patches);
for i = 1:mesh.boundary.length
    mesh.boundary.patches.(patch_names{i}) = setfield(mesh.boundary.patches.(patch_names{i}),'nFaces',str2num(mesh.boundary.patches.(patch_names{i}).nFaces));
    mesh.boundary.patches.(patch_names{i}) = setfield(mesh.boundary.patches.(patch_names{i}),'startFace',str2num(mesh.boundary.patches.(patch_names{i}).startFace));
end

%% output management
fprintf('mesh structure has been loaded.\n');

end