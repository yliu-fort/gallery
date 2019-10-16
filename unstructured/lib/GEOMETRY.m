function [ GeometryProperties ] = GEOMETRY( meshfile,assertion )
%% Author: Yuxuan Liu
% Meshfile should contain grid point coordinates stored as variable p 
% and a matrix representing the set of triangles that make up the
% triangulation which stored as variable t.
% Mesh type should be fully triangle.

if nargin < 2
    assertion = 'assert on';
end

%% load meshfile

addpath mesh
filename = meshfile;
ts = tic;
meshdata = load (filename);
verts = getfield (meshdata,'p');
link_cell_to_node =getfield (meshdata,'t');
clear meshdata
telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for extracting data from file.\n']);

%% Geometry information

nnodes = size(verts,1);
ncells = size(link_cell_to_node,1);

%% Link face to node -> faces
ts = tic;

link_face_to_node = [];
for icell = 1:ncells
    fcycle = [[link_cell_to_node(icell,end) link_cell_to_node(icell,1:end-1)];link_cell_to_node(icell,:)]';
    link_face_to_node = [link_face_to_node;fcycle];
end

% mark repeated surfaces
for iface = 1:size(link_face_to_node)
    rept = all(ismember(link_face_to_node,link_face_to_node(iface,:)),2);
    rept(iface) = 0;
    link_face_to_node(rept == 1,:) = 0;
end

% remove repeated surfaces
link_face_to_node(all(link_face_to_node,2) == 0,:) = [];
clear rept

%link_face_to_node = edges(TR); % !sorted list may break connectivity for
%unknown reason.

telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for creating link from face to node.\n']);
%% Geometry information

nfaces = size(link_face_to_node,1);
nface = 3*ones(ncells,1); % only tri-type mesh

%% (mark for removal)Link face to cell
ts = tic;

for iface = 1:nfaces
    [local_cell_1,~] = find(link_cell_to_node == link_face_to_node(iface,1));
    [local_cell_2,~] = find(link_cell_to_node == link_face_to_node(iface,2));
        local_cell = [intersect(local_cell_1,local_cell_2); 0 ]; % add zero
        local_cell = local_cell([1 2]); % 2D
    link_face_to_cell(iface,:) = local_cell';
end
% Existing of zero means this face located in boundary

telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for creating link from face to cell.\n']);

%% Determine owner and neighbour for each face
owner = zeros(nfaces,1);
neighbour = zeros(nfaces,1);
is_bface = zeros(nfaces,1);

for ifc = 1:nfaces
    owner(ifc,1) = link_face_to_cell(ifc,1);
    neighbour(ifc,1) = link_face_to_cell(ifc,2);
    is_bface(ifc,1) = 0;
    if neighbour(ifc,1) == 0
        neighbour(ifc,1) = 0;
        is_bface(ifc,1) = 1;
    end
end

nbfaces = size(nonzeros(is_bface),1);

%% (mark for removal)Link from bface to face

[link_bface_to_face,~] = find(is_bface);

%% Check if there is bface with no association to any faces
if strcmp(assertion ,'assert on')
Assertion1 = min(link_bface_to_face) ~= 0;
assert( Assertion1,...
    'GEOProcessor:Incomplete connectivity between faces and bfaces.')
clear Assertion1
end

%% Link cell to face
ts = tic;

for icell = 1:ncells
    [local_face,~] = find(owner == icell|neighbour == icell);
    link_cell_to_face(icell,:) = local_face';
end

telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for creating link from cell to face.\n']);
%% Cell center coordinates

for icell = 1:ncells
    IC(icell,:) = mean(verts(link_cell_to_node(icell,:),:));
end

%% Face center coordinates

for iface = 1:nfaces
    FC(iface,:) = mean(verts(link_face_to_node(iface,:),:));
end

%% Cell circumcenter coordinates
%{
for icell = 1:ncells
    P = verts(link_cell_to_node(icell,:),:); % get nodes associated with this cell
    A1 = [P(:,1).^2+P(:,2).^2 P(:,2) [1;1;1]];
    A2 = [P(:,1) P(:,1).^2+P(:,2).^2 [1;1;1]];
    A = [P(:,1) P(:,2) [1;1;1]];
    CC(icell,:) = [det(A1)/2/det(A) det(A2)/2/det(A)];
end
%}
%% Face area

for iface = 1:nfaces
    c1 = verts(link_face_to_node(iface,1),:);
    c2 = verts(link_face_to_node(iface,2),:);
    areaf(iface,1) = sqrt(sum((c1 - c2).^2));
end

%% Cell volumes
% Triangle

for icell = 1:ncells
    a = areaf(link_cell_to_face(icell,1));
    b = areaf(link_cell_to_face(icell,2));
    c = areaf(link_cell_to_face(icell,3));
    vol(icell,1)=0.25*sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a));
end

%% Surface normal

for iface = 1:nfaces
    c1 = verts(link_face_to_node(iface,1),:);
    c2 = verts(link_face_to_node(iface,2),:);
    t = (c2 - c1)/areaf(iface);
    sn(iface,:) = [t(2) -t(1)]; %2D
end

%% (mark for removal)Determine surface normal direction for each cell

for icell = 1:ncells
    for ifc = 1:nface(icell)
        gface = link_cell_to_face(icell,ifc);
        ic1 = link_face_to_cell(gface,1);
        if (icell == ic1)
            snsign(icell,ifc) = 1;
        else
            snsign(icell,ifc) = -1;
        end
    end
end

%% (05/28/2018) Pack geometry information
GeometryProperties.ncells = ncells;
GeometryProperties.nfaces = nfaces;
GeometryProperties.nbfaces = nbfaces;
GeometryProperties.nnodes = nnodes;
GeometryProperties.IC = IC;
GeometryProperties.FC = FC;
GeometryProperties.verts = verts;
GeometryProperties.sn = sn;
GeometryProperties.snsign = snsign;
GeometryProperties.vol = vol;
GeometryProperties.areaf = areaf;
GeometryProperties.nface = nface;
GeometryProperties.link_cell_to_face = link_cell_to_face;
GeometryProperties.link_face_to_cell = link_face_to_cell;
GeometryProperties.link_face_to_node = link_face_to_node;
GeometryProperties.link_cell_to_node = link_cell_to_node;
GeometryProperties.link_bface_to_face = link_bface_to_face;
GeometryProperties.owner = owner;
GeometryProperties.neighbour = neighbour;
GeometryProperties.is_bface = is_bface;

end

