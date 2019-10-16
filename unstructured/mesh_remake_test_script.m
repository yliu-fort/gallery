%% load meshfile


%% read mesh -> points
addpath mesh
%filename =  'rectangle_denser';
ts = tic;
meshdata = load (filename);
verts = getfield (meshdata,'p');
link_cell_to_node =getfield (meshdata,'t');
clear meshdata
telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for extracting data from file.\n']);

TR = triangulation(link_cell_to_node,verts);
%% Geometry information
nnodes = size(verts,1);
ncells = size(link_cell_to_node,1);

%% Link face to node -> faces
ts = tic;
link_face_to_node = edges(TR); % !sorted list will break connectivity

telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for creating link from face to node.\n']);
%% Geometry information

nfaces = size(link_face_to_node,1);

%% make index ...
ts = tic;
% link cell to node
j = [1:ncells 1:ncells 1:ncells];
i = link_cell_to_node(:);
link_cell_to_node = sparse(i,j,true);

% link face to node
j = [1:nfaces 1:nfaces];
i = link_face_to_node(:);
link_face_to_node = sparse(i,j,true);

% link face to cell
link_face_to_cell = double(link_face_to_node')*double(link_cell_to_node) == 2;

telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for creating index.\n']);
%% Determine owner and neighbour for each face -> owner & neighbour
owner = zeros(nfaces,1);
neighbour = zeros(nfaces,1);
is_bface = zeros(nfaces,1);

ts = tic;

transpose = link_face_to_cell';
for ifc = 1:nfaces
    cell_list = find(transpose(:,ifc));
    owner(ifc,1) = cell_list(1);
    if max(size(cell_list)) == 2
    neighbour(ifc,1) = cell_list(2);
    is_bface(ifc,1) = 0;
    else
        neighbour(ifc,1) = 0;
        is_bface(ifc,1) = 1;
    end
end
clear transpose
nbfaces = size(nonzeros(is_bface),1);

telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for creating owner and neighbour lists.\n']);
%% (mark for removal)Link from bface to face

[link_bface_to_face,~] = find(is_bface);

%% Cell center coordinates

ts = tic;
%IC(:,1) = full((sum(link_cell_to_node.*verts(:,1))./(sum(link_cell_to_node)))');
%IC(:,2) = full((sum(link_cell_to_node.*verts(:,2))./(sum(link_cell_to_node)))');

IC = zeros(ncells,2);
for icell = 1:ncells
    node_list = find(link_cell_to_node(:,icell));
    IC(icell,:) = mean(verts(node_list,:));
end

%% Face center coordinates

%FC(:,1) = full((sum(link_face_to_node.*verts(:,1))./(sum(link_face_to_node)))');
%FC(:,2) = full((sum(link_face_to_node.*verts(:,2))./(sum(link_face_to_node)))');

FC = zeros(nfaces,2);
for ifc = 1:nfaces
    node_list = find(link_face_to_node(:,ifc));
    FC(ifc,:) = mean(verts(node_list,:));
end

%% Face area
areaf = zeros(nfaces,1);

for ifc = 1:nfaces
    node_list = find(link_face_to_node(:,ifc));

    c1 = verts(node_list(1),:);
    c2 = verts(node_list(2),:);
    areaf(ifc,1) = sqrt(sum((c1 - c2).^2));
end

%% Cell volumes
% Triangle
vol = zeros(ncells,1);

for icell = 1:ncells
    face_list = find(link_face_to_cell(:,icell));

    a =areaf(face_list(1));
    b =areaf(face_list(2));
    c =areaf(face_list(3));
    vol(icell,1)=0.25*sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a));
end

%% Surface normal
sn = zeros(nfaces,2);

for ifc = 1:nfaces
    node_list = find(link_face_to_node(:,ifc));

    c1 = verts(node_list(1),:);
    c2 = verts(node_list(2),:);
    t = (c2 - c1)/areaf(ifc);
    sn(ifc,:) = [t(2) -t(1)]; %2D
end

%% make sn point from owner to neighbour
lf = zeros(nfaces,2);

for ifc = 1:nfaces
    if is_bface(ifc)
        lf(ifc,:) = FC(ifc,:) - IC(owner(ifc),:);
    else
        lf(ifc,:) = IC(neighbour(ifc),:) - IC(owner(ifc),:);
    end
end

snsign = (dot(sn,lf,2)>0) - (dot(sn,lf,2)<0);
sn = sn.*snsign;

%%
lf_weighted = lf./sqrt(sum(lf.^2,2)); 
df = dot(sn,lf,2);
assert(min(df) > 0, 'normal distance between owner and neighbour can not be negative or 0.\n')

telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for computing extra parameters.\n']);

%% (mark for removal)Determine surface normal direction for each cell
%{
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
%}
%% (05/28/2018) Pack geometry information
GeometryProperties.ncells = ncells;
GeometryProperties.nfaces = nfaces;
GeometryProperties.nbfaces = nbfaces;
GeometryProperties.nnodes = nnodes;
GeometryProperties.IC = IC;
GeometryProperties.FC = FC;
GeometryProperties.verts = verts;
GeometryProperties.sn = sn;
GeometryProperties.vol = vol;
GeometryProperties.areaf = areaf;
GeometryProperties.link_face_to_cell = link_face_to_cell;
GeometryProperties.link_face_to_node = link_face_to_node;
GeometryProperties.link_cell_to_node = link_cell_to_node;
GeometryProperties.link_bface_to_face = link_bface_to_face;
GeometryProperties.owner = owner;
GeometryProperties.neighbour = neighbour;
GeometryProperties.is_bface = is_bface;
GeometryProperties.TR = TR;
GeometryProperties.lf = lf;
GeometryProperties.lf_weighted = lf_weighted;
GeometryProperties.df = df;