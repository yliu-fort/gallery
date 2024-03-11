%% load meshfile
%% Convert 3d mesh to 2d (will eventually be removed)
mesh = loadMesh('mesh',filename);

%% get point index

point_index_2d = find(abs(mesh.points.field.z) < 1e-6);
point_index_toremove = find(abs(mesh.points.field.z) > 1e-6);

link_cell_to_node = mesh.faces.field.index;
[~,cell_index_toremove] = find(link_cell_to_node(point_index_toremove,:) ~= 0);
cell_index_toremove = unique(cell_index_toremove);

%% remove unused index from table

link_cell_to_node(point_index_toremove,:) = [];
link_cell_to_node(:,cell_index_toremove) = [];

%%

verts = [mesh.points.field.x(point_index_2d) mesh.points.field.y(point_index_2d)];

%% clean vars

clear point_index_2d point_index_toremove cell_index_toremove

%% read mesh -> points
%addpath mesh
%filename =  'rectangle_denser';
%ts = tic;
%meshdata = load (filename);
%verts = getfield (meshdata,'p');
%link_cell_to_node =getfield (meshdata,'t');
%clear meshdata
%telapsed = toc(ts);
%fprintf([num2str(telapsed),' sec elasped for extracting data from file.\n']);

%TR = triangulation(link_cell_to_node,verts);
%% Geometry information
nnodes = size(verts,1);
ncells = size(link_cell_to_node,2);

%% Link face to node -> faces
ts = tic;
edges = [];
for icell = 1:ncells
    node_list = find(link_cell_to_node(:,icell));
    if numel(node_list) == 3
        edge_list = [node_list node_list([2:end 1])];
    else
        if  numel(node_list) == 4
            
            tmp_points = verts(node_list,:);
            
            tmp_tf_1 = tmp_points(2,:) - tmp_points(1,:);
            tmp_tf_3 = tmp_points(3,:) - tmp_points(1,:);
            tmp_tf_4 = tmp_points(4,:) - tmp_points(1,:);
            tmp_sn = [tmp_tf_1(2) -tmp_tf_1(1)];
            tmp_1 = sign(dot(tmp_sn,tmp_tf_3))*sign(dot(tmp_sn,tmp_tf_4));
            if tmp_1 < 0
                edge_list =  [node_list node_list([3 4 2 1])];
            else
                
                tmp_tf_2 = tmp_points(3,:) - tmp_points(2,:);
                tmp_tf_1 = tmp_points(1,:) - tmp_points(2,:);
                tmp_tf_4 = tmp_points(4,:) - tmp_points(2,:);
                tmp_sn = [tmp_tf_2(2) -tmp_tf_2(1)];
                tmp_1 = sign(dot(tmp_sn,tmp_tf_1))*sign(dot(tmp_sn,tmp_tf_4));
                if tmp_1 < 0
                    edge_list =  [node_list node_list([2 4 1 3])];
                else
                    edge_list =  [node_list node_list([2 3 4 1])];
                end
            end
        else
            fprintf('Not implemented.\n');
        end
        
    end
    
    edges{icell,1} = edge_list;
end
link_face_to_node = cat(1, edges{:});
link_face_to_node = unique(sort(link_face_to_node,2),'rows');

clear edge_list tmp* edges
telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for creating link from face to node.\n']);

%% Check if faces are recognized correctly
%scatter(FC(:,1),FC(:,2),[],FC(:,1),'.'),axis equal
%hold on
%scatter(IC(:,1),IC(:,2),[],IC(:,1)),axis equal

%% Geometry information

nfaces = size(link_face_to_node,1);

%% make index ...
ts = tic;
% link cell to node
%j = [1:ncells 1:ncells 1:ncells];
%i = link_cell_to_node(:);
%link_cell_to_node = sparse(i,j,true);

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

tmp_transpose = link_face_to_cell';
for ifc = 1:nfaces
    tmp_cell_list = find(tmp_transpose(:,ifc));
    owner(ifc,1) = tmp_cell_list(1);
    if max(size(tmp_cell_list)) == 2
    neighbour(ifc,1) = tmp_cell_list(2);
    is_bface(ifc,1) = 0;
    else
        neighbour(ifc,1) = 0;
        is_bface(ifc,1) = 1;
    end
end
clear transpose
nbfaces = size(nonzeros(is_bface),1);



clear tmp_*
telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for creating owner and neighbour lists.\n']);

%% sort owner and neighbour list (renumber faces) @ 06/09/2018
[is_bface,tmp_index] = sort(is_bface);
owner = owner(tmp_index);
neighbour = neighbour(tmp_index);

link_face_to_node = link_face_to_node(:,tmp_index);
link_face_to_cell = link_face_to_cell(tmp_index,:);

clear tmp_*

%% (mark for removal)Link from bface to face

[link_bface_to_face,~] = find(is_bface);

%% Cell center coordinates
ts = tic;
IC = zeros(ncells,2);

for icell = 1:ncells
    node_list = find(link_cell_to_node(:,icell));
    IC(icell,:) = mean(verts(node_list,:));
end

%% Face center coordinates
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

%% Cell volumes (need further improve)
% Triangle
vol = zeros(ncells,1);

for icell = 1:ncells
    node_list = find(link_cell_to_node(:,icell));
    
    tmp_points = verts(node_list,:);
    tmp_tris = delaunay(tmp_points(:,1),tmp_points(:,2));

    tmp_tria = [];
    for itri = 1:size(tmp_tris,1)
        
    tmp_verts = tmp_points(tmp_tris(itri,:),:);
    tmp_edges =sqrt(sum((tmp_verts - tmp_verts([2:end 1],:)).^2,2));

    a =tmp_edges(1);
    b =tmp_edges(2);
    c =tmp_edges(3);
    s = (a + b + c)/2;
    
    tmp_tria(itri) = sqrt(s*(s - a)*(s - b)*(s - c));
    end
    
    vol(icell,1)=sum(tmp_tria);
end
clear tmp* a b c s

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

%% Surface normal * area

Sf = sn.*areaf;

%%
lf_weighted = lf./sqrt(sum(lf.^2,2)); 
df = dot(sn,lf,2);

assert(min(df) > 0, 'normal distance between owner and neighbour can not be negative or 0.\n')

telapsed = toc(ts);
fprintf([num2str(telapsed),' sec elasped for computing extra parameters.\n']);

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
GeometryProperties.lf = lf;
GeometryProperties.lf_weighted = lf_weighted;
GeometryProperties.df = df;
GeometryProperties.Sf = Sf;

clear -regexp ^[^Gfm]*