function PhysicalProperties = Preprocessing(filename, GeometryProperties )

%% Read User-defined physical property data

nu = 0.0001;

%% Read geometric data

ncells = GeometryProperties.ncells;
nfaces = GeometryProperties.nfaces;
nbfaces = GeometryProperties.nbfaces;

%% (mark for removal)Compute for geometry info only use meshdata
% df : compute by geometry info, no need to update

%nbcell = GeometricInfo( GeometryProperties ); % removed 05/30/2018

%% Find nearest cell at specific location
% for setting fixed point
%IND = Nearest( [0.0 0.0],IC,ncells);
[~,setRef.refCell] = max(sum(GeometryProperties.IC,2));
%setRef.refCell = 1;
setRef.refValue = 0;

%% Setting boundary

[ WALL,MOVINGWALL,INLET,OUTLET,FREESLIP,u_wall, v_wall ,u_inlet,v_inlet]...
    = BoundarySettings(filename, GeometryProperties);

%% Mount different type boundaries for passing into solver

dirchletbc = sparse([WALL;MOVINGWALL;INLET],1,1,nbfaces,1);
neumannbc = sparse([OUTLET;FREESLIP],1,1,nbfaces,1);

%% Initialization

U.field = zeros(ncells,2) + [u_inlet v_inlet];

%u = zeros(ncells,1); % + u_inlet;
%v = zeros(ncells,1);
p = zeros(ncells,1);
phi = zeros(nfaces,1);

U.boundary.field= zeros(nbfaces,2);
%U.boundary.y = zeros(nbfaces,1);
%u_bc = zeros(nbfaces,1);
%v_bc = zeros(nbfaces,1);

%NZMAX = ncells + size(GeometryProperties.owner,1)+size(GeometryProperties.neighbour,1)+1; % allocate space for sparse matrix

% Set boundary values
U.boundary.field(INLET,1) = u_inlet;
U.boundary.field(INLET,2) = v_inlet;
%U.boundary.y(INLET) = v_inlet;
U.boundary.field(MOVINGWALL,1) = u_wall;
U.boundary.field(MOVINGWALL,2) = v_wall;
%U.boundary.y(MOVINGWALL) = v_wall;
U.boundary.dirchlet = dirchletbc;
U.boundary.neumann = neumannbc;
% further implementation for moving wall required
% further implementation for free-slip required

%% Output

PhysicalProperties = {...
    nu;...
    U;...
    p;...
    phi;...
    INLET;...
    OUTLET;...
    WALL; ...
    MOVINGWALL;...
    setRef...
    };
end

