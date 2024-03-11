function PhysicalProperties = gpuPreprocessing(filename, GeometryProperties )

%% Read User-defined physical property data

nu = 1.8e-05;

%% Read geometric data

ncells = GeometryProperties.ncells;
nfaces = GeometryProperties.nfaces;
nbfaces = GeometryProperties.nbfaces;

%% Find nearest cell at specific location
[~,setRef.refCell] = max(sum(GeometryProperties.IC,2));
setRef.refValue = 0;

%% Setting boundary
[ WALL,MOVINGWALL,INLET,OUTLET,FREESLIP,u_wall, v_wall ,u_inlet,v_inlet]...
    = BoundarySettings(filename, GeometryProperties);

%% Mount different type boundaries for passing into solver
dirchletbc = full(sparse([WALL;MOVINGWALL;INLET],1,1,nbfaces,1));
neumannbc = full(sparse([OUTLET;FREESLIP],1,1,nbfaces,1));

%% Initialization

% initial fields
U.field = zeros(ncells,2) + [u_inlet v_inlet];
p = zeros(ncells,1,'gpuArray');
phi = zeros(nfaces,1,'gpuArray');

% initial boundary field
U.boundary.field = zeros(nbfaces,2);

% Set boundary values
U.boundary.field(INLET,1) = u_inlet;
U.boundary.field(INLET,2) = v_inlet;

U.boundary.field(MOVINGWALL,1) = u_wall;
U.boundary.field(MOVINGWALL,2) = v_wall;

U.boundary.dirchlet = dirchletbc;
U.boundary.neumann = neumannbc;
% further implementation for moving wall required
% further implementation for free-slip required

%% gpu support
U.field = gpuArray(U.field);
U.boundary.field = gpuArray(U.boundary.field);
U.boundary.dirchlet = gpuArray(U.boundary.dirchlet);
U.boundary.neumann = gpuArray(U.boundary.neumann);

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

