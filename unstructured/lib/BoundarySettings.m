function [ WALL,MOVINGWALL,INLET,OUTLET,FREESLIP,u_wall, v_wall ,u_inlet,v_inlet] = BoundarySettings( filename, GeometryProperties )
% Setting boundaries here.
% only 1 boundary class can be assigned for moving wall, and the boundary
% is supposed to be straight line. i.e lid-driven problem
% free slip requires further implementation.
% default settings for cavity flow, pipe flow, circular flow and duct flow
%% Capture&Sort boundary

[ WBC,EBC,SBC,NBC,INT ] = CatchBoundary( filename,1e-7, GeometryProperties);
%% Set boundary conditions here

WALL = [INT;SBC;EBC;WBC];
MOVINGWALL = [NBC]; u_wall = 0.0;   v_wall = 0.0; % requires velocity decomposition
INLET = [];         u_inlet = 0.0;  v_inlet = 0.0; % velocity inlet
OUTLET = []; % pressure outlet
FREESLIP = [];

%% Some default settings

if strfind(filename,'cavity')
    WALL = [INT;SBC;EBC;WBC];
    MOVINGWALL = [NBC];
    u_wall = 1.0;
elseif strfind(filename, 'cylinder')
    WALL = [INT;SBC;NBC];
    MOVINGWALL = [];        u_wall = 0.0;   v_wall = 0.0; % requires velocity decomposition
    INLET = [WBC];          u_inlet = 1.0;  v_inlet = 0.0; % velocity inlet
    OUTLET = [EBC]; % pressure outlet
    FREESLIP = [];
elseif strfind(filename, 'pitzDaily')
    WALL = [INT;NBC;SBC];
    MOVINGWALL = [];        u_wall = 0.0;   v_wall = 0.0; % requires velocity decomposition
    INLET = [WBC];          u_inlet = 10.0;  v_inlet = 0.0; % velocity inlet
    OUTLET = [EBC]; % pressure outlet
    FREESLIP = [];
end
end

