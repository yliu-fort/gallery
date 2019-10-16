%% Incompressible flow Solver

clear; clc; close all;clear mex;
addpath ./lib ./solver ./lib/readmesh ./numeric
addpath ./lib/mexfunction

%% Import geometry information

filename =  'cavity';
eval('polymesh_remake_test_script')
% some demos: cavity, cylinder

%% Preprocessing

PhysicalProperties = Preprocessing(filename, GeometryProperties );

%% Iteration settings

IterSettings = IterSetting('default');

%% Start Computing

%mex ./lib/mexfunction/interpPhi.c -outdir ./lib/mexfunction
icoSolver(IterSettings,PhysicalProperties,GeometryProperties);

%% Post_processing

%eval(Postprocessing)
