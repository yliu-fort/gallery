clear;clc;close all;

% Read the JSON file
InArgs = jsondecode(fileread('wave_profile.json'));
globalArgs = jsondecode(fileread('config.json'));

% Specify the file path
[tb, header_string, variable_names, meshsize] = read_dumpfile("db.dat");
%tb1 = read_dumpfile("db1.dat");

%%
dx = 0.0001;
dy = 0.0001;
mask = tb.gpy < max(tb.gpy(tb.f>0.001))*4 & tb.isleaf > 0;

pot_S22_px = 0*tb.gpx;
pot_S22_mx = 0*tb.gpx;
pot_S22_py = 0*tb.gpy;
pot_S22_my = 0*tb.gpy;

pot_S22_px(mask) = generate_potential_fun(InArgs.bsk.A,InArgs.bsk.Tw, ...
    InArgs.bsk.h,InArgs.bsk.rev_period,tb.gpx(mask).'+dx,tb.gpy(mask).');

pot_S22_mx(mask) = generate_potential_fun(InArgs.bsk.A,InArgs.bsk.Tw, ...
    InArgs.bsk.h,InArgs.bsk.rev_period,tb.gpx(mask).'-dx,tb.gpy(mask).');

pot_S22_py(mask) = generate_potential_fun(InArgs.bsk.A,InArgs.bsk.Tw, ...
    InArgs.bsk.h,InArgs.bsk.rev_period,tb.gpx(mask).',tb.gpy(mask).'+dy);

pot_S22_my(mask) = generate_potential_fun(InArgs.bsk.A,InArgs.bsk.Tw, ...
    InArgs.bsk.h,InArgs.bsk.rev_period,tb.gpx(mask).',tb.gpy(mask).'-dy);

ux = (pot_S22_px - pot_S22_mx)/2/dx;
uy = (pot_S22_py - pot_S22_my)/2/dy;


%tb.ux = ux.*(tb.f>0.001);
%tb.uy = uy.*(tb.f>0.001);
tb.ux(tb.isleaf==0) = nan;
tb.uy(tb.isleaf==0) = nan;

clf, figure
scatter(tb.gpx, tb.gpy, [], tb.uy, '.')
colorbar


%% Read the header string from line 1-3
file_id = fopen('db2.dat', 'w');
fprintf(file_id, "%s\n", header_string{1}{1});
fprintf(file_id, "%s\n", header_string{1}{2});
fprintf(file_id, "%s", strjoin(variable_names));
writetable(tb, 'db2.dat','Delimiter',' ','WriteMode','append','WriteVariableNames',false);
fclose(file_id);

%%
function [tb, header_string, variable_names, meshsize] = read_dumpfile(file_path)
% Specify the file path
%file_path = 'db.dat';

% Open the file for reading
file_id = fopen(file_path, 'r');

% Read the header string from line 1
header_string = textscan(file_id, '%s', 2, 'delimiter', '\n', 'headerlines', 0);

% Read the variable names from line 3
variable_names = textscan(file_id, '%s', 1, 'delimiter', '\n', 'headerlines', 0);
variable_names = strsplit(variable_names{1}{1}, ' ');
meshsize = variable_names((end-3):end);
meshsize = arrayfun(@str2num, string(meshsize));

% Close the file
fclose(file_id);

% Read the table data skipping the first three lines
tb = readtable(file_path, 'Delimiter', ' ', 'HeaderLines', 3, 'ReadVariableNames', false);

% Assign the variable names to the table
tb.Properties.VariableNames = replace(variable_names(1:end-4),'.','');

end