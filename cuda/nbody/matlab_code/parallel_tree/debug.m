clear; clc; close all;

%%
rng(77777)
softenFactor = 1e-3;
graviFactor = 1;
dt = 1e-3;

n = 128;
s = 1;
p.x = s*rand(n,1);
p.y = s*rand(n,1);
p.z = s*rand(n,1);
p.w = ones(n,1);

rp = 33;
p.x(rp)=p.x(rp-1);
p.y(rp)=p.y(rp-1);
p.z(rp)=p.z(rp-1);

rp = 36;
p.x(rp)=p.x(rp+1);
p.y(rp)=p.y(rp+1);
p.z(rp)=p.z(rp+1);

p.vx = zeros(n,1);
p.vy = zeros(n,1);
p.vz = zeros(n,1);

p.ax = zeros(n,1);
p.ay = zeros(n,1);
p.az = zeros(n,1);

tree = {};

%% Calc morton code
morton = uint32(morton3D(p.x,p.y,p.z));

%% Sort (radix will be best)
[morton,sortedID] = sort(morton);

%% Gen hierarchy
[tree.connectivity,tree.bbox,tree.segment,tree.isoctreenode, tree.level] = ...
    generateHierachy(morton);

%% Calc treeMass
tree.mass = zeros(n-1,4);
for i = 1:n-1
    if(tree.isoctreenode(i))
        for j = (tree.segment(i,1)+1):(tree.segment(i,2)+1)
            idx = sortedID(j);
            tree.mass(i,1) = tree.mass(i,1) + p.x(idx)*p.w(idx);
            tree.mass(i,2) = tree.mass(i,2) + p.y(idx)*p.w(idx);
            tree.mass(i,3) = tree.mass(i,3) + p.z(idx)*p.w(idx);
            tree.mass(i,4) = tree.mass(i,4) +          p.w(idx);
        end
        if(tree.mass(i,4) ~= 0)
            tree.mass(i,1) = tree.mass(i,1)/tree.mass(i,4);
            tree.mass(i,2) = tree.mass(i,2)/tree.mass(i,4);
            tree.mass(i,3) = tree.mass(i,3)/tree.mass(i,4);
        end
    end
end

%% Calc accel

[acc.x,acc.y,acc.z] = traverse(uint32(sortedID-1),p.x,p.y,p.z,p.w,...
    tree.connectivity,tree.mass, tree.bbox, tree.isoctreenode, ...
    tree.segment, tree.level, softenFactor, graviFactor);

%% Update velocity and position

p.ax = acc.x;
p.ay = acc.y;
p.az = acc.z;

p.vx = p.vx + p.ax*dt;
p.vy = p.vy + p.ay*dt;
p.vz = p.vz + p.az*dt;

p.x = p.x + p.vx*dt;
p.y = p.y + p.vy*dt;
p.z = p.z + p.vz*dt;

%% Resize bounding box & calc gravityFactor

p.x = p.x - min(p.x);
p.y = p.y - min(p.y);
p.z = p.z - min(p.z);
dist = max(max(abs(max(p.x)), abs(max(p.y))),abs(max(p.z)));
if(dist == 0), disp("Scale no longer works.");end

% scale spatial resolution
p.x = p.x/dist;
p.y = p.y/dist;
p.z = p.z/dist;
p.vx = p.vx/dist;
p.vy = p.vy/dist;
p.vz = p.vz/dist;
graviFactor = graviFactor/dist^2;

% scale temporal resolution
dt = 1e-2/(1e1 + sqrt(max(p.vx.^2 + p.vy.^2 + p.vz.^2)));

%% Some plots

plot(tree.segment(tree.isoctreenode > 0,:).',...
    [tree.level(tree.isoctreenode > 0)/3 tree.level(tree.isoctreenode > 0)/3].'...
    ,'.-')

% draw tree connectivity
plot(tree.connectivity, [tree.level tree.level],'.-')

