clear; clc; close all;

%%
rng(55555)

nl = 20; % 10:1024 13:8192 16:65536

softenFactor = 1e-4;
graviFactor = (1.5e-1*(2^((16-nl)/3)))^2;
dt = 1e-6;

n = 2^nl;
p.x = rand(n,1);
p.y = rand(n,1);
p.z = rand(n,1);
p.w = ones(n,1);
p.w(1:64:end) = -63;

t.x = p.x - 0.5;
t.y = p.y - 0.5;
t.z = p.z - 0.5;
t.w = sqrt(1e-6 + t.x.^2 + t.y.^2 + t.z.^2);

expand = 225;
p.vx = zeros(n,1) + expand*t.x;
p.vy = zeros(n,1) + expand*t.y;
p.vz = zeros(n,1) + expand*t.z;

p.ax = zeros(n,1);
p.ay = zeros(n,1);
p.az = zeros(n,1);

tree = {};

t = 0;
t_final = 1e3;

offset.x = 0;offset.y = 0; offset.z = 0;
fig = scatter3(p.x, p.y, p.z, 10, '.','k');axis equal
th=annotation('textbox',[.05 .05 .2 .1],'String',[],'FitBoxToText','on');

%%
iter = 0;
while t < t_final
    tic
    t = t + dt;
    
    %% Calc morton code
    
    morton = uint32(morton3D(p.x,p.y,p.z));
    
    %% Sort (radix will be best)
    
    [morton,sortedID] = sort(morton);
    
    %% Gen hierarchy
    
    [tree.connectivity,tree.bbox,tree.segment,tree.isoctreenode, tree.level]...
        = generateHierachy(morton);
    
    %% Calc treeMass
    
    tree.mass = calc_weights(sortedID,p.x,p.y,p.z,p.w,...
    tree.isoctreenode,tree.segment);

    
    %% Calc accel
    
    [acc.x,acc.y,acc.z] = traverse(uint32(sortedID-1),p.x,p.y,p.z,p.w,...
        tree.connectivity,tree.mass, tree.bbox, tree.isoctreenode, ...
        tree.segment, tree.level, softenFactor, graviFactor);
    
    %% Update velocity and position
    
    p.vx = p.vx + 0.5*p.ax*dt;
    p.vy = p.vy + 0.5*p.ay*dt;
    p.vz = p.vz + 0.5*p.az*dt;
    
    p.x = p.x + p.vx*dt;
    p.y = p.y + p.vy*dt;
    p.z = p.z + p.vz*dt;
    
    p.ax = acc.x;
    p.ay = acc.y;
    p.az = acc.z;
    
    p.vx = p.vx + 0.5*p.ax*dt;
    p.vy = p.vy + 0.5*p.ay*dt;
    p.vz = p.vz + 0.5*p.az*dt;
    
    
    %% Resize bounding box & calc gravityFactor
    
    offset.x = offset.x + min(p.x)/sqrt(graviFactor);
    offset.y = offset.y + min(p.y)/sqrt(graviFactor);
    offset.z = offset.z + min(p.z)/sqrt(graviFactor);
    
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
    p.ax = p.ax/dist;
    p.ay = p.ay/dist;
    p.az = p.az/dist;
    
    graviFactor = graviFactor/dist^2;
    
    % scale temporal resolution
    dt = 4e-2/(4e-3 + sqrt(median(p.vx.^2 + p.vy.^2 + p.vz.^2)));
    
    
    %% Plot
    fig.XData = p.x;
    fig.YData = p.y;
    fig.ZData = p.z;
    %fig.SizeData = 25*abs(p.w).^(1/3);
    %fig.SizeData = 40;
    %fig.CData = sqrt(p.ax.^2 + p.ay.^2 + p.az.^2)*dt^2;
    %axis equal, axis([0 1 0 1 0 1]/sqrt(graviFactor) + [offset.x offset.x offset.y offset.y offset.z offset.z])
    th.String = ['Scale = ' num2str(1.0/sqrt(graviFactor)) ', Time = ' num2str(t) ', dt = ' num2str(dt) ', nobj = ' num2str(n)];
    drawnow;
    
    %% Report
    iter = iter + 1;
    %if(iter == 109) break;end
    %data.time(i) = t;
    %data.clustering(i) = var(tree.mass(:,4)./(tree.bbox(:,4).^3),'omitnan');
    disp(["iter = " num2str(iter)])
    toc
end

% Unused code
%{
%% Merge&sort (optional)
    
    [morton,~,oIdx] = unique(morton, 'stable');
    
    p.x = p.x.*p.w;
    p.y = p.y.*p.w;
    p.z = p.z.*p.w;
    
    p.vx = p.vx.*p.w;
    p.vy = p.vy.*p.w;
    p.vz = p.vz.*p.w;
    
    p.ax = p.ax.*p.w;
    p.ay = p.ay.*p.w;
    p.az = p.az.*p.w;
    
    p.w = accumarray(oIdx,p.w);
    p.x = accumarray(oIdx,p.x)./p.w;
    p.y = accumarray(oIdx,p.y)./p.w;
    p.z = accumarray(oIdx,p.z)./p.w;
    
    p.vx = accumarray(oIdx,p.vx)./p.w;
    p.vy = accumarray(oIdx,p.vy)./p.w;
    p.vz = accumarray(oIdx,p.vz)./p.w;
    
    p.ax = accumarray(oIdx,p.ax)./p.w;
    p.ay = accumarray(oIdx,p.ay)./p.w;
    p.az = accumarray(oIdx,p.az)./p.w;
    
    % Remove point with zero mass
    p.x = p.x(p.w~=0);p.y = p.y(p.w~=0);p.z = p.z(p.w~=0);
    p.vx = p.vx(p.w~=0);p.vy = p.vy(p.w~=0);p.vz = p.vz(p.w~=0);
    p.ax = p.ax(p.w~=0);p.ay = p.ay(p.w~=0);p.az = p.az(p.w~=0);
    
    morton = morton(p.w~=0);%sortedID = sortedID(p.w~=0);
    p.w = p.w(p.w~=0);
    
    n = numel(morton);

%}