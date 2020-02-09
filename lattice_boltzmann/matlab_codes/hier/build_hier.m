%% Build hier data structure
clear; clc; close all;

% Note 01/31/2020
% staggered profile observed in fine grid which agrees with
% dicussion in paper, this staggering can be cured by apply
% averaging to coarse grid. not sure how it impacts the stability.
%% Global parameter
% Global resolution
scale = 64;ratio = [2 1];
refine_level = 3;

% The kinematic viscosity
csSqrInv = 3.0;
csSqr = 1.0/csSqrInv;
c = sqrt(csSqr);
%tau = 0.9; % L2err = 0.0017 (H=13) 7.1686e-05 (H = 55)
%tau = 1.0; % L2err = 0.0038 (H=13)
%tau = sqrt(3/16)+1/2; % L2err = 8.2703e-09
tau = 0.502;

% Nondimensionalization
length_ui = 1.0; % m
nu_ui = 1.0/10000; % m^2/sec
time_length_ui = 7200*0.040; % sec
u_max_ui = 1.0;

% D2Q9
ndir = 9;
% The lattice weights
w0 = 4.0/9.0; % zero weight
ws = 1.0/9.0; % adjacent weight
wd = 1.0/36.0; % diagonal weight
% Arrays of the lattice weights and direction components
wi = ([w0 ws ws ws ws wd wd wd wd]); % gpuArray
dirx = ([0 1 0 -1  0 1 -1 -1  1]);
diry = ([0 0 1  0 -1 1  1 -1 -1]);
bi   = ([1 4 5  2  3 8  9  6  7]);

% The fluid density
rho0 = 1.0;

% initial functor
%init = @taylor_green;
init = @wake;

% hier
hier = {};
%% Level i
for i = 1:refine_level
    data.level = 2^(i-1);
    data.ratio = ratio;
    data.NX = data.ratio(1)*data.level*scale+2;data.NY = data.ratio(2)*data.level*scale+2;
    data.numElements = data.NX*data.NY;
    [data.X, data.Y] = ndgrid(...
        data.ratio(1)*[0 linspace(0+1/2/(data.NX-2),1-1/2/(data.NX-2),data.NX-2) 1],...
        data.ratio(2)*[0 linspace(0+1/2/(data.NY-2),1-1/2/(data.NY-2),data.NY-2) 1]);
    
    data.Cl = length_ui/(data.NY - 2)/1.0; % boundary node take 2 * dx/2
    data.tau = (tau-0.5)*data.level + 0.5;
    data.nu = csSqr*(data.tau-0.5); % nu = cs^2*(tao - 0.5)
    data.Ct =  data.nu/nu_ui*data.Cl^2; % dtstar = 1.0
    data.NSTEPS = floor(time_length_ui / data.Ct);
    data.PLOTGAP = floor(0.02 / data.Ct)+1;
    %PLOTGAP = 1;
    %NSTEPS = 64*scale*scale; % specific timestep for validation
    data.Cu = data.Cl/data.Ct;
    data.u_max = u_max_ui/data.Cu; % should below 0.1 or 0.03, max < 0.4
    %Ca = Cl^2/Ct;
    
    % and the corresponding relaxation parameter
    % The maximum flow speed
    %Reg = 0.18; % Stability condition, max < O(10)
    %tau = 0.5 + u_max/(csSqr*Reg)% tao/dt >= 0.5 + O(1/8) * u_max
    data.Reg = data.u_max/(data.tau - 0.5)*csSqrInv;
    data.Re = data.Reg / data.Cl;
    data.rho0 = rho0/(4^(i-1));
    
    % initialise rho  ux  uy fields.
    [data.rho, data.ux, data.uy] = feval(init, 0, data.NX, data.NY, data.rho0, data.nu, data.u_max);
    %  initialise f[2] as equilibrium for rho  ux  uy
    %data.f = init_equilibrium(data.rho, data.ux, data.uy, wi, dirx, diry, data.NX, data.NY, ndir);
    data.f=compute_equilibrium(data.rho, data.ux, data.uy, dirx, diry, wi);
    
    % interface: must sync with coarse/fine grid when streaming
    data.gridding = zeros(data.NX*data.NY,1);
    
    % visualization
    xld = data.X-1/2/(data.NX-2)*data.ratio(1);
    yld = data.Y-1/2/(data.NY-2)*data.ratio(2);

    xrd = data.X+1/2/(data.NX-2)*data.ratio(1);
    yrd = data.Y-1/2/(data.NY-2)*data.ratio(2);
    
    xrt = data.X+1/2/(data.NX-2)*data.ratio(1);
    yrt = data.Y+1/2/(data.NY-2)*data.ratio(2);
    
    xlt = data.X-1/2/(data.NX-2)*data.ratio(1);
    ylt = data.Y+1/2/(data.NY-2)*data.ratio(2);

    xld(end,:) = xld(end,:)+1/2/(data.NX-2)*data.ratio(1);
    xlt(end,:) = xlt(end,:)+1/2/(data.NX-2)*data.ratio(1);
    yld(:,end) = yld(:,end)+1/2/(data.NY-2)*data.ratio(2);
    yrd(:,end) = yrd(:,end)+1/2/(data.NY-2)*data.ratio(2);
    xrd(1,:) = xrd(1,:)-1/2/(data.NX-2)*data.ratio(1);
    xrt(1,:) = xrt(1,:)-1/2/(data.NX-2)*data.ratio(1);
    yrt(:,1) = yrt(:,1)-1/2/(data.NY-2)*data.ratio(2);
    ylt(:,1) = ylt(:,1)-1/2/(data.NY-2)*data.ratio(2);
    
    data.pX = [xld(:)';  xlt(:)'; xrt(:)'; xrd(:)'];
    data.pY = [yld(:)';  ylt(:)'; yrt(:)'; yrd(:)'];
    
    % push to hier
    hier{i} = data;
    data = [];
end

%% local refinement
%rlvl=2;
%hier{rlvl}.gridding = 10.*double(((hier{rlvl}.X(:) > 0.1)&(hier{rlvl}.X(:) < 0.9)&(hier{rlvl}.Y(:) > 0.1)&(hier{rlvl}.Y(:) < 0.9)));
%rlvl=refine_level;
%hier{rlvl}.gridding = 10.*double(((hier{rlvl}.X(:)-0.5).^2 + (hier{rlvl}.Y(:)-0.5).^2 < 0.25^2));

%rlvl=refine_level;
%hier{rlvl}.gridding = 10.*double(((hier{rlvl}.X(:) < 0.15)|(hier{rlvl}.X(:) > 0.85)|(hier{rlvl}.Y(:) < 0.15)|(hier{rlvl}.Y(:) > 0.85)));
rlvl=2;
hier{rlvl}.gridding = 10.*double(((hier{rlvl}.X(:)-0.5).^2 + (hier{rlvl}.Y(:)-0.5).^2 < 0.2^2));

rlvl=refine_level;
hier{rlvl}.gridding = 10.*double(((hier{rlvl}.Y(:) >= 0.2)&(hier{rlvl}.Y(:) < 0.3))|((hier{rlvl}.Y(:) >= 0.7)&(hier{rlvl}.Y(:) < 0.8)));
%rlvl=2;
%hier{rlvl}.gridding = 10.*double((hier{rlvl}.Y(:) > 0.15)&(hier{rlvl}.Y(:) < 0.85));


% mark coarest level as fluid cell
hier{1}.gridding = 10 + hier{1}.gridding;

%% Build connectivity
for i = 1:refine_level-1
    % level 1,2,3 -> to_fine
    hier{i}.interface_to_fine_connectivity = [];
    [sx,sy] = ndgrid([1 2:hier{i+1}.NX-1 hier{i+1}.NX],[1 2:hier{i+1}.NY-1 hier{i+1}.NY]);
    sx = sx(:);sy = sy(:);
    hier{i}.interface_to_fine_connectivity=...
        sub2ind([hier{i}.NX,hier{i}.NY],floor(sx/2)+1,floor(sy/2)+1);
    
    % level 2,3,4 -> to_coarse
    hier{i+1}.interface_to_coarse_connectivity = [];
    subindex = [0 0;1 0;0 1;1 1];
    for ii = 1:4
        [sx,sy] = ndgrid([1 (2+subindex(ii,1)):2:(hier{i+1}.NX-2+subindex(ii,1)) hier{i+1}.NX],...
            [1 (2+subindex(ii,2)):2:(hier{i+1}.NY-2+subindex(ii,2)) hier{i+1}.NY]);
        sx = sx(:);sy = sy(:);
        hier{i+1}.interface_to_coarse_connectivity = [hier{i+1}.interface_to_coarse_connectivity ...
            sub2ind([hier{i+1}.NX,hier{i+1}.NY],sx,sy)];
    end
    
end

%% mark cells
%0: unused 1: to_coarse 10: fluid 100: to_fine 1000:occupied
% alternative approach: bool isOccupied, bool isUpdate, bool isInterface
% from finest mesh
for i = refine_level:-1:1
    % pass 1
    if(i < refine_level)
        % pass 2
        % sweep subcells
        subcellCount = 0*hier{i}.X;
        subindex = [0 0;1 0;0 1;1 1];
        for ii = 1:4
            [sx,sy] = ndgrid([1 (2+subindex(ii,1)):2:(hier{i+1}.NX-2+subindex(ii,1)) hier{i+1}.NX],...
                [1 (2+subindex(ii,2)):2:(hier{i+1}.NY-2+subindex(ii,2)) hier{i+1}.NY]);
            subcellCount = subcellCount + hier{i+1}.gridding(sub2ind([hier{i+1}.NX,hier{i+1}.NY],sx,sy));
        end
        so = fix(subcellCount/1000); % occupied
        sp = fix((subcellCount - so*1000)/100); % to_fine
        sq = fix((subcellCount - so*1000 - sp*100)/10); % fluid
        sr = subcellCount - so*1000 - sp*100 - sq*10; % to_coarse
        st = 8-so-sp-sq-sr; %unused
        % may have some extra function to determine the local refinement
        %gridding = 0*gridding + 10; % all mark as fluid
        hier{i}.gridding(((so+sp+sq) == 4)) = 1000; % mark as occupied
        hier{i}.gridding(sr>0) = 100; % mark as interface_to_fine
        
    end
    % pass 2
    if(i > 1)
        % sweep neighbours
        % add fluid cell surround to_fine interface cell to avoid conflict
        cellCount = 0*hier{i}.X;
        for ii = 2:9
            [sx,sy] = ndgrid([1 (2+dirx(ii)):(hier{i}.NX-1+dirx(ii)) hier{i}.NX],...
                [1 (2+diry(ii)):(hier{i}.NY-1+diry(ii)) hier{i}.NY]);
            cellCount = cellCount + hier{i}.gridding(sub2ind([hier{i}.NX,hier{i}.NY],sx,sy));
        end
        cellCount = cellCount(:);
        o = fix(cellCount/1000); % occupied
        p = fix((cellCount - o*1000)/100); % to_fine
        hier{i}.gridding((p>0)&(hier{i}.gridding==0)) = 10;
        
        % sweep neighbours
        % mark to_coarse interface cells
        cellCount = 0*hier{i}.X;
        for ii = 2:9
            [sx,sy] = ndgrid([1 (2+dirx(ii)):(hier{i}.NX-1+dirx(ii)) hier{i}.NX],...
                [1 (2+diry(ii)):(hier{i}.NY-1+diry(ii)) hier{i}.NY]);
            cellCount = cellCount + hier{i}.gridding(sub2ind([hier{i}.NX,hier{i}.NY],sx,sy));
        end
        cellCount = cellCount(:);
        o = fix(cellCount/1000); % occupied
        p = fix((cellCount - o*1000)/100); % to_fine
        q = fix((cellCount - o*1000 - p*100)/10); % fluid
        r = cellCount - o*1000 - p*100 - q*10; % to_coarse
        t = 8-o-p-q-r; %unused
        flag = (o>0|p>0|q>0)&(t>0)&(hier{i}.gridding<=1);
        hier{i}.gridding(flag) = 1;
        % both f/u neighbour and self-not fluid cell -> to_coarse
        
    end
    % pass 3
    if(i < refine_level)
        % all childs of coarse interface cell must be marked as to_coarse interface
        hier{i+1}.gridding(hier{i+1}.interface_to_coarse_connectivity(hier{i}.gridding==100,1)) = 1;
        hier{i+1}.gridding(hier{i+1}.interface_to_coarse_connectivity(hier{i}.gridding==100,2)) = 1;
        hier{i+1}.gridding(hier{i+1}.interface_to_coarse_connectivity(hier{i}.gridding==100,3)) = 1;
        hier{i+1}.gridding(hier{i+1}.interface_to_coarse_connectivity(hier{i}.gridding==100,4)) = 1;
        
    end
    
end

%% Replace unused values to NaN -> detect false numerical propagation
for i = 1:refine_level
    hier{i}.f(hier{i}.gridding==0|hier{i}.gridding==1|hier{i}.gridding==1000,:) = NaN;
    hier{i}.ux(hier{i}.gridding==0|hier{i}.gridding==1|hier{i}.gridding==1000) = NaN;
    hier{i}.uy(hier{i}.gridding==0|hier{i}.gridding==1|hier{i}.gridding==1000) = NaN;
end

%% Draw grid
%draw_hier;

%% cell marker
% mark interface cells
% Hier-info -> 0: unused 1: to_coarse 10: fluid 100: to_fine 1000:occupied
% 1000:occupied by coarse
% sum%100/10 = p, sum%10 == q -> p interface + q fluid
% sum/100 = r -> r interface_to_fine
% t = 8-p-q-r -> t unused cell
% sweep 1: finest mesh, check surrounding cell marks, if fluid|to_coarse + unused
% mark, mark as interface to_coarse cell
% sweep 2: coarse mesh, check next finer mesh, if all fluid, mark as unused cell
% if all unused, mark as fluid cell
% if interface presented, mark as interface to_fine
% then perform sweep 1 in current level to mark to_coarse cell
% sweep 3: coarsest mesh, only perform sweep 2

% pass 1
% sweep neighbours
%{
cellCount = 0*X1;
for i = 2:9
    [sx,sy] = ndgrid([1 (2+dirx(i)):(NX1-1+dirx(i)) NX1],...
        [1 (2+diry(i)):(NY1-1+diry(i)) NY1]);
    cellCount = cellCount + gridding1(sub2ind([NX1,NY1],sx,sy));
end
cellCount = cellCount(:);
p = floor(cellCount/100); % to_fine
q = floor((cellCount - p*100)/10); % fluid
r = cellCount - p*100 - q*10; % to_coarse
t = 8-p-q-r; %unused
gridding1(((q>0)|(p>0))&(t>0)&(gridding1~=10)) = 1;
% both f/u neighbour and self-not fluid cell -> to_coarse

% pass 2
% sweep subcells
subcellCount = 0*X;
subindex = [0 0;1 0;0 1;1 1];
for i = 1:4
    [sx,sy] = ndgrid([1 (2+subindex(i,1)):2:(NX1-2+subindex(i,1)) NX1],...
        [1 (2+subindex(i,2)):2:(NY1-2+subindex(i,2)) NY1]);
    subcellCount = subcellCount + gridding1(sub2ind([NX1,NY1],sx,sy));
end
sp = floor(subcellCount/100); % to_fine
sq = floor((subcellCount - sp*100)/10); % fluid
sr = subcellCount - sp*100 - sq*10; % to_coarse
st = 4-sp-sq-sr; %unused
% may have some extra function to determine the local refinement
%gridding = 0*gridding + 10; % all mark as fluid
gridding(sq==4) = 1000; % mark as occupied
gridding(sr>0) = 100; % mark as interface_to_fine
%}
% Build connectivity
%{
[sx,sy] = ndgrid([1 2:NX1-1 NX1],[1 2:NY1-1 NY1]);
sx = sx(:);sy = sy(:);
interface_to_fine_connectivity=...
    sub2ind([NX,NY],floor(sx/2)+1,floor(sy/2)+1);

interface_to_coarse_connectivity = [];
subindex = [0 0;1 0;0 1;1 1];
for i = 1:4
    [sx,sy] = ndgrid([1 (2+subindex(i,1)):2:(NX1-2+subindex(i,1)) NX1],...
        [1 (2+subindex(i,2)):2:(NY1-2+subindex(i,2)) NY1]);
    sx = sx(:);sy = sy(:);
    interface_to_coarse_connectivity = [interface_to_coarse_connectivity ...
        sub2ind([NX1,NY1],sx,sy)];
end
%}

% all childs of coarse interface cell must be marked as interface
%gridding1(interface_to_coarse_connectivity(gridding==100,1)) = 1;
%gridding1(interface_to_coarse_connectivity(gridding==100,2)) = 1;
%gridding1(interface_to_coarse_connectivity(gridding==100,3)) = 1;
%gridding1(interface_to_coarse_connectivity(gridding==100,4)) = 1;

% optional: remove unused interface cells
%gridding1(interface_to_coarse_connectivity(gridding==10,1)) = 0;
%gridding1(interface_to_coarse_connectivity(gridding==10,2)) = 0;
%gridding1(interface_to_coarse_connectivity(gridding==10,3)) = 0;
%gridding1(interface_to_coarse_connectivity(gridding==10,4)) = 0;

% optional: mark f as NaN in unused cells
% to confirm no particles outside is able to propagate into the region
%f(gridding==0|gridding==1000,:,:) = NaN;
%f1(gridding1==0|gridding1==1000,:,:) = NaN;
%ux(gridding==0|gridding==1000,:,:) = NaN;
%ux1(gridding1==0|gridding1==1000,:,:) = NaN;
%uy(gridding==0|gridding==1000,:,:) = NaN;
%uy1(gridding1==0|gridding1==1000,:,:) = NaN;