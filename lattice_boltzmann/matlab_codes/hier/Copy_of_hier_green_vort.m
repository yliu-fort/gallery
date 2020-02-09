clear; clc; close all;

% Note 01/31/2020
% staggered profile observed in fine grid which agrees with
% dicussion in paper, this staggering can be cured by apply
% averaging to coarse grid. not sure how it impacts the stability.
%% Global parameter
% Global resolution
scale = 16;

% The kinematic viscosity
csSqrInv = 3.0;
csSqr = 1.0/csSqrInv;
c = sqrt(csSqr);
%tau = 0.9; % L2err = 0.0017 (H=13) 7.1686e-05 (H = 55)
%tau = 1.0; % L2err = 0.0038 (H=13)
%tau = sqrt(3/16)+1/2; % L2err = 8.2703e-09
tau = 0.55

% Nondimensionalization
length_ui = 1.0; % m
nu_ui = 1.0/100; % m^2/sec
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

%% Hier-grid

% Level 0
level = 1; ratio = [1 1];
NX = ratio(1)*level*scale+2;NY = ratio(2)*level*scale+2;
numElements = NX*NY;
[X, Y] = ndgrid([0 linspace(0+1/2/(NX-2),ratio(1)-1/2/(NX-2),NX-2) ratio(1)],...
    [0 linspace(0+1/2/(NY-2),ratio(2)-1/2/(NY-2),NY-2) ratio(2)]);

Cl = length_ui/(NY - 2)/1.0; % boundary node take 2 * dx/2
nu = csSqr*(tau-0.5); % nu = cs^2*(tao - 0.5)
Ct =  nu/nu_ui*Cl^2; % dtstar = 1.0
NSTEPS = floor(time_length_ui / Ct);
PLOTGAP = floor(0.02 / Ct)+1;
%PLOTGAP = 1;
%NSTEPS = 64*scale*scale; % specific timestep for validation
Cu = Cl/Ct;
u_max = u_max_ui/Cu % should below 0.1 or 0.03, max < 0.4
%Ca = Cl^2/Ct;

% and the corresponding relaxation parameter
% The maximum flow speed
%Reg = 0.18; % Stability condition, max < O(10)
%tau = 0.5 + u_max/(csSqr*Reg)% tao/dt >= 0.5 + O(1/8) * u_max
Reg = u_max/(tau - 0.5)*csSqrInv
Re = Reg / Cl;

% initialise rho  ux  uy fields.
[rho, ux, uy] = taylor_green(0, NX, NY, rho0, nu, u_max);
%  initialise f[2] as equilibrium for rho  ux  uy
f = init_equilibrium(rho, ux, uy, wi, dirx, diry, NX, NY, ndir);

% interface: must sync with coarse/fine grid when streaming
gridding = zeros(NX*NY,1);

% Level 1
level = 2; ratio = [1 1];
NX1 = ratio(1)*level*scale+2;NY1 = ratio(2)*level*scale+2;
numElements1 = NX1*NY1;
[X1, Y1] = ndgrid([0 linspace(0+1/2/(NX1-2),ratio(1)-1/2/(NX1-2),NX1-2) ratio(1)],...
    [0 linspace(0+1/2/(NY1-2),ratio(2)-1/2/(NY1-2),NY1-2) ratio(2)]);

tau1 = (tau-0.5)*level + 0.5
Cl1 = length_ui/(NY1 - 2)/1.0; % boundary node take 2 * dx/2
nu1 = csSqr*(tau1-0.5); % nu = cs^2*(tao - 0.5)
Ct1 =  nu1/nu_ui*Cl1^2; % dtstar = 1.0

Cu1 = Cl1/Ct1;
u_max1 = u_max_ui/Cu1 % should below 0.1 or 0.03, max < 0.4

Reg1 = u_max1/(tau1 - 0.5)*csSqrInv
Re1 = Reg1 / Cl1;

% initialise rho  ux  uy fields.
rho0_1 = rho0/4;
[rho1, ux1, uy1] = taylor_green(0, NX1, NY1, rho0_1, nu1, u_max1);
%  initialise f[2] as equilibrium for rho  ux  uy
f1 = init_equilibrium(rho1, ux1, uy1, wi, dirx, diry, NX1, NY1, ndir);
% mapping, fine to coarse floor(uv/2)+1 [1,ni]
%coarse to fine (uv-1)*2 + 0/1 [1,ni-1]
[mcfX1,mcfY1] = ndgrid(ceil((2:2*NX-1)/2),ceil((2:2*NY-1)/2));
mapping_coarse_to_fine1 = sub2ind([NX,NY],mcfX1,mcfY1);
%gridding1 = 10.*((X1(:)-0.5).^2 + (Y1(:)-0.5).^2 > 0.3^2);
%gridding1 = 10.*((X1(:) < 0.2)|(X1(:) > 0.8)|(Y1(:) < 0.2)|(Y1(:) > 0.8));
gridding1 = 10.*((X1(:) > 0.25)&(X1(:) < 0.75)&(Y1(:) > 0.25)&(Y1(:) < 0.75));
gridding1 = 10.*((X1(:)-0.5).^2 + (Y1(:)-0.5).^2 < 0.1^2);

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
gridding = 0*gridding + 10; % all mark as fluid
gridding(sq==4) = 1000; % mark as occupied
gridding(sr>0) = 100; % mark as interface_to_fine

% Build connectivity
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

% all childs of coarse interface cell must be marked as interface
gridding1(interface_to_coarse_connectivity(gridding==100,1)) = 1;
gridding1(interface_to_coarse_connectivity(gridding==100,2)) = 1;
gridding1(interface_to_coarse_connectivity(gridding==100,3)) = 1;
gridding1(interface_to_coarse_connectivity(gridding==100,4)) = 1;

% optional: remove unused interface cells
%gridding1(interface_to_coarse_connectivity(gridding==10,1)) = 0;
%gridding1(interface_to_coarse_connectivity(gridding==10,2)) = 0;
%gridding1(interface_to_coarse_connectivity(gridding==10,3)) = 0;
%gridding1(interface_to_coarse_connectivity(gridding==10,4)) = 0;

% optional: mark f as NaN in unused cells
% to confirm no particles outside is able to propagate into the region
f(gridding==0|gridding==1000,:,:) = NaN;
f1(gridding1==0|gridding1==1000,:,:) = NaN;
ux(gridding==0|gridding==1000,:,:) = NaN;
ux1(gridding1<=1|gridding1==1000,:,:) = NaN;
uy(gridding==0|gridding==1000,:,:) = NaN;
uy1(gridding1<=1|gridding1==1000,:,:) = NaN;

% Check
%draw_hier

%% Visualization
figure;
%subplot(1,2,1)
hold on
h2 = pcolor(X);h2.XData = X;h2.YData = Y;h2.ZData = 0*X;h2.EdgeColor='none';h2.FaceColor='interp';
h21 = pcolor(X1);h21.XData = X1;h21.YData = Y1;h21.ZData = 0*X1;h21.EdgeColor='none';h21.FaceColor='interp';
colormap(desaturated_rainbow);
axis equal;axis tight;axis([0 ratio(1) 0 ratio(2)]);
%subplot(1,2,2)
%h3 = pcolor(X);h3.XData = X;h3.YData = Y;h3.ZData = 0*X;h3.EdgeColor='none';h3.FaceColor='interp';
%h31 = pcolor(X1);h31.XData = X1;h31.YData = Y1;h31.ZData = 0*X1;h31.EdgeColor='none';h31.FaceColor='interp';
%axis equal;axis tight;axis([0 ratio(1) 0 ratio(2)]);
%colormap gray
%l2 = plot(Y(1,:),Y(1,:),'-+');
fprintf('NSTEPS=%d, FRAMEGAP=%d\n',NSTEPS,PLOTGAP);
fprintf('Umax=%f\n',u_max);
fprintf('tau=%f\n',tau);
fprintf('Reg=%f\n',Reg);
fprintf('Re=%f\n',Re);

% buf to store temporary field
buf = zeros(NX*NY,ndir); % gpuArray

% buffer indicator
s.read = 1;
s.write = 2;

% mass conservation
mass_conservation_ref = sum(rho(gridding==10|gridding==100)) + sum(rho1(gridding1==10));

%%  main simulation loop; take NSTEPS time steps
convergence = Inf;
n = 0;
while((~isfinite(convergence) || convergence > 1e-15))
    
    % record statistics
    uhist = ux;
    vhist = uy;
    
    % read buffer
    buf = f(:,:,s.read);
    buf1 = f1(:,:,s.read);
    
    %fprintf("conservation: %f",sum(rho(gridding==10|gridding==100)) + sum(rho1(gridding1==10)));
    % integrate on hier grid
    % step 1: collision on both grids
    % step 2: redistribution from coarse to fine (prolongation)
    % step 3: streaming on both grids
    % step 4a:collide on fine grid
    % step 4b:streaming on both grids
    % step 5: redistribution from fine to coarse (restriction)
    % sequence
    % 1.C-collide 1. F-redis   2. F-collide 2. F-stream  
    % 2.F-collide 2. F-stream  1. C-stream  1. C-redist
    
    % collide on both grids
    rho = sum(buf,2);
    ux  = sum(bsxfun(@times,buf,dirx),2)./rho;
    uy  = sum(bsxfun(@times,buf,diry),2)./rho;
    feq=compute_equilibrium(rho, ux, uy, dirx, diry, wi);
    buf = (1.0-1/tau)*buf + (1/tau)*feq;

    % redistribution, operate->interface_cell on fine grid
    % Here's the issue, mass conservation is not achieved exactly as
    % particles from unused region propagate into the region.
    % we should only propagate particles from fluid region
    buf1(gridding1==1,:) = buf(interface_to_fine_connectivity(gridding1==1),:)/4;
    
    % F - compute rho and u
    rho1 = sum(buf1,2);
    ux1  = sum(bsxfun(@times,buf1,dirx),2)./rho1;
    uy1  = sum(bsxfun(@times,buf1,diry),2)./rho1;
    
    % F - collide
    feq1=compute_equilibrium(rho1, ux1, uy1, dirx, diry, wi);
    buf1(gridding1==10,:) = (1.0-1/tau1)*buf1(gridding1==10,:) + (1/tau1)*feq1(gridding1==10,:);
    
    % F - streaming on fine grid (swap buffer required)
    buf1 = streaming(buf1, dirx, diry, NX1, NY1, ndir);
    
    % F - compute rho and u
    rho1(gridding1==10) = sum(buf1(gridding1==10,:),2);
    ux1(gridding1==10)  = sum(bsxfun(@times,buf1(gridding1==10,:),dirx),2)./rho1(gridding1==10);
    uy1(gridding1==10)  = sum(bsxfun(@times,buf1(gridding1==10,:),diry),2)./rho1(gridding1==10);
    
    % F - collide
    feq1=compute_equilibrium(rho1, ux1, uy1, dirx, diry, wi);
    buf1(gridding1==10,:) = (1.0-1/tau1)*buf1(gridding1==10,:) + (1/tau1)*feq1(gridding1==10,:);
    
    % F - streaming on fine grid (swap buffer required)
    buf1 = streaming(buf1, dirx, diry, NX1, NY1, ndir);
    
    % streaming on coarse grid
    buf = streaming(buf, dirx, diry, NX, NY, ndir);
    
    % a cell marker is needed in fine interface cells
    % to mark in which direction the particle is updated in current level
    % rest parts will be updated from next coarser level after streaming

    isUpdated = repmat(gridding1==10,1,9);
    isUpdated1 = streaming(isUpdated, dirx, diry, NX1, NY1, ndir);
    isUpdated2 = streaming(isUpdated1, dirx, diry, NX1, NY1, ndir);
    isUpdated = (isUpdated1|isUpdated2)&(~isUpdated);
    
    for i = 2:9
    buf1((gridding1==1)&(~isUpdated(:,i)),i) =...
        buf(interface_to_fine_connectivity((gridding1==1)&(~isUpdated(:,i))),i)/4;
    end
    
    % redistribution, operate->interface_cell on coarse grid
    rbuf = buf1(interface_to_coarse_connectivity(gridding==100,1),:)  ...
        + buf1(interface_to_coarse_connectivity(gridding==100,2),:)  ...
        + buf1(interface_to_coarse_connectivity(gridding==100,3),:)  ...
        + buf1(interface_to_coarse_connectivity(gridding==100,4),:) ;
    
    buf(gridding==100,:) = rbuf;
    
    % Output
    f(:,:,s.write) = buf;
    f1(:,:,s.write) = buf1;
    
    % swap buffer
    s.read = mod(s.read,2)+1; % self-increment
    s.write = mod(s.read,2)+1; % offset
    
    % Advance
    n = n+1;
    %if(n == 1) break;end
    % visualization
    if(mod(n,PLOTGAP) == 0 || n == NSTEPS)
        % Convergence
        conv = sum((ux - uhist).*(ux - uhist) + (uy - vhist).*(uy - vhist),'omitnan');
        convergence = sqrt(conv/sum(uhist.^2 + vhist.^2,'omitnan'));

        % conservation
        mass_conservation = sum(rho(gridding==10|gridding==100)) + sum(rho1(gridding1==10));
        fprintf('res: %e | mass conservation: %.15f>\n', convergence, mass_conservation/mass_conservation_ref);
        
        % Visualization
        rh = reshape((rho),NX,NY);
        uh = reshape((Cu*ux),NX,NY);
        vh = reshape((Cu*uy),NX,NY);
        h2.CData = sqrt(uh.^2 + vh.^2);
        
        rh1 = reshape((rho1),NX1,NY1);
        uh1 = reshape((Cu1*ux1),NX1,NY1);
        vh1 = reshape((Cu1*uy1),NX1,NY1);
        h21.CData = sqrt(uh1.^2 + vh1.^2);
        
        %[dvdx,~] = gradient(vh);[dudy,~] = gradient(uh);
        %h3.CData = dvdx - dudy;
        %h3.CData = max(abs(dudy),abs(dvdx));
        %l2.YData = uh(floor(NX/2),:);
        %h2.CData = vh;
        %h2.CData = rh-rho0;
        %subplot(1,2,1)
        %axis equal;axis tight;axis([0 ratio(1) 0 ratio(2)]);
        %subplot(1,2,2)
        %axis equal;axis tight;axis([0 ratio(1) 0 ratio(2)]);
        %title(sprintf('%d',n/PLOTGAP))%,caxis([0 2])
        drawnow;
        
    end
end



