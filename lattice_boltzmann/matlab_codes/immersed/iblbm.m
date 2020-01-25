clear; clc; close all;

%%
% Compile ptx code first
%system('nvcc -ptx d2q9.cu')
% if using atomicAdd in kernel function processComputeIbNode then
% sm_60, sm_61, sm_70 flag required
system('nvcc -ptx -arch=sm_60 d2q9.cu')

disp("Unknown issue when using scale >= 8 or large NX,NY, device crashed")
% guess it is related to unknown memory management issue between multiple
% kernel functions(PTXs), should be a bug

% Reset device
reset(gpuDevice(1));

% Initialization
scale = 4; ratio = [4 1];
NX = ratio(1)*32*scale+2;NY = ratio(2)*32*scale+2;
%NX = 5;NY = 11;
numElements = NX*NY;
[X, Y] = ndgrid(ratio(1)*[0 linspace(0+1/2/(NX-2),1-1/2/(NX-2),NX-2) 1],...
    ratio(2)*[0 linspace(0+1/2/(NY-2),1-1/2/(NY-2),NY-2) 1]);

% The kinematic viscosity
csSqrInv = 3.0;
csSqr = 1.0/csSqrInv;
c = sqrt(csSqr);
%tau = 0.9; % L2err = 0.0017 (H=13) 7.1686e-05 (H = 55)
%tau = 1.0; % L2err = 0.0038 (H=13)
%tau = sqrt(3/16)+1/2; % L2err = 8.2703e-09
tau = 0.55;

% Nondimensionalization
length_ui = 1.0; % m
nu_ui = 1.0/1000; % m^2/sec
time_length_ui = 7200*0.040; % sec
u_max_ui = 1.0;

Cl = length_ui/(NY - 2)/1.0; % boundary node take 2 * dx/2
nu = csSqr*(tau-0.5); % nu = cs^2*(tao - 0.5)
Ct =  nu/nu_ui*Cl^2; % dtstar = 1.0
NSTEPS = floor(time_length_ui / Ct);
PLOTGAP = floor(0.02 / Ct)+1;
%PLOTGAP = 1;
%NSTEPS = 64*scale*scale; % specific timestep for validation
Cu = Cl/Ct;
u_max = u_max_ui/Cu; % should below 0.1 or 0.03, max < 0.4
%Ca = Cl^2/Ct;

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

% and the corresponding relaxation parameter
% The maximum flow speed
%Reg = 0.18; % Stability condition, max < O(10)
%tau = 0.5 + u_max/(csSqr*Reg)% tao/dt >= 0.5 + O(1/8) * u_max
Reg = u_max/(tau - 0.5)*csSqrInv;
Re = Reg / Cl;
% The fluid density
rho0 = 1.0;
uwall = 0;
uc = u_max;
dp = 8*nu*uc*Cl^2*(NX-1);
drho = -dp*csSqrInv;

Fx = 0;
%Fx = -8*nu*uc*Cl^2;
Fy = 0;

% Boundary
north=false(NX,NY);
north(:,end) = true;
north = north(:);% gpuArray

south=false(NX,NY);
south(:,1) = true;
south = south(:);% gpuArray

north2=false(NX,NY);
north2(:,end-1) = true;
north2 = north2(:);% gpuArray

south2=false(NX,NY);
south2(:,2) = true;
south2 = south2(:);% gpuArray

west=false(NX,NY);
west(1,:) = true;
west = west(:);

west2=false(NX,NY);
west2(2,:) = true;
west2 = west2(:);

east=false(NX,NY);
east(end,:) = true;
east = east(:);

east2=false(NX,NY);
east2(end-1,:) = true;
east2 = east2(:);

mask =  north | south | east | west;

% periodic pressure drop
%east = east & ~(north|south);east2 = east2 & ~(north|south);
%west = west & ~(north|south);west2 = west2 & ~(north|south);


%cylinder=((X-center(1)).^2 + (Y-center(2)).^2) < radius^2;
%cylinder = cylinder(:);
%mask = mask | cylinder;

% Immersed Boundary
plist = [];
numPoints = 0;

% Circle 1
center = [0.5 0.5]*NY;
radius = 0.1/pi*3.75*NY;
dr = 2.0; % segmentation length, default = 1 lattice unit
np = round(2*pi*radius/dr); % num of L-point
dr = 2*pi*radius/np;

% Gather->lagrangian point list
% start from phase = 0
numPoints = numPoints + np;
for i = 0:np-1
    phi = i/np*2*pi;
    p = center + (radius).*[sin(phi) cos(phi)];
    plist = [plist;p];
end

% Circle 2
center = [0.8 0.7]*NY;
radius = 0.05/pi*3.75*NY;
dr = 2.0; % segmentation length, default = 1 lattice unit
np = round(2*pi*radius/dr); % num of L-point
dr = 2*pi*radius/np;

numPoints = numPoints + np;
for i = 0:np-1
    phi = i/np*2*pi;
    p = center + (radius).*[sin(phi) cos(phi)];
    plist = [plist;p];
end

% Circle 3
center = [0.8 0.3]*NY;
radius = 0.05/pi*3.75*NY;
dr = 2.0; % segmentation length, default = 1 lattice unit
np = round(2*pi*radius/dr); % num of L-point
dr = 2*pi*radius/np;

numPoints = numPoints + np;
for i = 0:np-1
    phi = i/np*2*pi;
    p = center + (radius).*[sin(phi) cos(phi)];
    plist = [plist;p];
end

% Lattice coordinate offset
plist = plist + [0.5 0.5];
r_star = zeros(numPoints,1);
u_star = zeros(numPoints,1);
v_star = zeros(numPoints,1);

% initialise rho  ux  uy fields.
[rho, ux, uy, alpha, fx, fy] = plain(NX, NY, rho0, uwall,uc, 1, mask);
%  initialise f[2] as equilibrium for rho  ux  uy
f = init_equilibrium(rho, ux, uy, wi, dirx, diry, NX, NY, ndir, 1);

% Visualization
h2 = pcolor(X);h2.XData = X;h2.YData = Y;h2.ZData = 0*X;h2.EdgeColor='none';h2.FaceColor='interp';
colormap(desaturated_rainbow);
%l2 = plot(Y(1,:),Y(1,:),'-+');
fprintf('NSTEPS=%d, FRAMEGAP=%d\n',NSTEPS,PLOTGAP);
fprintf('Umax=%f\n',u_max);
fprintf('tau=%f\n',tau);
fprintf('Reg=%f\n',Reg);
fprintf('Re=%f\n',Re);

% buf to store temporary field
buf = zeros(NX*NY,ndir); % gpuArray

% Transfer to gpu
wi = gpuArray(wi); % gpuArray
dirx = gpuArray(dirx);
diry = gpuArray(diry);
bi   = gpuArray(bi);
buf = gpuArray(buf); % gpuArray
north = gpuArray(north);
south = gpuArray(south);
east = gpuArray(east);
west = gpuArray(west);
north2 = gpuArray(north2);
south2 = gpuArray(south2);
east2 = gpuArray(east2);
west2 = gpuArray(west2);
mask = gpuArray(mask);

% Assign boundaries to patches
%p_ind1 = [find(east);find(west)];
%p_ind2 = [find(west2);find(east2)];

% IB-LBM
plist = gpuArray(plist);
r_star = gpuArray(r_star);
u_star = gpuArray(u_star);
v_star = gpuArray(v_star);

%% Load kernels
cudaFilename = 'd2q9.cu';
ptxFilename = 'd2q9.ptx'; %,parallel.gpu.ptxext

% load streaming kernel
streaming = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processStreaming' );
% Make sure we have sufficient blocks to cover all of the locations
streaming.ThreadBlockSize = [streaming.MaxThreadsPerBlock,1,1];
streaming.GridSize = [ceil(numElements/streaming.MaxThreadsPerBlock),1];

% load collision kernel
collision = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processCollision' );
% Make sure we have sufficient blocks to cover all of the locations
collision.ThreadBlockSize = streaming.ThreadBlockSize;
collision.GridSize = streaming.GridSize;

% load compute_rho_and_u kernel
compute_rho_and_u = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processCompute' );
% Make sure we have sufficient blocks to cover all of the locations
compute_rho_and_u.ThreadBlockSize = streaming.ThreadBlockSize;
compute_rho_and_u.GridSize = streaming.GridSize;

% load compute_ibnode kernel (immersed boundary)
compute_ibnode = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processComputeIbNode' );
% Make sure we have sufficient blocks to cover all of the locations
%compute_ibnode.ThreadBlockSize = [16,1,1]; % 16 neighbours of L-point
%compute_ibnode.GridSize = [1,1];
compute_ibnode.ThreadBlockSize = [streaming.MaxThreadsPerBlock,1,1];
compute_ibnode.GridSize = [ceil(numPoints/streaming.MaxThreadsPerBlock),1];

% load interp_ibnode kernel (immersed boundary)
interp_ibnode = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processInterpIbNode' );
% Make sure we have sufficient blocks to cover all of the locations
interp_ibnode.ThreadBlockSize = [streaming.MaxThreadsPerBlock,1,1];
interp_ibnode.GridSize = [ceil(numPoints/streaming.MaxThreadsPerBlock),1];
%interp_ibnode.ThreadBlockSize = streaming.ThreadBlockSize;
%interp_ibnode.GridSize = streaming.GridSize;

% load compute_ibnode kernel (legacy, immersed boundary)
compute_ibnodeS = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processComputeIbNodeS' );
% Make sure we have sufficient blocks to cover all of the locations
compute_ibnodeS.ThreadBlockSize = [16,1,1]; % 16 neighbours of L-point
compute_ibnodeS.GridSize = [1,1];

% load bc kernel
correct_noslip = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processNoslip' );
% Make sure we have sufficient blocks to cover all of the locations
correct_noslip.ThreadBlockSize = streaming.ThreadBlockSize;
correct_noslip.GridSize = streaming.GridSize;

correct_periodic = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processPeriodic' );
% Make sure we have sufficient blocks to cover all of the locations
correct_periodic.ThreadBlockSize = streaming.ThreadBlockSize;
correct_periodic.GridSize = streaming.GridSize;

correct_outlet = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processNRoutlet' );
% Make sure we have sufficient blocks to cover all of the locations
correct_outlet.ThreadBlockSize = streaming.ThreadBlockSize;
correct_outlet.GridSize = streaming.GridSize;

correct_Routlet = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processRoutlet' );
% Make sure we have sufficient blocks to cover all of the locations
correct_Routlet.ThreadBlockSize = streaming.ThreadBlockSize;
correct_Routlet.GridSize = streaming.GridSize;

% Set constants
setConstantMemory(streaming, ...
    'tileX', int32(NX), 'tileY', int32(NY), 'nElem', int32(numElements),...
    'nDir', int32(ndir), 'dirX', int32(dirx), 'dirY', int32(diry),...
    'bi', int32(bi-1),'wi',wi);
setConstantMemory(streaming,'cs',[csSqrInv csSqr sqrt(csSqrInv) sqrt(csSqr)]);
setConstantMemory(streaming,'tau',[tau 1/tau 1.0-1/tau]);
%setConstantMemory(streaming,'constant_Fx',wi.*dirx*Fx*csSqrInv);
%setConstantMemory(streaming,'constant_Fy',wi.*diry*Fy*csSqrInv);

% entropy constraint
setConstantMemory(streaming,'enableEntropyConstraint',true);
setConstantMemory(streaming,'stableDeviation',0.3); % magic number
setConstantMemory(streaming,'alphaMin',0.2); % magic number

%%  main simulation loop; take NSTEPS time steps
convergence = Inf;
n = 0;
while((~isfinite(convergence) || convergence > 1e-10))
    
    % record statistics
    uhist = ux;
    vhist = uy;
    
    % swap buffer
    s.read = mod(n,2)+1;
    s.write = mod(s.read,2)+1;
   
    % compute rho and u
    [rho, ux, uy] = feval(compute_rho_and_u,rho,ux,uy,f(:,:,s.read),mask);
    
    % Immersed Boundary - direct forcing
    fx = 0*fx;fy = 0*fy; % Erase previous result
    
    % 6-DOF (in lattice unit)
    %usolid = uc; % translation
    usolid = 0.66*uc*sin(2*pi*n/2000); % sinusoidal
    vsolid = 0.34*uc*sin(2*pi*n/2000 + pi/2); % sinusoidal
    plist = plist + [usolid vsolid];
        
    for m = 1:1
        % Interp Eulerian Vel to Lagrangian
        %[r_star,u_star,v_star] = feval(interp_ibnode,r_star,u_star,v_star,rho,ux,uy,fx,fy,plist(:,1),plist(:,2),numPoints);
        
        Fu = griddedInterpolant(gather(reshape(ux + 0.5*fx./rho,[],NY)));
        Fv = griddedInterpolant(gather(reshape(uy + 0.5*fy./rho,[],NY)));
        u_star = 0*u_star + Fu(gather(plist+1));
        v_star = 0*v_star + Fv(gather(plist+1));
        % Compute&Distribute IB force
        % Caution: low-efficiency implementation
        % Require atomic optimisation
        %for i = 1:numPoints
        %    [fx,fy] = feval(compute_ibnodeS,fx,fy,u_star(i) - usolid,v_star(i) - vsolid,plist(i,1),plist(i,2),i-1,dr);
        %end
        
        % Atomic version
        [fx,fy] = feval(compute_ibnode,fx,fy,u_star-usolid,v_star-vsolid,plist(:,1),plist(:,2),dr,numPoints);
    end
    %fprintf("Max force [%f, %f],Min force [%f, %f]\n",max(fx),max(fy),min(fx),min(fy));
    % Immersed Boundary - direct forcing - end

    % Collision
    %[f(:,:,s.read), alpha] = feval(collision, f(:,:,s.read), alpha,f(:,:,s.read), rho, ux, uy, mask);
    [f(:,:,s.read), alpha] = feval(collision, f(:,:,s.read),alpha, f(:,:,s.read), rho, ux, uy,fx,fy, mask);
    
    % periodic pressure drop
    %a1 = (1 - drho./rho(west2));rhop1 = rho(west2)-mean(rho(west2))+rho0+drho/2;
    %a2 = (1 + drho./rho(east2));rhop2 = rho(east2)-mean(rho(east2))+rho0-drho/2;
    
    %f(p_ind1,:,s.read) = feval(correct_periodic,...
    %    f(p_ind1,:,s.read), f(p_ind2,:,s.read), [a1;a2], [rhop1;rhop2],...
    %    rho(p_ind2), ux(p_ind2), uy(p_ind2), numel(p_ind1));
    
    % periodic without pressure drop
    %f(east,:,s.read) = f(west2,:,s.read);
    %f(west,:,s.read) = f(east2,:,s.read);
    
    % noslip: leaking
    f(:,:,s.write) = feval(correct_noslip, f(:,:,s.write), f(:,:,s.read), rho, ux, uy, north|south|west);
    
    % CBC outlet
    f(:,:,s.read) = feval(correct_outlet, f(:,:,s.read), rho, ux, uy, east);
    
    % R outlet (for comparison)
    %f(:,:,s.read) = feval(correct_Routlet,f(:,:,s.read), f(:,:,s.read), rho, ux, uy, east);
    
    % Streaming
    f(:,:,s.write) = feval(streaming, f(:,:,s.read), f(:,:,s.read), rho, ux, uy, mask);
    
    % Advance
    n = n+1;
    %min_alpha = min(min_alpha, min(alpha(:)));
    
    if(mod(n,PLOTGAP) == 0 || n == NSTEPS)
        % Convergence
        conv = sum((ux - uhist).*(ux - uhist) + (uy - vhist).*(uy - vhist));
        convergence = sqrt(conv/sum(uhist.^2 + vhist.^2));
        fprintf('%e>\n', convergence);
        
        % Visualization
        rh = reshape(gather(rho),NX,NY);
        uh = reshape(gather(Cu*ux),NX,NY);
        vh = reshape(gather(Cu*uy),NX,NY);
        fxh = reshape(gather(fx),NX,NY);
        fyh = reshape(gather(fy),NX,NY);
        [dvdx,~] = gradient(vh);[dudy,~] = gradient(uh);
        h2.CData = sqrt(uh.^2 + vh.^2);
        %h2.CData = dvdx - dudy;
        %l2.YData = uh(floor(NX/2),:);
        %h2.CData = vh;
        %h2.CData = rh;
        axis equal;axis tight;axis([0 ratio(1) 0 ratio(2)]);
        %quiver(plist(:,1),plist(:,2),u_star,v_star)
        title(sprintf('%d',n/PLOTGAP))%,caxis([-1 1])
        drawnow;
        
    end
end

function [r,u,v, alpha, fx, fy] = plain(NX,NY,rho0,uw,uc, IsgpuArray, mask)

if(nargin < 6)
    IsgpuArray = 0;
end

r = zeros(NX,NY);
u = zeros(NX,NY);
v = zeros(NX,NY);
alpha = zeros(NX,NY)+2.0;
fx = zeros(NX,NY);
fy = zeros(NX,NY);

r = r + rho0;

% Input bc
u([1 end],:) = uc;

r = r(:);
u = u(:);
v = v(:);
alpha = alpha(:);
fx = fx(:);
fy = fy(:);

% Set multiphase internal field
%G = -6 rhog = 0.056 rhol = 2.659
%r(~mask) = rand(numel(r(~mask)),1)*(2.659-0.056)+0.056;
%r(~mask) = 1.0+1e-2*rand(numel(r(~mask)),1);

% wetting angle
%r(mask) = 1;
%u(~mask) = uc;

if(IsgpuArray)
    r = gpuArray(r);
    u = gpuArray(u);
    v = gpuArray(v);
    alpha = gpuArray(alpha);
    fx = gpuArray(fx);
    fy = gpuArray(fy);
end

end

function buf = init_equilibrium(rho, ux, uy, wi, dirx, diry, NX, NY, NDIR, IsgpuArray)

if(nargin < 10)
    IsgpuArray = 0;
end

buf =zeros(NX*NY, NDIR, 2); % gpuArray

if(IsgpuArray)
    buf = gpuArray(buf);
end

f=ux*dirx+uy*diry;
f=(3+4.5*f).*f;
f=bsxfun(@minus,f,1.5*(ux.^2+uy.^2));
f=bsxfun(@times,1+f,wi);
f=bsxfun(@times,f,rho);

buf(:,:,1) = f;

end


function feq = compute_equilibrium(rho, ux, uy, dirx, diry, wi)
feq=ux*dirx+uy*diry;
feq=(3+4.5*feq).*feq;
feq=bsxfun(@minus,feq,1.5*(ux.^2+uy.^2));
feq=bsxfun(@times,1.0+feq,wi);
feq=bsxfun(@times,feq,rho);
end
