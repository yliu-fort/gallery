clear; clc; close all;

%%
scale = 4; ratio = [3 1];
NX = ratio(1)*32*scale+2;NY = ratio(2)*32*scale+2;
%NX = 5;NY = 11;
numElements = NX*NY;
[X, Y] = ndgrid([0 linspace(0+1/2/(NX-2),ratio(1)-1/2/(NX-2),NX-2) ratio(1)],...
    [0 linspace(0+1/2/(NY-2),ratio(2)-1/2/(NY-2),NY-2) ratio(2)]);

% The kinematic viscosity
csSqrInv = 3.0;
csSqr = 1.0/csSqrInv;
c = sqrt(csSqr);
%tau = 0.9; % L2err = 0.0017 (H=13) 7.1686e-05 (H = 55)
%tau = 1.0; % L2err = 0.0038 (H=13)
%tau = sqrt(3/16)+1/2; % L2err = 8.2703e-09
tau = 0.501;

% Nondimensionalization
length_ui = 1.0; % m
nu_ui = 1.0/40000; % m^2/sec
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

cylinder=((X-0.5).^2 + (Y-0.5).^2) < 0.1^2;
cylinder = cylinder(:);
mask = mask | cylinder ;

% initialise rho  ux  uy fields.
[rho, ux, uy] = plain(NX, NY, rho0, uwall,uc, 1, mask);
%  initialise f[2] as equilibrium for rho  ux  uy
f = init_equilibrium(rho, ux, uy, wi, dirx, diry, NX, NY, ndir, 1);

% Visualization
h2 = pcolor(X);h2.XData = X;h2.YData = Y;h2.ZData = 0*X;h2.EdgeColor='none';h2.FaceColor='interp';
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
p_ind1 = [find(east);find(west)];
p_ind2 = [find(west2);find(east2)];

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
setConstantMemory(streaming,'constant_Fx',wi.*dirx*Fx*csSqrInv);
setConstantMemory(streaming,'constant_Fy',wi.*diry*Fy*csSqrInv);

% entropy constraint
setConstantMemory(streaming,'enableEntropyConstraint',true);
setConstantMemory(streaming,'stableDeviation',0.3); % magic number
setConstantMemory(streaming,'alphaMin',0.5); % magic number

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
    
    % Collision
    %[f(:,:,s.read), alpha] = feval(collision, f(:,:,s.read), alpha,f(:,:,s.read), rho, ux, uy, mask);
    f(:,:,s.read) = feval(collision, f(:,:,s.read),f(:,:,s.read), rho, ux, uy, mask);
    
    % noslip
    f(:,:,s.write) = feval(correct_noslip, f(:,:,s.write), f(:,:,s.read), rho, ux, uy, north|south|west);

    % CBC outlet
    f(:,:,s.read) = feval(correct_outlet, f(:,:,s.read), rho, ux, uy, east);

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
        %h2.CData = sqrt(uh.^2 + vh.^2);
        [dvdx,~] = gradient(vh);[dudy,~] = gradient(uh);
        h2.CData = dvdx - dudy;
        %l2.YData = uh(floor(NX/2),:);
        %h2.CData = vh;
        %h2.CData = rh-rho0;
        axis equal;axis tight;axis([0 ratio(1) 0 ratio(2)]);
        title(sprintf('%d',n/PLOTGAP))%,caxis([-1 1])
        drawnow;
        
    end
end
rho = reshape(rho, NX, NY);
ux = reshape(ux, NX, NY);
uy = reshape(uy, NX, NY);

function [r,u,v] = plain(NX,NY,rho0,uw,uc, IsgpuArray, mask)

if(nargin < 6)
    IsgpuArray = 0;
end

r = zeros(NX,NY);
u = zeros(NX,NY);
v = zeros(NX,NY);

r = r + rho0;

% Moving wall on north
u(:,NY) = uw;

% Inlet west
u([1 end],:) = uc;
%r(1,:) = 0.5774/(0.5774+uc);

% perturbation
%r(floor(NX/2),floor(NY/2)) = rho0 + 0.01;

r = r(:);
u = u(:);
v = v(:);

% initial internal field
%u(~mask) = uc;% + 1e-13*rand(sum(~mask),1);
%v(floor(NX*NY/4)) = 1e-13;v(floor(NX*NY/4)+1) = -1e-13;
%r(~mask) = rho0+1e-7*rand(sum(~mask),1);


if(IsgpuArray)
    r = gpuArray(r);
    u = gpuArray(u);
    v = gpuArray(v);
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
