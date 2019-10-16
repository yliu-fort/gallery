clear; clc; close all;

%%
scale = 4; ratio = [3 1 1];
NX = ratio(1)*16*scale+2;NY = ratio(2)*16*scale+2;NZ = ratio(3)*16*scale+2;
numElements = NX*NY*NZ;

% The kinematic viscosity
csSqrInv = 3.0;
csSqr = 1.0/csSqrInv;
c = sqrt(csSqr);
tau = 0.51; % taustar

% Nondimensionalization
length_ui = 1.0; % m
nu_ui = 1.0/2000; % m^2/sec
time_length_ui = 150.0; % sec
u_max_ui = 1.0;

Cl = length_ui/(NY - 2)/1.0; % boundary node take 2 * dx/2
nu = csSqr*(tau-0.5); % nu = cs^2*(tao - 0.5)
Ct =  nu/nu_ui*Cl^2; % dtstar = 1.0
NSTEPS = floor(time_length_ui / Ct);
%NSTEPS = 64*scale*scale; % specific timestep for validation
Cu = Cl/Ct;
u_max = u_max_ui/Cu; % should below 0.1 or 0.03, max < 0.4
%Ca = Cl^2/Ct;

% D3Q27
ndir = 27;
% The lattice weights
w0 = 8.0/27.0; % zero weight
ws = 2.0/27.0; % adjacent weight
wd = 1.0/54.0; % diagonal weight
wc = 1.0/216.0; % diagonal weight
% Arrays of the lattice weights and direction components
wi = ([w0 ws ws ws ws ws ws wd wd wd wd wd wd wd wd wd wd wd wd wc wc wc wc wc wc wc wc]); % gpuArray
dirx = ([0 1 -1  0  0  0  0  1 -1  1 -1 1 -1  1 -1 0  0  0  0 1 -1  1 -1 -1  1 -1  1]);
diry = ([0 0  0  1 -1  0  0  1 -1 -1  1 0  0  0  0 1 -1  1 -1 1 -1 -1  1  1 -1 -1  1]);
dirz = ([0 0  0  0  0  1 -1  0  0  0  0 1 -1 -1  1 1 -1 -1  1 1 -1  1 -1  1 -1  1 -1]);
bi   = ([1 3 2 5 4 7 6 9 8 11 10 13 12 15 14 17 16 19 18 21 20 23 22 25 24 27 26]);

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

%  compute Taylor-Green flow at t=0
%  to initialise rho  ux  uy fields.
[rho, ux, uy, uz] = plain(NX, NY,NZ, rho0, uwall,uc, 1);
%  initialise f[2] as equilibrium for rho  ux  uy
f = init_equilibrium(rho, ux, uy, uz, wi, dirx, diry,dirz, NX, NY,NZ, ndir, 1);

[X, Y, Z] = ndgrid(...
    [0 linspace(0+1/2/(NX-2),1-1/2/(NX-2),NX-2) 1].*ratio(1),...
    [0 linspace(0+1/2/(NY-2),1-1/2/(NY-2),NY-2) 1].*ratio(2),...
    [0 linspace(0+1/2/(NZ-2),1-1/2/(NZ-2),NZ-2) 1].*ratio(3));

% Boundary
north=false(NX,NY,NZ);
north(:,end,:) = true;
north = north(:);% gpuArray

south=false(NX,NY,NZ);
south(:,1,:) = true;
south = south(:);% gpuArray

east=false(NX,NY,NZ);
east(end,:,:) = true;
east = east(:);

west=false(NX,NY,NZ);
west(1,:,:) = true;
west = west(:);

front=false(NX,NY,NZ);
front(:,:,1) = true;
front = front(:);

back=false(NX,NY,NZ);
back(:,:,end) = true;
back = back(:);

cylinder=(((X-0.35).^2 + (Y-0.5).^2 + (Z-0.5).^2) < (0.2^2));
cylinder = cylinder(:);


% Visualization
%subplot(2,2,1)
%h1 = quiver3([],[],[],[],[],[],2,'b');h1.XData = X(:);h1.YData = Y(:);h1.ZData = Z(:);axis equal;axis tight
%subplot(2,2,2)
%h2 = pcolor(X);h2.XData = X;h2.YData = Y;h2.ZData = 0*X;h2.EdgeColor='none';h2.FaceColor='interp';
%h2 = scatter([],[],[],[],'.');h2.XData = X(:);h2.YData = Y(:);axis equal;axis tight
%s1 = semilogy([1]);hold on, s2 = semilogy([1]);
%s3 = semilogy([1]);hold on,s4 = semilogy([1]);

fprintf('NSTEPS=%d\n',NSTEPS);
fprintf('Umax=%f\n',u_max);
fprintf('tau=%f\n',tau);
fprintf('Reg=%f\n',Reg);
fprintf('Re=%f\n',Re);

% buf to store temporary field
buf = zeros(NX*NY*NZ,ndir); % gpuArray

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
front = gpuArray(front);
back = gpuArray(back);
cylinder = gpuArray(cylinder);

mask = east | west | north | south | front | back | cylinder;
noslip = west | north | south | front | back | cylinder;
outlet = east;

%% Load kernels
cudaFilename = 'd3q27.cu';
ptxFilename = 'd3q27.ptx'; %,parallel.gpu.ptxext
streaming = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processStreaming' );
setConstantMemory(streaming, ...
    'tileX', int32(NX), 'tileY', int32(NY), 'tileZ', int32(NZ), 'nElem', int32(numElements),...
	'nDir', int32(ndir), 'dirX', int32(dirx), 'dirY', int32(diry), 'dirZ', int32(dirz),...
    'bi', int32(bi-1),'wi',wi);

% Make sure we have sufficient blocks to cover all of the locations
streaming.ThreadBlockSize = [streaming.MaxThreadsPerBlock,1,1];
streaming.GridSize = [ceil(numElements/streaming.MaxThreadsPerBlock),1];

%
collision = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processCollision' );
% Make sure we have sufficient blocks to cover all of the locations
collision.ThreadBlockSize = streaming.ThreadBlockSize;
collision.GridSize = streaming.GridSize;

%
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

correct_Routlet = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processRoutlet' );
% Make sure we have sufficient blocks to cover all of the locations
correct_Routlet.ThreadBlockSize = streaming.ThreadBlockSize;
correct_Routlet.GridSize = streaming.GridSize;

setConstantMemory(streaming,'cs',[csSqrInv csSqr sqrt(csSqrInv) sqrt(csSqr)]);
setConstantMemory(streaming,'tau',[tau 1/tau 1.0-1/tau]);
setConstantMemory(streaming,'constant_Fx',wi.*dirx*0*csSqrInv);
setConstantMemory(streaming,'constant_Fy',wi.*diry*0*csSqrInv);
setConstantMemory(streaming,'constant_Fz',wi.*dirx*0*csSqrInv);

% entropy constraint
setConstantMemory(streaming,'enableEntropyConstraint',false);
setConstantMemory(streaming,'stableDeviation',10); % magic number
setConstantMemory(streaming,'alphaMin',0.2); % magic number

%%  main simulation loop; take NSTEPS time steps
convergence = Inf;
n = 0;
while((~isfinite(convergence) || convergence > 1e-10) && n < NSTEPS)
    
    % record statistics
    uhist = ux;
    vhist = uy;
    whist = uz;
    
    % swap buffer
    selector.read = mod(n,2)+1;
    selector.write = mod(selector.read,2)+1;
    
    % compute rho and u
    [rho, ux, uy, uz] = feval(compute_rho_and_u,rho,ux,uy,uz,f(:,:,selector.read),mask);
    
    % Collision
    f(:,:,selector.read) = feval(collision, f(:,:,selector.read), f(:,:,selector.read), rho, ux, uy,uz, mask);
    
    % noslip
    f(:,:,selector.write) = feval(correct_noslip, f(:,:,selector.write), f(:,:,selector.read), rho, ux, uy, uz, noslip);

    % CBC outlet
    %f(:,:,s.read) = feval(correct_outlet, f(:,:,s.read), rho, ux, uy, east);

    % R outlet (for comparison)
    f(:,:,selector.read) = feval(correct_Routlet,f(:,:,selector.read), f(:,:,selector.read), rho, ux, uy, uz, outlet);
    
    % Streaming
    f(:,:,selector.write) = feval(streaming, f(:,:,selector.read), f(:,:,selector.read), rho, ux, uy,uz, mask);

    % Advance in time
    n = n+1;
        
    if(mod(n,scale*8) == 0 || n == NSTEPS)
        % Convergence
        conv = sum((ux - uhist).*(ux - uhist) + (uy - vhist).*(uy - vhist) + (uz - whist).*(uz - whist));
        convergence = sqrt(conv/sum(uhist.^2 + vhist.^2 + whist.^2));  
        fprintf('%e>\n', convergence);
        
        % Visualization
        uh = gather(Cu*ux);vh = gather(Cu*uy);wh = gather(Cu*uz);
        
        clf;fid = figure(1);
        set(fid,'Units','pixels','Position',[1 1 1920 1080]);
        
        %h1.UData = uh-Cu*uc;h1.VData = vh;h1.WData = wh;
        umag = reshape(sqrt(uh.^2 + vh.^2 + wh.^2), NX, NY, NZ);
        umag = permute(umag,[2 1 3]);
        h=slice(umag,floor(3*NX/4),floor(NY/2),floor(NZ/3));
        h(1).EdgeColor='none';
        h(2).EdgeColor='none';
        h(3).EdgeColor='none';
        % lambda2
        hold on
        [dudx,dudy,dudz]=gradient(reshape(uh,NX,NY,NZ));
        [dvdx,dvdy,dvdz]=gradient(reshape(vh,NX,NY,NZ));
        [dwdx,dwdy,dwdz]=gradient(reshape(wh,NX,NY,NZ));
        l2 = arrayfun(@lambda2,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz);
        l2 = permute(l2,[2 1 3]);
        
        p = patch(isosurface(l2, 0.002));isonormals(l2, p)
        p.FaceColor = 'none';p.EdgeColor = 'blue';
        camlight; lighting phong
        
        daspect([1 1 1]),view(3),axis equal
        colormap desaturated_rainbow
        drawnow;
        
        saveas(fid,['pic/result.' sprintf('%06d',n) '.png'])

    end
end
rho = reshape(rho, NX, NY, NZ);
ux = reshape(ux, NX, NY, NZ);
uy = reshape(uy, NX, NY, NZ);
uz = reshape(uz, NX, NY, NZ);

function [r,u,v,w] = plain(NX,NY,NZ,rho0,uw,uc, IsgpuArray)

if(nargin < 7)
    IsgpuArray = 0;
end

r = zeros(NX,NY,NZ);
u = zeros(NX,NY,NZ);
v = zeros(NX,NY,NZ);
w = zeros(NX,NY,NZ);

r = r + rho0;

% Moving wall on north
u(:,:,NZ) = uw;

% Inlet on west
u(1,:,:) = uc;
u(:,1,:) = uc;
u(:,NY,:) = uc;
u(:,:,1) = uc;
u(:,:,NZ) = uc;

r = r(:);
u = u(:);
v = v(:);
w = w(:);

if(IsgpuArray)
    r = gpuArray(r);
    u = gpuArray(u);
    v = gpuArray(v);
    w = gpuArray(w);
end

end

function buf = init_equilibrium(rho, ux, uy, uz, wi, dirx, diry,dirz, NX, NY, NZ, NDIR, IsgpuArray)

if(nargin < 11)
    IsgpuArray = 0;
end

buf =zeros(NX*NY*NZ, NDIR, 2); % gpuArray

if(IsgpuArray)
    buf = gpuArray(buf);
end

f=ux*dirx+uy*diry+uz*dirz;
f=(3+4.5*f).*f;
f=bsxfun(@minus,f,1.5*(ux.^2+uy.^2+uz.^2));
f=bsxfun(@times,1+f,wi);
f=bsxfun(@times,f,rho);

buf(:,:,1) = f;

end

function prop = compute_flow_properties(r, u, v, uxat, uyat)

% prop must point to space for 4 doubles:
% 0: energy
% 1: L2 error in rho
% 2: L2 error in ux
% 3: L2 error in uy

E = sum(r.*(u.*u + v.*v));

%drho = (r - rhoat)./rhoat*100;
%drho(~isfinite(drho)) = 0;

conv = sum((u - uxat).^2 + (v - uyat).^2); 
prop(1) = norm(u - uxat, 2)/sqrt(sum(uxat.^2));
prop(2) = norm(v - uyat, 2)/sqrt(sum(uyat.^2));
prop(3) = sqrt(conv/sum(uxat.^2 + uyat.^2));
end
