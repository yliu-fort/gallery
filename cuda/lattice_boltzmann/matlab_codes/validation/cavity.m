clear; clc; close all;

%%
scale = 4;
NX = 32*scale+2;NY = 32*scale+2;
numElements = NX*NY;

% The kinematic viscosity
csSqrInv = 3.0;
csSqr = 1.0/csSqrInv;
c = sqrt(csSqr);
tau = 0.8; % taustar

% Nondimensionalization
length_ui = 1.0; % m
nu_ui = 1.0/100; % m^2/sec
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
uwall = u_max;

%  compute Taylor-Green flow at t=0
%  to initialise rho  ux  uy fields.
[rho, ux, uy] = plain(NX, NY, rho0, uwall,0, 1);
%  initialise f[2] as equilibrium for rho  ux  uy
f = init_equilibrium(rho, ux, uy, wi, dirx, diry, NX, NY, ndir, 1);

% Boundary
north=false(NX,NY);
north(:,end) = true;
north = north(:);% gpuArray

south=false(NX,NY);
south(:,1) = true;
south = south(:);% gpuArray

east=false(NX,NY);
east(end,:) = true;
east = east(:);

west=false(NX,NY);
west(1,:) = true;
west = west(:);

north2=false(NX,NY);
north2(:,end-1) = true;
north2 = north2(:);% gpuArray

south2=false(NX,NY);
south2(:,2) = true;
south2 = south2(:);% gpuArray

east2=false(NX,NY);
east2(end-1,:) = true;
east2 = east2(:);

west2=false(NX,NY);
west2(2,:) = true;
west2 = west2(:);

mask = east | west | north | south;

% Visualization
[X, Y] = ndgrid([0 linspace(0+1/2/(NX-2),1-1/2/(NX-2),NX-2) 1],...
    [0 linspace(0+1/2/(NY-2),1-1/2/(NY-2),NY-2) 1]);
subplot(2,2,1)
h1 = quiver([],[],[],[],1.5,'b');h1.XData = X(:);h1.YData = Y(:);axis equal;axis tight
subplot(2,2,2)
%h2 = pcolor(X);h2.XData = X;h2.YData = Y;h2.ZData = 0*X;h2.EdgeColor='none';h2.FaceColor='interp';
%h2 = scatter([],[],[],[],'.');h2.XData = X(:);h2.YData = Y(:);axis equal;axis tight
%s1 = semilogy([1]);hold on, s2 = semilogy([1]);
s3 = semilogy([1]);hold on,s4 = semilogy([1]);

fprintf('NSTEPS=%d\n',NSTEPS);
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

% Validation
T1 = load('u_verticalline.mat');
T2 = load('v_horizontalline.mat');
ind.u = T1.index(:,2);
ind.v = T2.index(:,2);
uRef = T1.Re_100;
vRef = T2.Re_100;

line1=false(NX,NY);
line1(floor(NX/2),:) = true;
line1 = gpuArray(line1(:));

line2=false(NX,NY);
line2(:,floor(NY/2)) = true;
line2 = gpuArray(line2(:));

subplot(2,2,3)
l1 = plot(X(:,1),Cu*ux(line1));hold on,plot(ind.u,uRef,'+')
subplot(2,2,4)
l2 = plot(Y(1,:),Cu*uy(line2));hold on,plot(ind.v,vRef,'+')

%% Load kernels
cudaFilename = 'd2q9.cu';
ptxFilename = 'd2q9.ptx'; %,parallel.gpu.ptxext
streaming = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename, 'processStreaming' );
setConstantMemory(streaming, ...
    'tileX', int32(NX), 'tileY', int32(NY), 'nElem', int32(numElements),...
	'nDir', int32(ndir), 'dirX', int32(dirx), 'dirY', int32(diry),...
    'bi', int32(bi-1),'wi',wi);
setConstantMemory(streaming,'cs',[csSqrInv csSqr sqrt(csSqrInv) sqrt(csSqr)]);
setConstantMemory(streaming,'tau',[tau 1/tau 1.0-1/tau]);
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

%%  main simulation loop; take NSTEPS time steps
convergence = Inf;
n = 0;
while((~isfinite(convergence) || convergence > 1e-10) && n < NSTEPS)
    
    % record statistics
    uhist = ux;
    vhist = uy;
    
    % swap buffer
    selector.read = mod(n,2)+1;
    selector.write = mod(selector.read,2)+1;
    
    % compute rho and u
    [rho, ux, uy] = feval(compute_rho_and_u,rho,ux,uy,f(:,:,selector.read),mask);
    
    % Collision
    f(:,:,selector.read) = feval(collision, f(:,:,selector.read), f(:,:,selector.read), rho, ux, uy, mask);
    
    % Streaming
    f(:,:,selector.write) = feval(streaming, f(:,:,selector.read), f(:,:,selector.read), rho, ux, uy, mask);

    % Advance in time
    n = n+1;
        
    if(mod(n,scale*8) == 0 || n == NSTEPS)
        % Convergence
        conv = sum((ux - uhist).*(ux - uhist) + (uy - vhist).*(uy - vhist));
        convergence = sqrt(conv/sum(uhist.^2 + vhist.^2));  
        fprintf('%e>\n', convergence);
        
        % Visualization
        uh = gather(Cu*ux);vh = gather(Cu*uy);
        h1.UData = uh;h1.VData = vh;
        h2.CData = sqrt(uh.*uh + vh.*vh);
        l1.YData = uh(line1);
        l2.YData = vh(line2);
        
        prop = compute_flow_properties(1.0, ...
        interp1(X(:,1) ,uh(line1),ind.u),...
        interp1(Y(1,:)',vh(line2),ind.v),...
        uRef, vRef)';
    
        if(n == 0)
        %s1.YData = prop(1,1);s2.YData = prop(2,1);
        s3.YData = prop(3,1);s4.YData = gather(convergence);
        else
        %s1.YData = [s1.YData prop(1,n+1)];s2.YData = [s2.YData prop(2,n+1)];
        s3.YData = [s3.YData prop(3,1)];
        s4.YData = [s4.YData gather(convergence)];
        end
        %axis equal;axis tight;axis([0 1 0 1]);
        %title(sprintf('%d',n))
        drawnow;
    end
end
rho = reshape(rho, NX, NY);
ux = reshape(ux, NX, NY);
uy = reshape(uy, NX, NY);

function [r,u,v] = plain(NX,NY,rho0,uw,uc, IsgpuArray)

if(nargin < 6)
    IsgpuArray = 0;
end

r = zeros(NX,NY);
u = zeros(NX,NY);
v = zeros(NX,NY);

r = r + rho0;

% Moving wall on north
u(:,NY) = uw;

r = r(:);
u = u(:)+uc;
v = v(:);

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