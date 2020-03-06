clear; clc; close all;

%% Magic numbers
L = 10;
eta = L/91; % 991: unstable when t > 5
nk = 64;
E0 = 3.2;
lambda = 0.5;

LI = L/3;
TL = 0.2*LI/sqrt(E0);

rng(6356435)
theta = 2*pi*rand(nk,1);

k1 = 2*pi/L;
k_Nk = 2*pi/eta;
alpha = (L/eta)^(1/(nk-1));
k = k1*alpha.^((1:nk).'-1);
E = E0*L*(k*L).^(-5/3);
dk = gradient(k);

a = sqrt(E.*dk);
b = a;

omega = lambda*sqrt(k.^3.*E).';
dt = pi/8/omega(end);
dt = 1e-2;

%% Vectors
An = (a.*[ cos(theta) -sin(theta)]).';
Bn = (b.*[-cos(theta)  sin(theta)]).';
Kn = (k.*[ sin(theta)  cos(theta)]).';


u = @(x,y,t)(sum(An(1,:).*cos( x*Kn(1,:) + y*Kn(2,:) + omega.*t ) ...
      + Bn(1,:).*sin( x*Kn(1,:) + y*Kn(2,:) + omega.*t ),2));
v = @(x,y,t)(sum(An(2,:).*cos( x*Kn(1,:) + y*Kn(2,:) + omega.*t ) ...
      + Bn(2,:).*sin( x*Kn(1,:) + y*Kn(2,:) + omega.*t ),2));
  
%% Configuration

nx = 64;
ny = 64;
%[X, Y]= meshgrid((0:nx-1)/k_Nk+1, (0:ny-1)/k_Nk+1);
[X, Y]= meshgrid((0:nx-1)/k(end)+1, (0:ny-1)/k(end)+1);
x = X(:);
y = Y(:);
n = nx*ny;

pos(1).x = x;
pos(1).y = y;
time(1) = 0;

%% Setup simulation
t = 0;
i = 1;
tol = 1e-2;
s = 1;
%while (t < 100)
    
    % 2nd-order runge-kutta
    %a1 = [u(x,y,t) v(x,y,t)];
    %a2 = [u(x+0.5*dt*a1(:,1),y+0.5*dt*a1(:,2),t+0.5*dt) v(x+0.5*dt*a1(:,1),y+0.5*dt*a1(:,2),t+0.5*dt)];
    %dx = dt*a2(:,1);
    %dy = dt*a2(:,2);
    
    % 4th-order runge-kutta
    %a1 = dt*[u(x,y,t) v(x,y,t)];
    %a2 = dt*[u(x+0.5*a1(:,1),y+0.5*a1(:,2),t+0.5*dt) v(x+0.5*a1(:,1),y+0.5*a1(:,2),t+0.5*dt)];
    %a3 = dt*[u(x+0.5*a2(:,1),y+0.5*a2(:,2),t+0.5*dt) v(x+0.5*a2(:,1),y+0.5*a2(:,2),t+0.5*dt)];
    %a4 = dt*[u(x+a3(:,1),y+a3(:,2),t+dt) v(x+a3(:,1),y+a3(:,2),t+dt)];
    %dxdt = 1.0/6.0*(a1(:,1) + 2*a2(:,1) + 2*a3(:,1) + a4(:,1));
    %dydt = 1.0/6.0*(a2(:,2) + 2*a2(:,2) + 2*a3(:,2) + a4(:,2));
    
    % 4/5th Runge-kutta adaptive
    while true
        
    a1 = dt*[u(x,y,t) v(x,y,t)];
    a2 = dt*[u(x+0.25*a1(:,1),y+0.25*a1(:,2),t+0.25*dt)...
             v(x+0.25*a1(:,1),y+0.25*a1(:,2),t+0.25*dt)];
    a3 = dt*[u(x+3.0/32.0*a1(:,1)+9.0/32.0*a2(:,1),y+3.0/32.0*a1(:,2)+9.0/32.0*a2(:,2),t+3.0/8.0*dt)...
             v(x+3.0/32.0*a1(:,1)+9.0/32.0*a2(:,1),y+3.0/32.0*a1(:,2)+9.0/32.0*a2(:,2),t+3.0/8.0*dt)];
    a4 = dt*[u(x+1932.0/2197.0*a1(:,1)+7200.0/2197.0*a2(:,1)+7296.0/2197.0*a3(:,1),y+1932.0/2197.0*a1(:,2)+7200.0/2197.0*a2(:,2)+7296.0/2197.0*a3(:,2),t+12.0/13.0*dt)...
             v(x+1932.0/2197.0*a1(:,1)+7200.0/2197.0*a2(:,1)+7296.0/2197.0*a3(:,1),y+1932.0/2197.0*a1(:,2)+7200.0/2197.0*a2(:,2)+7296.0/2197.0*a3(:,2),t+12.0/13.0*dt)];
    a5 = dt*[u(x+439.0/216.0*a1(:,1)-8.0*a2(:,1)+3680.0/513.0*a3(:,1)-845.0/4104.0*a4(:,1),...
               y+439.0/216.0*a1(:,2)-8.0*a2(:,2)+3680.0/513.0*a3(:,2)-845.0/4104.0*a4(:,2),t+dt) ...
             v(x+439.0/216.0*a1(:,1)-8.0*a2(:,1)+3680.0/513.0*a3(:,1)-845.0/4104.0*a4(:,1),...
               y+439.0/216.0*a1(:,2)-8.0*a2(:,2)+3680.0/513.0*a3(:,2)-845.0/4104.0*a4(:,2),t+dt)];
    a6 = dt*[u(x-8.0/27.0*a1(:,1)+2.0*a2(:,1)-3544.0/2565.0*a3(:,1)+1859.0/4104.0*a4(:,1)-11.0/40.0*a5(:,1),...
               x-8.0/27.0*a1(:,2)+2.0*a2(:,2)-3544.0/2565.0*a3(:,2)+1859.0/4104.0*a4(:,2)-11.0/40.0*a5(:,2),t+0.5*dt) ...
             v(x-8.0/27.0*a1(:,1)+2.0*a2(:,1)-3544.0/2565.0*a3(:,1)+1859.0/4104.0*a4(:,1)-11.0/40.0*a5(:,1),...
               x-8.0/27.0*a1(:,2)+2.0*a2(:,2)-3544.0/2565.0*a3(:,2)+1859.0/4104.0*a4(:,2)-11.0/40.0*a5(:,2),t+0.5*dt)];  
    
    dxdt = (25.0/216.0)*a1(:,1) + (1408.0/2565.0)*a3(:,1) + (2197.0/4101.0)*a4(:,1) - 0.2*a5(:,1);
    dydt = (25.0/216.0)*a1(:,2) + (1408.0/2565.0)*a3(:,2) + (2197.0/4101.0)*a4(:,2) - 0.2*a5(:,2);
           
    % Residual estimation
    ds = 0.002777777777778*a1-0.029941520467836*a3-0.029591504049595*a4-0.02*a5+0.036363636363636*a6;
    ds = max(abs(ds(:)));
    s = 0.84*(tol/ds)^0.25;

        if ( s < 1 )
            % h is too large
            dt = 0.5*dt;
        elseif ( s < 2 )
            % h is appropriate
            break;
        else % s >= 2
            % h is too small
            dt = 2*dt;
        end

    end
    
    t = t + dt;
    x = x + dxdt;
    y = y + dydt;    
    
    fprintf("Simulation time = %f, dt = %e, tolerance = %e\n", t, dt, ds);
    
    % Output
    %if(t > i*dt*10)
    i = i + 1;
    pos(i).x = x;
    pos(i).y = y;
    time(i) = t;
    
    % mean distance square
    d2 = [];
    for j = 1:nx-1
        d2 = [d2;(pos(i).x([1:ny]+(j-1)*ny)-pos(i).x([1:ny]+j*ny)).^2+(pos(i).y([1:ny]+(j-1)*ny)-pos(i).y([1:ny]+j*ny)).^2];
    end
    for j = 1:ny-1
        d2 = [d2;(pos(i).x([1:nx:end]+(j-1))-pos(i).x([1:nx:end]+j)).^2+(pos(i).y([1:nx:end]+(j-1))-pos(i).y([1:nx:end]+j)).^2];
    end
    d2b(i) = mean(d2);

    %end
    if(numel(d2b) > 16384),d2b = d2b(1:2:end);end

%end

%% Visualization
for i = 1:10:numel(pos)
    
    % Displacement
    x = x + dxdt;
    y = y + dydt;
    scatter(pos(i).x,pos(i).y,10,1:n,'.'),axis equal,axis([-10 10 -10 10])
    title(time(i))
    drawnow
   % pause
end    

%% Plot
loglog(time, d2b,'.')
hold on,loglog(time, time.^(3)/250)

%% Distance
%{
d2b = [];
for j = 1:numel(time)
    d2 = [];
    for i = 1:nx-1
        d2 = [d2;(pos(j).x([1:ny]+(i-1)*ny)-pos(j).x([1:ny]+i*ny)).^2+(pos(j).y([1:ny]+(i-1)*ny)-pos(j).y([1:ny]+i*ny)).^2];
    end
    for i = 1:ny-1
        d2 = [d2;(pos(j).x([1:nx:end]+(i-1))-pos(j).x([1:nx:end]+i)).^2+(pos(j).y([1:nx:end]+(i-1))-pos(j).y([1:nx:end]+i)).^2];
    end
    d2b(j) = mean(d2);
end
%}
