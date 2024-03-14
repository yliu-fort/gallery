clear; clc; close all;

%% Inviscid Burgers
t0 = 0;
tf = 10.0;
idim = 1001;
ntmaxi = 1250;
flux = zeros(idim+2,1);
u = zeros(idim,1);

dx = 1.0/idim;

% smooth profile
xpos = -1.0;

% Test 1
for i = 1:idim
    xpos = xpos + 2.0/idim;
    u(i) = exp(-8.0*xpos*xpos);
end

% Test 2
xleft = 0.1*1.0;
xmiddl = 0.5*1.0;
xright = 0.9*1.0;

for i = 1:idim
    xpos = (i - 1.0) * dx;
    if(xpos < xleft)
        u(i) = -1.0;
    end
    if(xpos >= xleft && xpos <= xmiddl)
        u(i) =  1.0;
    end
    if(xpos > xmiddl && xpos <= xright)
        u(i) =  0.0;
    end
    if(xpos > xright)
        u(i) = -1.0;
    end
end


%% Main program
t = t0;
sol = [];
for i = 1:ntmaxi
    u = BoundaryCorrection(u);
    dt = cflcon(0.9, u, dx, t, tf);
    t = t+dt;
    flux = fluxes(u);
    u = update(u, flux, dt, dx);
    sol = [sol u];
    
    if(t >= tf) break;end
end

%% Plot
surf(sol.','EdgeColor','none');
view(2), axis tight;
colormap gray;

%% Animate
h = plot(sol(:,1),'.-');
for i = 1:size(sol,2)
    h.YData = sol(:,i);
    drawnow
end

%% Helper functions
function u = BoundaryCorrection(u)
u(1) = u(end-1);
u(end)=u(2);

end

function dt = cflcon(cflcoe, u, dx, t, tf)
dt = cflcoe * dx / max(abs(u));

if((t+dt) > tf)
    dt = tf - t;
end

end

function u = update(u, flux, dt, dx)
u(2:end-1) = u(2:end-1) - dt/dx*diff(flux);

end

function flux = fluxes(u)
ul = u(1:end-1);
ur = u(2:end);

% Call riemann solver for inviscid burger
ustar = arrayfun(@riemann,ul,ur);

% Compute Godunov intercell flux
flux = 0.5*ustar.*ustar;

end

function ustar = riemann(ul,ur)
if(ul > ur)
    S = 0.5*(ul+ur);
    
    if(S >= 0.0)
        ustar = ul;
    else
        ustar = ur;
    end
else
    if(ul >= 0.0)
        ustar = ul;
    end
    
    if(ur <= 0.0)
        ustar = ur;
    end
    
    if(ul <= 0.0 && ur >= 0.0)
        ustar = 0.0;
    end
end

end