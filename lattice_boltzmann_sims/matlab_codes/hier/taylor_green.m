function [r,u,v] = taylor_green(t, NX, NY, rho0, nu, u_max)

r = zeros(NX, NY);
u = zeros(NX, NY);
v = zeros(NX, NY);

for j = 1:NY
    for i = 1:NX
        
        kx = 2.0*pi/NX;
        ky = 2.0*pi/NY;
        td = 1.0/(nu*(kx*kx+ky*ky));
        X = i+0.5;
        Y = j+0.5;
        ux = -u_max*sqrt(ky/kx)*cos(kx*X)*sin(ky*Y)*exp(-1.0*t/td);
        uy = u_max*sqrt(kx/ky)*sin(kx*X)*cos(ky*Y)*exp(-1.0*t/td);
        
        P = -0.25*rho0*u_max*u_max*( (ky/kx)*cos(2.0*kx*X) ...
            +(kx/ky)*cos(2.0*ky*Y) )*exp(-2.0*t/td);
        
        rho = rho0+3.0*P;
        r(i,j) = rho;
        u(i,j) = ux;
        v(i,j) = uy;
    end
end

u = u(:);
v = v(:);
r = r(:);

end