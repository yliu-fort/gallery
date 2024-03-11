function [r,u,v] = wake(t, NX, NY, rho0, nu, u_max)

r = zeros(NX, NY);
u = zeros(NX, NY);
v = zeros(NX, NY);

for j = 1:NY
    for i = 1:NX
        
        ky = 2.0*pi/NY;
        td = 1.0/(nu*(ky*ky));
        X = i+0.5;
        Y = j+0.5;
        width = 2;
        ux = u_max*tanh(100*exp(-1.0*t/td)*(width/2 - abs(width*(2*Y/NY-1))));
        uy = 0;
        
        rho = rho0;
        r(i,j) = rho;
        u(i,j) = ux;
        v(i,j) = uy;
    end
end

u = u(:);
v = v(:);
r = r(:);

end