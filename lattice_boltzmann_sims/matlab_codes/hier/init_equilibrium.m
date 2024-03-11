function buf = init_equilibrium(rho, ux, uy, wi, dirx, diry, NX, NY, NDIR)

buf = zeros(NX*NY, NDIR, 2);

f=compute_equilibrium(rho, ux, uy, dirx, diry, wi);

buf(:,:,1) = f;
end