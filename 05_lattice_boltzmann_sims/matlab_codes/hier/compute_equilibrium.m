function feq = compute_equilibrium(rho, ux, uy, dirx, diry, wi)
feq=ux*dirx+uy*diry;
feq=(3+4.5*feq).*feq;
feq=bsxfun(@minus,feq,1.5*(ux.^2+uy.^2));
feq=bsxfun(@times,1.0+feq,wi);
feq=bsxfun(@times,feq,rho);
end