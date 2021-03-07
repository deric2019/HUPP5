function nmat=nmat_GRIN(n_max,n_clad,D_core,xmat,ymat)

% n_max: refractive index at core center
% n_clad: refractive index in cladding
% D_core: core diameter (i.e. diameter of graded index region)
% xmat, ymat: coordinate matrices

rmat=sqrt(xmat.^2+ymat.^2);
const=(n_max-n_clad)/(D_core/2)^2;
nmat=(rmat<(D_core/2)).*(n_max-const*rmat.^2)+(rmat>=(D_core/2))*n_clad;

figure(50)
mesh(xmat(1,:)*1e6,ymat(:,1)*1e6,nmat)
title('refractive index distribution')
xlabel('[µm]')
ylabel('[µm]')
colormap(jet)
drawnow