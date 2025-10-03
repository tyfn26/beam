% 8-node hexahedral (trilinear) solid element modal analysis (MATLAB)
% Kullanım: dosyayı kaydedip çalıştırın. Parametreleri başta değiştirebilirsiniz.
clear; close all; clc;

%% Model parametreleri (değiştirilebilir)
E = 210e9;         % Young's modulus (Pa)
nu = 0.3;          % Poisson's ratio
rho = 7800;        % density (kg/m^3)

L = 1.0;           % Length in z-direction (m)
a = 0.01; b = 0.01; % cross-section (x and y) (m) - ince çubuk

nel_x = 1; nel_y = 1; nel_z = 20; % mesh: 1x1x20 eleman (uzunluk boyunca bölünmüş)
% (uzunluk boyunca daha çok eleman axial konverjans için faydalı)

% derived
nnx = nel_x+1; nny = nel_y+1; nnz = nel_z+1;
nNodes = nnx*nny*nnz;
ndof = 3*nNodes;

% Generate nodal coordinates
xv = linspace(0,a,nnx);
yv = linspace(0,b,nny);
zv = linspace(0,L,nnz);

nodes = zeros(nNodes,3);
nodeID = @(ix,iy,iz) (iz-1)*nnx*nny + (iy-1)*nnx + ix;
for iz=1:nnz
    for iy=1:nny
        for ix=1:nnx
            id = nodeID(ix,iy,iz);
            nodes(id,:) = [xv(ix), yv(iy), zv(iz)];
        end
    end
end

% Connectivity (8-node hexahedron, node ordering consistent)
elems = zeros(nel_x*nel_y*nel_z,8);
e = 0;
for kz=1:nel_z
    for jy=1:nel_y
        for ix=1:nel_x
            e = e+1;
            n1 = nodeID(ix,   jy,   kz);
            n2 = nodeID(ix+1, jy,   kz);
            n3 = nodeID(ix+1, jy+1, kz);
            n4 = nodeID(ix,   jy+1, kz);
            n5 = nodeID(ix,   jy,   kz+1);
            n6 = nodeID(ix+1, jy,   kz+1);
            n7 = nodeID(ix+1, jy+1, kz+1);
            n8 = nodeID(ix,   jy+1, kz+1);
            elems(e,:) = [n1 n2 n3 n4 n5 n6 n7 n8];
        end
    end
end
nElems = size(elems,1);

%% Material matrix (isotropic linear elasticity, 3D)
lambda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));
D = [lambda+2*mu, lambda,      lambda,      0,0,0;
     lambda,      lambda+2*mu, lambda,      0,0,0;
     lambda,      lambda,      lambda+2*mu, 0,0,0;
     0,0,0, mu,0,0;
     0,0,0, 0,mu,0;
     0,0,0, 0,0,mu];

%% Gauss points for 2x2x2
gp = [-1/sqrt(3), 1/sqrt(3)];
gw = [1, 1];

% Initialize global K and M
K = sparse(ndof, ndof);
M = sparse(ndof, ndof);

% Loop elements: compute K_e and M_e
for ie=1:nElems
    elemNodes = elems(ie,:);
    Xe = nodes(elemNodes, :); % 8x3
    Ke = zeros(24,24);
    Me = zeros(24,24);
    
    % 2x2x2 Gauss integration
    for i=1:2
        for j=1:2
            for k=1:2
                xi = gp(i); eta = gp(j); zeta = gp(k);
                w = gw(i)*gw(j)*gw(k);
                % Shape functions N (8) and derivatives dN/dxi
                N = 1/8*[ (1-xi)*(1-eta)*(1-zeta);
                         (1+xi)*(1-eta)*(1-zeta);
                         (1+xi)*(1+eta)*(1-zeta);
                         (1-xi)*(1+eta)*(1-zeta);
                         (1-xi)*(1-eta)*(1+zeta);
                         (1+xi)*(1-eta)*(1+zeta);
                         (1+xi)*(1+eta)*(1+zeta);
                         (1-xi)*(1+eta)*(1+zeta)];
                dN_dxi = 1/8*[
                    [ -(1-eta)*(1-zeta), -(1-xi)*(1-zeta), -(1-xi)*(1-eta) ];
                    [  (1-eta)*(1-zeta),  (1+xi)*(1-zeta), -(1+xi)*(1-eta) ];
                    [  (1+eta)*(1-zeta),  (1+xi)*(1-zeta), -(1+xi)*(1+eta) ];
                    [ -(1+eta)*(1-zeta), -(1-xi)*(1-zeta), -(1-xi)*(1+eta) ];
                    [ -(1-eta)*(1+zeta), -(1-xi)*(1+zeta),  (1-xi)*(1-eta) ];
                    [  (1-eta)*(1+zeta),  (1+xi)*(1+zeta),  (1+xi)*(1-eta) ];
                    [  (1+eta)*(1+zeta),  (1+xi)*(1+zeta),  (1+xi)*(1+eta) ];
                    [ -(1+eta)*(1+zeta), -(1-xi)*(1+zeta),  (1-xi)*(1+eta) ];
                    ];
                % Jacobian
                J = dN_dxi' * Xe; % 3x3
                detJ = det(J);
                if detJ <= 0
                    error('Non-positive Jacobian det at element %d', ie);
                end
                invJ = inv(J);
                % derivatives dN/dX
                dN_dX = (invJ * dN_dxi')'; % 8x3
                % B matrix (6 x 24)
                B = zeros(6, 24);
                for a=1:8
                    idx = (a-1)*3 + (1:3);
                    dNx = dN_dX(a,1); dNy = dN_dX(a,2); dNz = dN_dX(a,3);
                    B(:, idx) = [ dNx,    0,    0;
                                   0,   dNy,    0;
                                   0,    0,   dNz;
                                   dNy,  dNx,   0;
                                   0,    dNz,  dNy;
                                   dNz,  0,    dNx ]';
                end
                Ke = Ke + (B' * D * B) * detJ * w;
                % Consistent mass: Me += rho * (N^T * N) * detJ * w
                NN = zeros(3,24);
                for a=1:8
                    Ni = N(a);
                    NN(:, (a-1)*3 + (1:3)) = Ni * eye(3);
                end
                Me = Me + rho * (NN' * NN) * detJ * w;
            end
        end
    end
    
    % Assemble into global
    gdof = zeros(24,1);
    for a=1:8
        gn = elemNodes(a);
        gdof((a-1)*3 + (1:3)) = (gn-1)*3 + (1:3);
    end
    K(gdof, gdof) = K(gdof, gdof) + Ke;
    M(gdof, gdof) = M(gdof, gdof) + Me;
end

%% Apply boundary conditions - FIXED-FIXED along z ends (all DOFs fixed at z=0 and z=L)
fixedNodes = [];
tol = 1e-9;
for n=1:nNodes
    if abs(nodes(n,3) - 0) < tol || abs(nodes(n,3) - L) < tol
        fixedNodes = [fixedNodes; n];
    end
end
fixedDofs = reshape([ (fixedNodes-1)*3 + 1, (fixedNodes-1)*3 + 2, (fixedNodes-1)*3 + 3 ]',[],1);

freeDofs = setdiff(1:ndof, fixedDofs);

Kff = K(freeDofs, freeDofs);
Mff = M(freeDofs, freeDofs);

%% Solve generalized eigenvalue problem
nModes = 6; % kaç mod isteniyor
opts.maxit = 500;
opts.tol = 1e-8;
% use eigs
[Phi, Omega2] = eigs(Kff, Mff, nModes, 'SM'); % SM: smallest magnitude eigenvalues (lowest freq)
omega = sqrt(abs(diag(Omega2))); % rad/s
freqs = omega / (2*pi);

%% Display modal frequencies
fprintf('FEM doğal frekanslar (Hz):\n');
for i=1:length(freqs)
    fprintf('%2d: %10.4f Hz\n', i, freqs(i));
end

%% Theoretical axial (1D rod) frequencies (fixed-fixed): f_n = n * c / (2L), c = sqrt(E/rho)
c = sqrt(E/rho);
nvals = (1:6)';
f_theory = nvals * c / (2*L);
fprintf('\n1D axial (fixed-fixed) teorik frekanslar (Hz):\n');
for i=1:length(nvals)
    fprintf('%2d: %10.4f Hz\n', nvals(i), f_theory(i));
end

%% Compare: show table (FEM vs Theory)
fprintf('\nKarşılaştırma (FEM vs 1D axial teorik):\n');
fprintf('Mode\tFEM(Hz)\t\tTheory(Hz)\tRelError(%%)\n');
for i=1:min(nModes, length(f_theory))
    relerr = 100*(freqs(i) - f_theory(i))/f_theory(i);
    fprintf('%d\t%8.3f\t%8.3f\t%8.3f\n', i, freqs(i), f_theory(i), relerr);
end

%% (Opsiyonel) Plot first mode shape (x displacement) as an example
% expand Phi to full dofs
Phi_full = zeros(ndof, size(Phi,2));
Phi_full(freeDofs, :) = Phi;
% plot z vs uz at centerline (x=a/2, y=b/2)
centerNodes = find(abs(nodes(:,1)-a/2)<1e-6 & abs(nodes(:,2)-b/2)<1e-6);
if ~isempty(centerNodes)
    uz = zeros(length(centerNodes),1);
    zvals = zeros(length(centerNodes),1);
    for i=1:length(centerNodes)
        n = centerNodes(i);
        zvals(i) = nodes(n,3);
        uz(i) = Phi_full((n-1)*3 + 1, 1); % x displacement of mode 1 (example)
    end
    [zvals, idxs] = sort(zvals);
    uz = uz(idxs);
    figure; plot(zvals, uz, '-o'); xlabel('z (m)'); ylabel('u_x (mode 1)'); title('Mode 1 x-displacement (centerline)');
end
