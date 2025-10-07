% run_hex_modal.m
% 8-node hexahedral (trilinear) solid element modal analysis
% Robust element orientation checking + full pipeline.
clear; close all; clc;

%% -------------------- Model parameters --------------------
E   = 210e9;      % Young's modulus (Pa)
nu  = 0.3;        % Poisson's ratio
rho = 7800;       % Density (kg/m^3)

% Geometry (meters)
a = 0.01;   % x-size
b = 0.01;   % y-size
L = 1.0;    % z-size (length)

% Mesh (number of elements in each direction)
nel_x = 1;
nel_y = 1;
nel_z = 20;

% Number of nodes in each direction
nnx = nel_x + 1;
nny = nel_y + 1;
nnz = nel_z + 1;

%% -------------------- Create nodes --------------------
xv = linspace(0, a, nnx);
yv = linspace(0, b, nny);
zv = linspace(0, L, nnz);

nNodes = nnx * nny * nnz;
nodes = zeros(nNodes, 3);

nodeID = @(ix,iy,iz) (iz-1)*nnx*nny + (iy-1)*nnx + ix;

for iz = 1:nnz
    for iy = 1:nny
        for ix = 1:nnx
            id = nodeID(ix,iy,iz);
            nodes(id,:) = [xv(ix), yv(iy), zv(iz)];
        end
    end
end

%% -------------------- Create connectivity (elems) --------------------
nElems = nel_x * nel_y * nel_z;
elems = zeros(nElems, 8);
e = 0;
for kz = 1:nel_z
    for jy = 1:nel_y
        for ix = 1:nel_x
            e = e + 1;
            % Standard positive-volume ordering (see comments below)
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

% Note on node ordering used above:
%   1:(-1,-1,-1), 2:(+1,-1,-1), 3:(+1,+1,-1), 4:(-1,+1,-1)
%   5:(-1,-1,+1), 6:(+1,-1,+1), 7:(+1,+1,+1), 8:(-1,+1,+1)
% This ordering usually yields positive Jacobian for regular meshes.

%% -------------------- Material matrix D --------------------
lambda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));
D = [lambda+2*mu, lambda,      lambda,      0,0,0;
     lambda,      lambda+2*mu, lambda,      0,0,0;
     lambda,      lambda,      lambda+2*mu, 0,0,0;
     0,0,0, mu,0,0;
     0,0,0, 0,mu,0;
     0,0,0, 0,0,mu];

%% -------------------- Gauss points & init global matrices --------------------
gp = [-1/sqrt(3), 1/sqrt(3)];
gw = [1, 1];

ndof = 3 * nNodes;
K = sparse(ndof, ndof);
M = sparse(ndof, ndof);

tol_det = 1e-12;  % threshold for acceptable detJ

%% -------------------- Element assembly with orientation check --------------------
fprintf('Starting element assembly (%d elements)...\n', nElems);

for ie = 1:nElems
    elemNodes = elems(ie,:);
    Xe = nodes(elemNodes, :); % 8x3

    % check duplicate nodes
    dup_found = false;
    for p=1:8
        for q=p+1:8
            if norm(Xe(p,:) - Xe(q,:)) < 1e-12
                fprintf('Element %d: duplicate nodes %d and %d (coords identical)\n', ie, p, q);
                dup_found = true;
            end
        end
    end
    if dup_found
        error('Mesh contains duplicate nodes in element %d. Durduruldu.', ie);
    end

    % Prepare candidate reorderings (try original first)
    candidateOrders = {};
    candidateNames  = {};
    candidateOrders{end+1} = elemNodes; candidateNames{end+1} = 'original';
    % Common fixes:
    ord = elemNodes; ord(5:8) = ord([8 7 6 5]); candidateOrders{end+1} = ord; candidateNames{end+1} = 'reverse top';
    ord = elemNodes; ord(1:4) = ord([4 3 2 1]); candidateOrders{end+1} = ord; candidateNames{end+1} = 'reverse bottom';
    ord = elemNodes([5 6 7 8 1 2 3 4]); candidateOrders{end+1} = ord; candidateNames{end+1} = 'swap top/bottom';
    ord = elemNodes([2 1 4 3 6 5 8 7]); candidateOrders{end+1} = ord; candidateNames{end+1} = 'reverse pairs';
    ord = elemNodes([8 7 6 5 4 3 2 1]); candidateOrders{end+1} = ord; candidateNames{end+1} = 'full reverse';
    ord = elemNodes; tmp=ord(2); ord(2)=ord(4); ord(4)=tmp; candidateOrders{end+1} = ord; candidateNames{end+1}='swap 2&4';
    ord = elemNodes; tmp=ord(5); ord(5)=ord(8); ord(8)=tmp; candidateOrders{end+1} = ord; candidateNames{end+1}='swap 5&8';

    chosenOrder = [];
    chosenName  = '';

    % Test each candidate ordering
    for c = 1:length(candidateOrders)
        Xe_try = nodes(candidateOrders{c}, :);
        minDet = inf;
        positive_all = true;
        for ii = 1:2
            for jj = 1:2
                for kk = 1:2
                    xi = gp(ii); eta = gp(jj); zeta = gp(kk);
                    dN_dxi = 1/8 * [
                        -(1-eta)*(1-zeta),  -(1-xi)*(1-zeta),  -(1-xi)*(1-eta);
                         (1-eta)*(1-zeta),   (1+xi)*(1-zeta),  -(1+xi)*(1-eta);
                         (1+eta)*(1-zeta),   (1+xi)*(1-zeta),  -(1+xi)*(1+eta);
                        -(1+eta)*(1-zeta),  -(1-xi)*(1-zeta),  -(1-xi)*(1+eta);
                        -(1-eta)*(1+zeta),  -(1-xi)*(1+zeta),   (1-xi)*(1-eta);
                         (1-eta)*(1+zeta),   (1+xi)*(1+zeta),   (1+xi)*(1-eta);
                         (1+eta)*(1+zeta),   (1+xi)*(1+zeta),   (1+xi)*(1+eta);
                        -(1+eta)*(1+zeta),  -(1-xi)*(1+zeta),   (1-xi)*(1+eta)
                    ];
                    J = dN_dxi' * Xe_try;
                    detJ = det(J);
                    if isnan(detJ)
                        detJ = -inf;
                    end
                    minDet = min(minDet, detJ);
                    if detJ <= tol_det
                        positive_all = false;
                        break;
                    end
                end
                if ~positive_all, break; end
            end
            if ~positive_all, break; end
        end
        if positive_all
            chosenOrder = candidateOrders{c};
            chosenName  = candidateNames{c};
            break;
        end
    end

    if isempty(chosenOrder)
        % none worked -> diagnostic and error
        % compute original detJ set for more info
        Xe_orig = Xe;
        dets = zeros(8,1); idx=0;
        for ii = 1:2
            for jj = 1:2
                for kk = 1:2
                    xi = gp(ii); eta = gp(jj); zeta = gp(kk);
                    dN_dxi = 1/8 * [
                        -(1-eta)*(1-zeta),  -(1-xi)*(1-zeta),  -(1-xi)*(1-eta);
                         (1-eta)*(1-zeta),   (1+xi)*(1-zeta),  -(1+xi)*(1-eta);
                         (1+eta)*(1-zeta),   (1+xi)*(1-zeta),  -(1+xi)*(1+eta);
                        -(1+eta)*(1-zeta),  -(1-xi)*(1-zeta),  -(1-xi)*(1+eta);
                        -(1-eta)*(1+zeta),  -(1-xi)*(1+zeta),   (1-xi)*(1-eta);
                         (1-eta)*(1+zeta),   (1+xi)*(1+zeta),   (1+xi)*(1-eta);
                         (1+eta)*(1+zeta),   (1+xi)*(1+zeta),   (1+xi)*(1+eta);
                        -(1+eta)*(1+zeta),  -(1-xi)*(1+zeta),   (1-xi)*(1+eta)
                    ];
                    J = dN_dxi' * Xe_orig;
                    idx = idx + 1;
                    dets(idx) = det(J);
                end
            end
        end
        fprintf('Element %d: orientation-fix failed. min detJ (original) = %.3e\n', ie, min(dets));
        fprintf('Node coords for element %d:\n', ie);
        for a=1:8
            fprintf('  %d: (%.6g, %.6g, %.6g)\n', a, Xe(a,1), Xe(a,2), Xe(a,3));
        end
        error('Element %d cannot be fixed automatically. Inspect mesh/connectivity.', ie);
    else
        if ~strcmp(chosenName, 'original')
            fprintf('Element %d: orientation fixed using -> %s\n', ie, chosenName);
            elems(ie,:) = chosenOrder; % update connectivity
            Xe = nodes(chosenOrder, :);
        end
    end

    % Now build Ke, Me with chosen Xe
    Ke = zeros(24,24);
    Me = zeros(24,24);
    for ii = 1:2
        for jj = 1:2
            for kk = 1:2
                xi = gp(ii); eta = gp(jj); zeta = gp(kk);
                w = gw(ii)*gw(jj)*gw(kk);

                N = 1/8 * [
                    (1-xi)*(1-eta)*(1-zeta);
                    (1+xi)*(1-eta)*(1-zeta);
                    (1+xi)*(1+eta)*(1-zeta);
                    (1-xi)*(1+eta)*(1-zeta);
                    (1-xi)*(1-eta)*(1+zeta);
                    (1+xi)*(1-eta)*(1+zeta);
                    (1+xi)*(1+eta)*(1+zeta);
                    (1-xi)*(1+eta)*(1+zeta)
                ];

                dN_dxi = 1/8 * [
                    -(1-eta)*(1-zeta),  -(1-xi)*(1-zeta),  -(1-xi)*(1-eta);
                     (1-eta)*(1-zeta),   (1+xi)*(1-zeta),  -(1+xi)*(1-eta);
                     (1+eta)*(1-zeta),   (1+xi)*(1-zeta),  -(1+xi)*(1+eta);
                    -(1+eta)*(1-zeta),  -(1-xi)*(1-zeta),  -(1-xi)*(1+eta);
                    -(1-eta)*(1+zeta),  -(1-xi)*(1+zeta),   (1-xi)*(1-eta);
                     (1-eta)*(1+zeta),   (1+xi)*(1+zeta),   (1+xi)*(1-eta);
                     (1+eta)*(1+zeta),   (1+xi)*(1+zeta),   (1+xi)*(1+eta);
                    -(1+eta)*(1+zeta),  -(1-xi)*(1+zeta),   (1-xi)*(1+eta)
                ];

                J = dN_dxi' * Xe;
                detJ = det(J);
                if detJ <= 0
                    error('Element %d: non-positive detJ encountered during assembly (%.3e).', ie, detJ);
                end
                invJ = inv(J);
                dN_dX = (invJ * dN_dxi')'; % 8x3

                B = zeros(6,24);
                NN = zeros(3,24);
                for a = 1:8
                    cols = (a-1)*3 + (1:3);
                    dNx = dN_dX(a,1); dNy = dN_dX(a,2); dNz = dN_dX(a,3);
                    % B matrix rows correspond to [eps_xx; eps_yy; eps_zz; eps_xy; eps_yz; eps_zx]
                    B(:, cols) = [ dNx,    0,    0;
                                   0,   dNy,    0;
                                   0,    0,   dNz;
                                   dNy,  dNx,   0;
                                   0,    dNz,  dNy;
                                   dNz,  0,    dNx ];
                    Ni = N(a);
                    NN(:, cols) = Ni * eye(3);
                end

                Ke = Ke + (B' * D * B) * detJ * w;
                Me = Me + rho * (NN' * NN) * detJ * w;
            end
        end
    end

    % Assemble to global K,M
    gdof = zeros(24,1);
    for a = 1:8
        gn = elems(ie,a);
        gdof((a-1)*3 + (1:3)) = (gn-1)*3 + (1:3);
    end
    K(gdof, gdof) = K(gdof, gdof) + Ke;
    M(gdof, gdof) = M(gdof, gdof) + Me;
end

fprintf('Assembly finished.\n');

%% -------------------- Apply boundary conditions (fixed-fixed at z=0 and z=L) --------------------
fixedNodes = [];
tol = 1e-9;
for n=1:nNodes
    if abs(nodes(n,3) - 0) < tol || abs(nodes(n,3) - L) < tol
        fixedNodes = [fixedNodes; n];
    end
end
fixedDofs = reshape([ (fixedNodes-1)*3 + 1, (fixedNodes-1)*3 + 2, (fixedNodes-1)*3 + 3 ]', [], 1);
freeDofs = setdiff(1:ndof, fixedDofs);

Kff = K(freeDofs, freeDofs);
Mff = M(freeDofs, freeDofs);

%% -------------------- Eigenvalue problem --------------------
nModes = 6;
fprintf('Solving eigenproblem for %d modes...\n', nModes);
try
    if size(Kff,1) > 300
        opts.tol = 1e-8; opts.maxit = 2000;
        [Phi, Omega2] = eigs(Kff, Mff, nModes, 'SM', opts);
    else
        % small problem -> dense solve more reliable
        [V,Dmat] = eig(full(Kff), full(Mff));
        [dvals, idxs] = sort(real(diag(Dmat)));
        Phi = V(:, idxs(1:nModes));
        Omega2 = diag(dvals(1:nModes));
    end
catch ME
    error('Eigenvalue solver failed: %s', ME.message);
end

omega = sqrt(abs(diag(Omega2)));
freqs = omega / (2*pi);

fprintf('FEM natural frequencies (Hz):\n');
for i=1:length(freqs)
    fprintf('%2d: %10.4f Hz\n', i, freqs(i));
end

%% -------------------- Theoretical 1D axial comparison (fixed-fixed) --------------------
c = sqrt(E / rho);
nvals = (1:6)';
f_theory = nvals * c / (2*L);

fprintf('\n1D axial (fixed-fixed) theoretical frequencies (Hz):\n');
for i=1:length(f_theory)
    fprintf('%2d: %10.4f Hz\n', i, f_theory(i));
end

fprintf('\nComparison (FEM vs Theory):\n');
fprintf('Mode\tFEM(Hz)\t\tTheory(Hz)\tRelErr(%%)\n');
for i=1:min(nModes, length(f_theory))
    relerr = 100*(freqs(i) - f_theory(i))/f_theory(i);
    fprintf('%d\t%8.3f\t%8.3f\t%8.3f\n', i, freqs(i), f_theory(i), relerr);
end

%% -------------------- Optional: plot first mode (axial displacement along centerline) --------------------
Phi_full = zeros(ndof, size(Phi,2));
Phi_full(freeDofs, :) = Phi;
centerNodes = find(abs(nodes(:,1)-a/2) < 1e-9 & abs(nodes(:,2)-b/2) < 1e-9);
if ~isempty(centerNodes)
    zvals = nodes(centerNodes, 3);
    uz = zeros(length(centerNodes),1);
    for i=1:length(centerNodes)
        n = centerNodes(i);
        uz(i) = Phi_full((n-1)*3 + 1, 1); % x-displacement of mode 1 as example
    end
    [zvals, ii] = sort(zvals);
    uz = uz(ii);
    figure; plot(zvals, uz, '-o'); xlabel('z (m)'); ylabel('u_x (mode 1)'); title('Mode 1 x-displacement (centerline)');
end

fprintf('\nScript finished successfully.\n');
