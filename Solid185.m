function [Ke, Me] = Hex8_Solid185_Assembly(ncoord, ndof, nnode, coords, nelem, maxnodes, ...
                                elident, lumpedmass, nelnodes, connect, materialprops)
% Hex8_Solid185_Assembly
% Advanced Hex8 element assembly with:
% - Selective Reduced Integration (SRI): volumetric 1x1x1, deviatoric 2x2x2
% - B-bar (assumed volumetric strain)
% - Simple incompatible-mode (enhanced strain) approximation
% - Hourglass stabilization
% - Triplet sparse assembly (memory efficient)
%
% Inputs:
%   coords (nnode x 3), connect (nelem x 8), materialprops = [E, nu, rho]
%   lumpedmass: true => lumped mass (row-sum), false => consistent
%
% Outputs:
%   Ke, Me : sparse global stiffness and mass matrices (neq x neq)
%
% Notes:
% - ndof must be 3
% - nelnodes must be 8
% - This code focuses on linear elasticity modal analysis

% --- params
E   = materialprops(1);
nu  = materialprops(2);
rho = materialprops(3);

if ndof~=3 || nelnodes~=8
    error('This routine assumes ndof=3 and nelnodes=8');
end

% Lame
lambda = E*nu/((1+nu)*(1-2*nu));
mu     = E/(2*(1+nu));

% Constitutive (Voigt) - not used directly for split (we build dev/vol parts)
C = [lambda+2*mu, lambda,     lambda,     0,0,0;
     lambda,     lambda+2*mu, lambda,     0,0,0;
     lambda,     lambda,     lambda+2*mu, 0,0,0;
     0,0,0, mu,0,0;
     0,0,0, 0,mu,0;
     0,0,0, 0,0,mu];

% Gauss points/types
% deviatoric: 2x2x2
[g2,w2] = gauss_1d(2);
% volumetric: 1x1x1 (center)
[g1,w1] = gauss_1d(1);

% node_loc standard
node_loc = [-1 -1 -1;
             1 -1 -1;
             1  1 -1;
            -1  1 -1;
            -1 -1  1;
             1 -1  1;
             1  1  1;
            -1  1  1];

% Prepare triplets
I_K=[]; J_K=[]; V_K=[];
I_M=[]; J_M=[]; V_M=[];

% Hourglass parameters
hg_gamma = 0.01; % coefficient (tunable), small

% Incompatible-mode: define small set of bubble modes (3 modes) for approximation
% We'll construct matrix G_im which maps incompatible DOF -> element strains
% For simplicity use three linear incompatible shape functions that vanish on boundary
% The implementation approximates effect as additional stiffness Ke_im = alpha * (G'*G)
alpha_im = mu * 1e-2; % scale factor for incompatible modes (tunable)

% Loop elements
for e=1:nelem
    con = connect(e,1:nelnodes);
    Xe = coords(con, :); % 8x3
    
    % element DOF map
    dofmap = zeros(1, nelnodes*ndof);
    for a=1:nelnodes
        dofmap((a-1)*ndof + (1:ndof)) = (con(a)-1)*ndof + (1:ndof);
    end
    
    % --- evaluate Bbar: compute volumetric strain shape (1x8) averaged at element center
    % At center xi=0,eta=0,zeta=0
    [N_c, dN_dxi_c] = hex8_shape(0,0,0, node_loc);
    Jc = dN_dxi_c * Xe;
    detJc = det(Jc);
    if detJc <= 0
        error('Element %d has non-positive detJ at center', e);
    end
    invJc = inv(Jc);
    dN_dx_c = invJc * dN_dxi_c; % 3x8
    % volumetric Bbar vector at center: sum_i dN_dx_c(:,i) * N weighting?
    % For Bbar we need averaged dN_dx for volumetric part — compute weights for 1x1x1 GP:
    % For hexahedra center only, Bvol = [dN_dx_c sums]
    % compute scalar shape sum per direction
    % Construct Bvol_bar (1x24) mapping to volumetric strain (trace)
    Bv_bar = zeros(1, nelnodes*ndof);
    for a=1:nelnodes
        idx = (a-1)*3 + (1:3);
        % volumetric strain = dNx/dx + dNy/dy + dNz/dz
        Bv_bar(idx) = [dN_dx_c(1,a), dN_dx_c(2,a), dN_dx_c(3,a)];
    end
    % We'll use this Bv_bar to compute volumetric contribution as Bv_bar'*Bv_bar scaled
    
    % --- Initialize element stiffness and mass (dense small matrices)
    nd = nelnodes*ndof;
    Ke_e = zeros(nd, nd);
    Me_e = zeros(nd, nd);
    
    % --- Deviatoric part: integrate over 2x2x2 gauss points (shear)
    for i=1:2
        xi = g2(i); wi = w2(i);
        for j=1:2
            eta = g2(j); wj = w2(j);
            for k=1:2
                zeta = g2(k); wk = w2(k);
                wt = wi*wj*wk;
                [N, dN_dxi] = hex8_shape(xi, eta, zeta, node_loc);
                J = dN_dxi * Xe;
                detJ = det(J);
                if detJ <= 0
                    error('Element %d has non-positive detJ at GP (%g,%g,%g)', e, xi,eta,zeta);
                end
                invJ = inv(J);
                dN_dx = invJ * dN_dxi; % 3x8
                
                % Build full B (6x24)
                B = zeros(6, nd);
                for a=1:nelnodes
                    idx = (a-1)*3 + (1:3);
                    dNa = dN_dx(:,a);
                    B(:, idx) = [ dNa(1)      0           0;
                                   0        dNa(2)        0;
                                   0          0        dNa(3);
                                   dNa(2)    dNa(1)      0;
                                   0        dNa(3)    dNa(2);
                                   dNa(3)      0      dNa(1) ];
                end
                
                % Split into volumetric and deviatoric components:
                % volumetric contribution uses Bv_bar (from center) not local Bv
                % compute local volumetric part BV_local (1x24): trace(B) mapping
                BV_local = zeros(1,nd); % dNx/dx + dNy/dy + dNz/dz
                for a=1:nelnodes
                    idx = (a-1)*3 + (1:3);
                    BV_local(idx) = [dN_dx(1,a), dN_dx(2,a), dN_dx(3,a)];
                end
                
                % deviatoric B matrix: remove volumetric part from B (Voigt)
                % Build projection: For strain vector e = [ex;ey;ez;gxy;gyz;gzx]
                % volumetric strain = (ex+ey+ez)/3 ; we remove that part
                % Form operator P_vol that maps BV_local to 6x24 contribution
                % Construct B_vol_matrix (6xnd) that corresponds to volumetric part:
                Bvol = zeros(6, nd);
                for a=1:nd
                    % BV_local is 1xnd, but volumetric contribution to ex,ey,ez equal to BV_local/3
                    % as columns mapping we set
                end
                % easier approach: form full constitutive split:
                % Use full C, but subtract volumetric part using Bv_bar later (SRI will handle)
                
                % Compute element deviatoric stiffness contribution using C_dev
                % Build C_dev = C - C_vol where C_vol corresponds to pure volumetric part:
                % volumetric part: C_vol = K_bulk * (1/3) * Jv, where Jv is matrix of ones in 3x3 block
                Kbulk = lambda + 2*mu/3; % bulk modulus approximation (for split)
                % Build 6x6 volumetric projector matrix Pv
                Pv = zeros(6,6);
                Pv(1,1)=1; Pv(1,2)=1; Pv(1,3)=1;
                Pv(2,1)=1; Pv(2,2)=1; Pv(2,3)=1;
                Pv(3,1)=1; Pv(3,2)=1; Pv(3,3)=1;
                Pv = Pv/3;
                C_vol = Kbulk * Pv;
                C_dev = C - C_vol;
                
                Ke_e = Ke_e + (B' * C_dev * B) * detJ * wt;
                
                % Consistent mass contribution (use Nmat)
                Nmat = zeros(3, nd);
                for a=1:nelnodes
                    idx = (a-1)*3 + (1:3);
                    Nmat(1, idx(1)) = N(a);
                    Nmat(2, idx(2)) = N(a);
                    Nmat(3, idx(3)) = N(a);
                end
                Me_e = Me_e + rho * (Nmat' * Nmat) * detJ * wt;
            end
        end
    end
    
    % --- Volumetric part (1x1x1 integration) using Bbar
    % Bv_bar (1xnd) computed earlier; form 6xnd volumetric Bbar matrix:
    Bbar_vol = zeros(6, nd);
    % fill first 3 rows (ex,ey,ez)
    for a=1:nd
        % Bv_bar holds [dNx,dNy,dNz] interleaved per DOF
        Bbar_vol(1,a) = Bv_bar(a) / 3; % ex gets one third of volumetric contribution per DOF
        Bbar_vol(2,a) = Bv_bar(a) / 3;
        Bbar_vol(3,a) = Bv_bar(a) / 3;
    end
    % volumetric stiffness
    % C_vol computed above
    % integrate at center
    Ke_e = Ke_e + (Bbar_vol' * C_vol * Bbar_vol) * detJc * w1;
    
    % --- Incompatible-mode (enhanced strain) approx.
    % Build simple incompatible shape G (nd x nmodes)
    nmodes = 3;
    G = zeros(nd, nmodes);
    % Define simple bubble-like modes (zero on faces approx):
    % mode1: linear in xi* (1 - zeta^2) mapped approximately via nodal positions
    % we'll use small heuristic: for node a, mode values = (xi_a + eta_a + zeta_a)/8 centered
    for a=1:nelnodes
        xi_a = node_loc(a,1); eta_a = node_loc(a,2); zeta_a = node_loc(a,3);
        gv1 = xi_a*(1 - eta_a^2)*(1 - zeta_a^2);
        gv2 = eta_a*(1 - xi_a^2)*(1 - zeta_a^2);
        gv3 = zeta_a*(1 - xi_a^2)*(1 - eta_a^2);
        idx = (a-1)*3 + (1:3);
        % distribute to three directional DOFs roughly
        G(idx,1) = gv1;      % mode1 in x-DOF pattern
        G(idx,2) = gv2;      % mode2 in y-DOF pattern
        G(idx,3) = gv3;      % mode3 in z-DOF pattern
    end
    % Form small incompatible stiffness approx: Ke_im = alpha_im * (G * G')
    Ke_im = alpha_im * (G * G');
    Ke_e = Ke_e + Ke_im;
    
    % --- Hourglass stabilization (simple form)
    % compute hourglass vectors q (nd x ngh), here ngh=4 typical for 8-node,
    % but implement simple 4-mode hourglass (one per face diagonal)
    hg_vecs = hourglass_vectors(node_loc); % returns nd x ngh
    ngh = size(hg_vecs,2);
    for h=1:ngh
        qh = hg_vecs(:,h);
        % compute hourglass stiffness contribution k_h = gamma * mu * (qh * qh') * (volume)
        Ke_e = Ke_e + hg_gamma * mu * (qh * qh') * detJc;
    end
    
    % --- Lump if requested
    if lumpedmass
        Me_e = diag(sum(Me_e,2));
    end
    
    % --- Append element triplets
    [ii,jj,vv] = find(Ke_e);
    I_K = [I_K; dofmap(ii)']; J_K = [J_K; dofmap(jj)']; V_K = [V_K; vv];
    [ii,jj,vv] = find(Me_e);
    I_M = [I_M; dofmap(ii)']; J_M = [J_M; dofmap(jj)']; V_M = [V_M; vv];
end

% Build sparse global matrices
neq = nnode * ndof;
Ke = sparse(I_K, J_K, V_K, neq, neq);
Me = sparse(I_M, J_M, V_M, neq, neq);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gp,w] = gauss_1d(n)
if n==1
    gp = 0; w = 2;
elseif n==2
    gp = [-1 1]/sqrt(3); w = [1 1];
elseif n==3
    gp = [-sqrt(3/5) 0 sqrt(3/5)]; w = [5/9 8/9 5/9];
else
    error('gauss_1d supports n=1,2,3');
end
end

function [N, dN_dxi] = hex8_shape(xi, eta, zeta, node_loc)
% node_loc: 8x3 with local +-1 coordinates as used earlier
N = zeros(1,8);
dN_dxi = zeros(3,8);
for a=1:8
    xi_a = node_loc(a,1); eta_a = node_loc(a,2); zeta_a = node_loc(a,3);
    N(a) = 1/8*(1 + xi*xi_a)*(1 + eta*eta_a)*(1 + zeta*zeta_a);
    dN_dxi(1,a) = 1/8 * xi_a * (1 + eta*eta_a) * (1 + zeta*zeta_a);
    dN_dxi(2,a) = 1/8 * (1 + xi*xi_a) * eta_a * (1 + zeta*zeta_a);
    dN_dxi(3,a) = 1/8 * (1 + xi*xi_a) * (1 + eta*eta_a) * zeta_a;
end
end

function hg = hourglass_vectors(node_loc)
% Produce simple hourglass vectors (nd x ngh)
% node_loc is 8x3 local coords. We'll create 4 basic hourglass modes:
% vector values alternate sign to produce zero integral over faces.
nd = 8*3;
hg = zeros(nd,4);
% For convenience compute per-node scalar functions and expand to dofs
% mode patterns inspired by typical hourglass basis
for a=1:8
    xi = node_loc(a,1); eta = node_loc(a,2); zeta = node_loc(a,3);
    % 4 patterns:
    s1 = xi*eta;
    s2 = xi*zeta;
    s3 = eta*zeta;
    s4 = xi*eta*zeta;
    idx = (a-1)*3 + (1:3);
    % distribute scalar to components with sign patterns
    hg(idx,1) = s1 * [1;1;1];
    hg(idx,2) = s2 * [1;1;1];
    hg(idx,3) = s3 * [1;1;1];
    hg(idx,4) = s4 * [1;1;1];
end
% normalize each vector
for k=1:4
    v = hg(:,k);
    if norm(v) > 0
        hg(:,k) = v / norm(v);
    end
end
end
