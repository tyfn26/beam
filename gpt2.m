ncoord = 3;
ndof = 3;
nnode = 8;
coords = [0 0 0;
          1 0 0;
          1 1 0;
          0 1 0;
          0 0 1;
          1 0 1;
          1 1 1;
          0 1 1];
nelem = 1;
maxnodes = 8;
elident = 1;
lumpedmass = false;
nelnodes = 8;
connect = [1 2 3 4 5 6 7 8];
materialprops = [210e9, 0.3, 7800];

[Me, Ke] = Me_Ke_Hex8(ncoord, ndof, nnode, coords, nelem, maxnodes, ...
                      elident, lumpedmass, nelnodes, connect, materialprops);

fprintf('Total mass = %.6f\n', sum(Me(:))/3);
function [Me, Ke] = Me_Ke_Hex8(ncoord, ndof, nnode, coords, nelem, maxnodes, ...
                                elident, lumpedmass, nelnodes, connect, materialprops)
% ============================================================
% Computes global stiffness and mass matrices for HEX8 elements
% Inputs:
%   ncoord        : number of spatial coordinates (3)
%   ndof          : degrees of freedom per node (3)
%   nnode         : total number of nodes
%   coords(nnode,3): nodal coordinates
%   nelem         : number of elements
%   maxnodes      : max nodes per element (8)
%   elident       : element type identifier (not used here)
%   lumpedmass    : true for lumped mass, false for consistent
%   nelnodes      : number of nodes per element (8)
%   connect(nelem,8): element connectivity
%   materialprops : [E, nu, rho]
%
% Outputs:
%   Ke : global stiffness matrix
%   Me : global mass matrix
% ============================================================

E   = materialprops(1);
nu  = materialprops(2);
rho = materialprops(3);

lambda = E*nu/((1+nu)*(1-2*nu));
mu     = E/(2*(1+nu));
C = [lambda+2*mu, lambda,     lambda,     0,0,0;
     lambda,     lambda+2*mu, lambda,     0,0,0;
     lambda,     lambda,     lambda+2*mu, 0,0,0;
     0,0,0, mu,0,0;
     0,0,0, 0,mu,0;
     0,0,0, 0,0,mu];

% Global matrix initialization
neq = nnode * ndof;
Ke = zeros(neq, neq);
Me = zeros(neq, neq);

% Gauss integration points
ngp = 2;
[gp, w] = gauss_1d(ngp);
node_loc = [-1 -1 -1;
             1 -1 -1;
             1  1 -1;
            -1  1 -1;
            -1 -1  1;
             1 -1  1;
             1  1  1;
            -1  1  1];

% === Loop over elements ===
for e = 1:nelem
    con = connect(e,1:nelnodes);
    Xe = coords(con,:);
    Ke_e = zeros(ndof*nelnodes);
    Me_e = zeros(ndof*nelnodes);

    for i = 1:ngp
        xi = gp(i); wi = w(i);
        for j = 1:ngp
            eta = gp(j); wj = w(j);
            for k = 1:ngp
                zeta = gp(k); wk = w(k);
                wt = wi*wj*wk;

                [N, dN_dxi] = hex8_shape(xi, eta, zeta, node_loc);
                J = dN_dxi * Xe;
                detJ = det(J);
                if detJ <= 0
                    error('Negative or zero det(J) in element %d', e);
                end
                invJ = inv(J);
                dN_dx = invJ * dN_dxi;

                % B-matrix
                B = zeros(6, ndof*nelnodes);
                for a = 1:nelnodes
                    idx = (a-1)*ndof + (1:ndof);
                    dNa = dN_dx(:,a);
                    B(:, idx) = [ dNa(1)      0           0;
                                   0        dNa(2)        0;
                                   0          0        dNa(3);
                                   dNa(2)    dNa(1)      0;
                                   0        dNa(3)    dNa(2);
                                   dNa(3)      0      dNa(1) ];
                end

                % N-matrix (for mass)
                Nmat = zeros(3, ndof*nelnodes);
                for a = 1:nelnodes
                    idx = (a-1)*3 + (1:3);
                    Nmat(1, idx(1)) = N(a);
                    Nmat(2, idx(2)) = N(a);
                    Nmat(3, idx(3)) = N(a);
                end

                Ke_e = Ke_e + (B' * C * B) * detJ * wt;
                Me_e = Me_e + rho * (Nmat' * Nmat) * detJ * wt;
            end
        end
    end

    % === Optional lumped mass ===
    if lumpedmass
        Me_e = diag(sum(Me_e,2));
    end

    % === Assembly into global ===
    dofmap = zeros(1, nelnodes*ndof);
    for a = 1:nelnodes
        dofmap((a-1)*ndof + (1:ndof)) = (con(a)-1)*ndof + (1:ndof);
    end

    Ke(dofmap, dofmap) = Ke(dofmap, dofmap) + Ke_e;
    Me(dofmap, dofmap) = Me(dofmap, dofmap) + Me_e;
end

end

% === Supporting Functions ===

function [gp, w] = gauss_1d(n)
if n==2
    gp = [-1 1]/sqrt(3);
    w  = [1 1];
elseif n==3
    gp = [-sqrt(3/5) 0 sqrt(3/5)];
    w  = [5/9 8/9 5/9];
else
    error('Only 2 or 3 Gauss points supported.');
end
end

function [N, dN_dxi] = hex8_shape(xi, eta, zeta, node_loc)
N = 1/8 * (1 + node_loc(:,1)*xi) .* (1 + node_loc(:,2)*eta) .* (1 + node_loc(:,3)*zeta);
dN_dxi = zeros(3,8);
for i = 1:8
    dN_dxi(1,i) = 1/8 * node_loc(i,1) * (1 + node_loc(i,2)*eta) * (1 + node_loc(i,3)*zeta);
    dN_dxi(2,i) = 1/8 * (1 + node_loc(i,1)*xi) * node_loc(i,2) * (1 + node_loc(i,3)*zeta);
    dN_dxi(3,i) = 1/8 * (1 + node_loc(i,1)*xi) * (1 + node_loc(i,2)*eta) * node_loc(i,3);
end
end
