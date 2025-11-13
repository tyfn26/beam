function [Ke, Me] = hex8_element(X, E, nu, rho)
% HEX8_ELEMENT  returns element stiffness and consistent mass (24x24 each)
% X: 8x3 nodal coordinates [x y z] in node order:
%   1:(-1,-1,-1), 2:(1,-1,-1), 3:(1,1,-1), 4:(-1,1,-1),
%   5:(-1,-1,1),  6:(1,-1,1),  7:(1,1,1),  8:(-1,1,1)
% E: Young's modulus
% nu: Poisson's ratio
% rho: density
%
% Usage: [Ke, Me] = hex8_element(X, E, nu, rho)

% --- Prepare
ngp = 2; % 2x2x2 Gauss
[g,w] = gauss_1d(ngp);

Ke = zeros(24,24);
Me = zeros(24,24);

% Elasticity matrix C (6x6 Voigt)
lambda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));
C = [lambda+2*mu, lambda,     lambda,     0, 0, 0;
     lambda,     lambda+2*mu, lambda,     0, 0, 0;
     lambda,     lambda,     lambda+2*mu, 0, 0, 0;
     0,0,0, mu, 0,0;
     0,0,0, 0, mu,0;
     0,0,0, 0,0, mu];

% Node local coords xi_i,eta_i,zeta_i
node_loc = [-1 -1 -1;
             1 -1 -1;
             1  1 -1;
            -1  1 -1;
            -1 -1  1;
             1 -1  1;
             1  1  1;
            -1  1  1];

% Loop Gauss points
for i=1:ngp
  xi = g(i);
  wi = w(i);
  for j=1:ngp
    eta = g(j);
    wj = w(j);
    for k=1:ngp
      zeta = g(k);
      wk = w(k);
      wt = wi*wj*wk;

      % Shape functions N (1x8) and derivatives dN/dxi (3x8)
      [N, dN_dxi] = hex8_shape(xi, eta, zeta, node_loc);

      % Jacobian (3x3)
      J = dN_dxi * X;        % dN_dxi (3x8) * X (8x3) -> 3x3
      detJ = det(J);
      if detJ <= 0
        error('Non-positive Jacobian determinant: detJ=%.4e', detJ)
      end
      invJ = inv(J);

      % dN/dx (3x8)
      dN_dx = invJ * dN_dxi;  % 3x8

      % Build B (6x24)
      B = zeros(6,24);
      for a = 1:8
        idx = (a-1)*3 + (1:3);
        dNa = dN_dx(:,a);  % [dNx; dNy; dNz]
        B(:, idx) = [ dNa(1)      0           0;
                       0        dNa(2)        0;
                       0          0        dNa(3);
                       dNa(2)    dNa(1)      0;
                       0        dNa(3)    dNa(2);
                       dNa(3)      0      dNa(1) ];
      end

      % N matrix for vector field (3x24)
      Nmat = zeros(3,24);
      for a=1:8
        idx = (a-1)*3 + (1:3);
        Nmat(1, idx(1)) = N(a);
        Nmat(2, idx(2)) = N(a);
        Nmat(3, idx(3)) = N(a);
      end

      % Integrate
      Ke = Ke + (B' * C * B) * detJ * wt;
      Me = Me + (rho * (Nmat' * Nmat)) * detJ * wt;
    end
  end
end

end

% -----------------------
% Helper functions below
% -----------------------
function [g, w] = gauss_1d(ngp)
% returns 1D Gauss points and weights
if ngp==1
  g = 0; w = 2;
elseif ngp==2
  g = [-1/sqrt(3); 1/sqrt(3)];
  w = [1; 1];
elseif ngp==3
  g = [-sqrt(3/5); 0; sqrt(3/5)];
  w = [5/9; 8/9; 5/9];
else
  error('gauss_1d: ngp>3 not implemented')
end
end

function [N, dN_dxi] = hex8_shape(xi, eta, zeta, node_loc)
% returns shape fn N (1x8) and derivatives dN_dxi (3x8) in parent coords
N = zeros(1,8);
dN_dxi = zeros(3,8);
for a=1:8
  xi_a  = node_loc(a,1);
  eta_a = node_loc(a,2);
  zeta_a= node_loc(a,3);
  N(a) = 1/8*(1 + xi*xi_a)*(1 + eta*eta_a)*(1 + zeta*zeta_a);
  % derivatives
  dN_dxi(1,a) = 1/8* xi_a*(1+eta*eta_a)*(1+zeta*zeta_a); % dN/dxi
  dN_dxi(2,a) = 1/8* (1+xi*xi_a)*eta_a*(1+zeta*zeta_a);  % dN/deta
  dN_dxi(3,a) = 1/8* (1+xi*xi_a)*(1+eta*eta_a)*zeta_a;   % dN/dzeta
end
end

% Define element node coords (a unit cube element)
X = [0 0 0;
     1 0 0;
     1 1 0;
     0 1 0;
     0 0 1;
     1 0 1;
     1 1 1;
     0 1 1];

E = 210e9; nu = 0.3; rho = 7850;
[Ke, Me] = hex8_element(X, E, nu, rho);

% Quick checks:
disp('Ke size:'); disp(size(Ke));
disp('Me size:'); disp(size(Me));
fprintf('Ke symmetric? %g, Me symmetric? %g\n', norm(Ke-Ke','fro'), norm(Me-Me','fro'));

% Lumped mass (row-sum)
Mlumped = diag(sum(Me,2));
