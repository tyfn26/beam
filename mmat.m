function [Me] = mmat(coord, rho, A, J, I)
    % Creates the consistent mass matrix for a 3D Timoshenko beam element.
    % Local DOF order: [u1, v1, w1, theta_x1, theta_y1, theta_z1, u2, v2, w2, theta_x2, theta_y2, theta_z2]
    Lambda=zeros(3,3);
x1=coord(1,1);
x2=coord(2,1);
y1=coord(1,2);
y2=coord(2,2);
z1=coord(1,3);
z2=coord(2,3);
L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +(z2-z1)*(z2-z1));
if x1 == x2 && y1 == y2
  if z2 > z1
    Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
  else
    Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
  end
  else
      CXx = (x2-x1)/L;
CYx = (y2-y1)/L;
CZx = (z2-z1)/L;
D = sqrt(CXx*CXx + CYx*CYx);
CXy = -CYx/D;
CYy = CXx/D;
CZy = 0;
CXz = -CXx*CZx/D;
CYz = -CYx*CZx/D;
CZz = D;
Lambda = [CXx CYx CZx ;CXy CYy CZy ;CXz CYz CZz];
endif

T = zeros(12);
T([1, 2, 3], [1, 2, 3]) = Lambda;
T([4,5,6], [4,5,6]) = Lambda;
T([7,8,9], [7,8,9]) = Lambda;
T([10,11,12], [10,11,12]) = Lambda;
TT = T';
    % Initialize 12x12 matrix
    Me = zeros(12);

    % Common factor for translational mass
    mass_factor = (rho * A * L) / 6;

    % 1. Axial & Translational Inertia (u, v, w DOFs)
    % Pattern is the same for all three translational directions
    m_translational = mass_factor * [2, 1; 1, 2];

    Me([1, 7], [1, 7]) = m_translational;  % u1, u2
    Me([2, 8], [2, 8]) = m_translational;  % v1, v2
    Me([3, 9], [3, 9]) = m_translational;  % w1, w2

    % 2. Rotational Inertia (theta_x, theta_y, theta_z DOFs)
    % a) Torsional (theta_x) - analogous to axial
    mass_torsion = (rho * J * L) / 6;
    m_torsion = mass_torsion * [2, 1; 1, 2];
    Me([4, 10], [4, 10]) = m_torsion;

    % b) Bending Rotations (Timoshenko formulation)
    % For bending in x-z plane (theta_y)
    mass_rot_y = (rho * I) / L;
    m_bend_y = mass_rot_y * [2/3, 1/3; 1/3, 2/3];
    Me([5, 11], [5, 11]) = m_bend_y;

    % For bending in x-y plane (theta_z)
    mass_rot_z = (rho * I) / L;
    m_bend_z = mass_rot_z * [2/3, 1/3; 1/3, 2/3];
    Me([6, 12], [6, 12]) = m_bend_z;
    Me = T*Me*TT;
end
