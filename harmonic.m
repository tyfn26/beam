% Harmonic Response with Base Acceleration and Stress Calculation for Brick Element
clear; close all; clc;

%% 1. Material Properties and Geometry
fprintf('Setting up brick element properties...\n');

% Material properties (Aluminum)
E = 70e9;          % Young's Modulus (Pa)
nu = 0.33;         % Poisson's Ratio
rho = 2700;        % Density (kg/m^3)
zeta = 0.02;       % Damping Ratio (2%)

% Brick element geometry
Lx = 0.1;          % Length (m)
Ly = 0.05;         % Width (m)
Lz = 0.02;         % Height (m)

% Base excitation parameters
base_accel = 1000; % mm/s² (input acceleration)
base_accel_m_s2 = base_accel / 1000; % Convert to m/s²

%% 2. Define Single Brick Element (8-node HEX8)
fprintf('Creating brick element...\n');

% Node coordinates for a single brick element
nodes = [0,    0,    0;     % Node 1
         Lx,   0,    0;     % Node 2
         Lx,   Ly,   0;     % Node 3
         0,    Ly,   0;     % Node 4
         0,    0,    Lz;    % Node 5
         Lx,   0,    Lz;    % Node 6
         Lx,   Ly,   Lz;    % Node 7
         0,    Ly,   Lz];   % Node 8

% Element connectivity
elementConnectivity = [1, 2, 3, 4, 5, 6, 7, 8];

nNodes = size(nodes, 1);
dofPerNode = 3; % x, y, z displacements
totalDof = nNodes * dofPerNode;

%% 3. Calculate Element Stiffness and Mass Matrices
fprintf('Calculating element matrices...\n');
[Ke, Me] = brick8Element(E, nu, rho, nodes, elementConnectivity);

% For single element, global matrices = element matrices
K_global = Ke;
M_global = Me;

%% 4. Apply Boundary Conditions and Base Excitation
fprintf('Applying boundary conditions...\n');

% Fixed face (nodes 1,2,3,4 - bottom face)
fixedNodes = [1, 2, 3, 4];
fixedDofs = [];
for n = fixedNodes
    fixedDofs = [fixedDofs, (n-1)*dofPerNode+1, (n-1)*dofPerNode+2, (n-1)*dofPerNode+3];
end
fixedDofs = unique(fixedDofs);

% Base excitation in Z-direction (vertical shaking)
influence_vector = zeros(totalDof, 1);
for i = 1:nNodes
    influence_vector((i-1)*dofPerNode + 3) = 1; % Z-direction DOFs
end

% Free DOFs
allDofs = 1:totalDof;
freeDofs = setdiff(allDofs, fixedDofs);

% Reduce matrices
K_red = K_global(freeDofs, freeDofs);
M_red = M_global(freeDofs, freeDofs);

% Effective force for base excitation
F_effective = -M_global * influence_vector * base_accel_m_s2;
F_red = F_effective(freeDofs);

%% 5. Calculate Natural Frequencies
fprintf('Calculating natural frequencies...\n');
[V, D] = eigs(K_red, M_red, 3, 'sm'); % Get first 3 modes
omega_n = sqrt(diag(D));
f_natural = omega_n / (2*pi);
fprintf('First natural frequency: %.2f Hz\n', f_natural(1));
fprintf('Second natural frequency: %.2f Hz\n', f_natural(2));
fprintf('Third natural frequency: %.2f Hz\n', f_natural(3));

% Use first mode for Rayleigh damping
alpha = 0; % Mass proportional damping
beta = 2 * zeta / omega_n(1); % Stiffness proportional damping
C_red = alpha * M_red + beta * K_red;

%% 6. Frequency Sweep with Stress Calculation
fprintf('Performing frequency sweep with stress calculation...\n');

freq_range = linspace(0, 1.5*f_natural(1), 300);
max_von_mises = zeros(1, length(freq_range));
max_principal = zeros(1, length(freq_range));
tip_displacement = zeros(1, length(freq_range));

% Precompute stress calculation matrices
[D, B_cell, detJ_cell, gp_weights] = precomputeStressMatrices(E, nu, nodes);

for i = 1:length(freq_range)
    omega = 2 * pi * freq_range(i);
    
    % Solve for complex displacement
    Z = (K_red - omega^2 * M_red) + 1i * omega * C_red;
    u_red_complex = Z \ F_red;
    
    % Reconstruct full displacement vector
    u_full = zeros(totalDof, 1);
    u_full(freeDofs) = u_red_complex;
    
    % Store tip displacement (node 8, Z-direction)
    tip_dof_z = (8-1)*dofPerNode + 3;
    tip_displacement(i) = abs(u_full(tip_dof_z)) * 1000; % mm
    
    % Calculate stresses at Gauss points
    [von_mises_stress, principal_stress] = calculateBrickStresses(u_full, D, B_cell, detJ_cell, gp_weights);
    
    max_von_mises(i) = max(von_mises_stress);
    max_principal(i) = max(principal_stress);
end

%% 7. Plot Results
fprintf('Plotting results...\n');

figure;

% Von Mises stress response
subplot(3, 1, 1);
plot(freq_range, max_von_mises / 1e6, 'r-', 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Max Von Mises Stress (MPa)');
title(sprintf('Maximum Stress (Base Acceleration: %d mm/s²)', base_accel));

% Mark resonant frequency
[peak_stress, peak_idx] = max(max_von_mises);
resonant_freq = freq_range(peak_idx);
hold on;
plot(resonant_freq, peak_stress/1e6, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
legend('Von Mises Stress', sprintf('Resonance: %.1f Hz', resonant_freq));

% Principal stress response
subplot(3, 1, 2);
plot(freq_range, max_principal / 1e6, 'b-', 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Max Principal Stress (MPa)');
title('Maximum Principal Stress');

% Displacement response
subplot(3, 1, 3);
plot(freq_range, tip_displacement, 'g-', 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Tip Displacement (mm)');
title('Tip Displacement (Node 8 Z-direction)');

%% 8. Display Key Results
fprintf('\n=== RESULTS SUMMARY ===\n');
fprintf('First natural frequency: %.2f Hz\n', f_natural(1));
fprintf('Resonant frequency: %.2f Hz\n', resonant_freq);
fprintf('Peak Von Mises stress: %.2f MPa\n', peak_stress/1e6);
fprintf('Peak principal stress: %.2f MPa\n', max(max_principal)/1e6);
fprintf('Peak tip displacement: %.2f mm\n', max(tip_displacement));

%% 9. Stress Distribution at Resonance
fprintf('\nCalculating stress distribution at resonance...\n');

% Get response at resonant frequency
omega_res = 2 * pi * resonant_freq;
Z = (K_red - omega_res^2 * M_red) + 1i * omega_res * C_red;
u_red_res = Z \ F_red;
u_full_res = zeros(totalDof, 1);
u_full_res(freeDofs) = u_red_res;

% Calculate stress distribution
[von_mises_res, principal_res, stress_locations] = calculateBrickStressDistribution(u_full_res, D, nodes);

figure;
scatter3(stress_locations(:,1), stress_locations(:,2), stress_locations(:,3), ...
         50, von_mises_res/1e6, 'filled');
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Von Mises Stress Distribution at Resonance (MPa)');
colormap('jet');

%% Helper Functions

function [Ke, Me] = brick8Element(E, nu, rho, nodes, conn)
    % Gauss points and weights (2x2x2)
    gp = [-1/sqrt(3), 1/sqrt(3)];
    w = [1, 1];
    
    nNodes = 8;
    dofPerNode = 3;
    elemDof = nNodes * dofPerNode;
    
    Ke = zeros(elemDof, elemDof);
    Me = zeros(elemDof, elemDof);
    
    % Constitutive matrix for 3D linear elasticity
    D = E/((1+nu)*(1-2*nu)) * [1-nu,   nu,   nu,        0,        0,        0;
                                 nu, 1-nu,   nu,        0,        0,        0;
                                 nu,   nu, 1-nu,        0,        0,        0;
                                  0,    0,    0, (1-2*nu)/2,      0,        0;
                                  0,    0,    0,        0, (1-2*nu)/2,      0;
                                  0,    0,    0,        0,        0, (1-2*nu)/2];
    
    % Loop over Gauss points
    for i = 1:length(gp)
        xi = gp(i);
        w_i = w(i);
        for j = 1:length(gp)
            eta = gp(j);
            w_j = w(j);
            for k = 1:length(gp)
                zeta_val = gp(k);
                w_k = w(k);
                
                % Get shape functions and derivatives
                [N, dN_dxi, dN_deta, dN_dzeta] = shapeFunc8Node(xi, eta, zeta_val);
                
                % Compute Jacobian matrix
                J = zeros(3, 3);
                for node = 1:nNodes
                    J(1,1) = J(1,1) + dN_dxi(node) * nodes(node, 1);
                    J(1,2) = J(1,2) + dN_dxi(node) * nodes(node, 2);
                    J(1,3) = J(1,3) + dN_dxi(node) * nodes(node, 3);
                    
                    J(2,1) = J(2,1) + dN_deta(node) * nodes(node, 1);
                    J(2,2) = J(2,2) + dN_deta(node) * nodes(node, 2);
                    J(2,3) = J(2,3) + dN_deta(node) * nodes(node, 3);
                    
                    J(3,1) = J(3,1) + dN_dzeta(node) * nodes(node, 1);
                    J(3,2) = J(3,2) + dN_dzeta(node) * nodes(node, 2);
                    J(3,3) = J(3,3) + dN_dzeta(node) * nodes(node, 3);
                end
                
                detJ = det(J);
                invJ = inv(J);
                
                % Compute derivatives w.r.t. global coordinates
                dN_dxyz = zeros(3, nNodes);
                for node = 1:nNodes
                    dN_dlocal = [dN_dxi(node); dN_deta(node); dN_dzeta(node)];
                    dN_dglobal = invJ * dN_dlocal;
                    dN_dxyz(1, node) = dN_dglobal(1);
                    dN_dxyz(2, node) = dN_dglobal(2);
                    dN_dxyz(3, node) = dN_dglobal(3);
                end
                
                % Form B matrix
                B = zeros(6, elemDof);
                for node = 1:nNodes
                    col_start = (node-1)*dofPerNode + 1;
                    B(1, col_start)   = dN_dxyz(1, node);
                    B(2, col_start+1) = dN_dxyz(2, node);
                    B(3, col_start+2) = dN_dxyz(3, node);
                    B(4, col_start)   = dN_dxyz(2, node);
                    B(4, col_start+1) = dN_dxyz(1, node);
                    B(5, col_start+1) = dN_dxyz(3, node);
                    B(5, col_start+2) = dN_dxyz(2, node);
                    B(6, col_start)   = dN_dxyz(3, node);
                    B(6, col_start+2) = dN_dxyz(1, node);
                end
                
                % Form N matrix for mass matrix
                N_mat = zeros(3, elemDof);
                for node = 1:nNodes
                    col_start = (node-1)*dofPerNode + 1;
                    N_mat(1, col_start)   = N(node);
                    N_mat(2, col_start+1) = N(node);
                    N_mat(3, col_start+2) = N(node);
                end
                
                % Add to matrices
                Ke = Ke + (B' * D * B) * w_i * w_j * w_k * detJ;
                Me = Me + (N_mat' * rho * N_mat) * w_i * w_j * w_k * detJ;
            end
        end
    end
end

function [D, B_cell, detJ_cell, gp_weights] = precomputeStressMatrices(E, nu, nodes)
    % Precompute matrices needed for stress calculation
    gp = [-1/sqrt(3), 1/sqrt(3)];
    w = [1, 1];
    nGauss = length(gp)^3;
    
    D = E/((1+nu)*(1-2*nu)) * [1-nu,   nu,   nu,        0,        0,        0;
                                 nu, 1-nu,   nu,        0,        0,        0;
                                 nu,   nu, 1-nu,        0,        0,        0;
                                  0,    0,    0, (1-2*nu)/2,      0,        0;
                                  0,    0,    0,        0, (1-2*nu)/2,      0;
                                  0,    0,    0,        0,        0, (1-2*nu)/2];
    
    B_cell = cell(1, nGauss);
    detJ_cell = zeros(1, nGauss);
    gp_weights = zeros(1, nGauss);
    
    idx = 1;
    for i = 1:length(gp)
        for j = 1:length(gp)
            for k = 1:length(gp)
                xi = gp(i); eta = gp(j); zeta_val = gp(k);
                [~, dN_dxi, dN_deta, dN_dzeta] = shapeFunc8Node(xi, eta, zeta_val);
                
                % Compute Jacobian
                J = zeros(3, 3);
                for node = 1:8
                    J(1,1) = J(1,1) + dN_dxi(node) * nodes(node, 1);
                    J(1,2) = J(1,2) + dN_dxi(node) * nodes(node, 2);
                    J(1,3) = J(1,3) + dN_dxi(node) * nodes(node, 3);
                    
                    J(2,1) = J(2,1) + dN_deta(node) * nodes(node, 1);
                    J(2,2) = J(2,2) + dN_deta(node) * nodes(node, 2);
                    J(2,3) = J(2,3) + dN_deta(node) * nodes(node, 3);
                    
                    J(3,1) = J(3,1) + dN_dzeta(node) * nodes(node, 1);
                    J(3,2) = J(3,2) + dN_dzeta(node) * nodes(node, 2);
                    J(3,3) = J(3,3) + dN_dzeta(node) * nodes(node, 3);
                end
                
                detJ = det(J);
                invJ = inv(J);
                
                % Compute global derivatives and form B matrix
                dN_dxyz = zeros(3, 8);
                B = zeros(6, 24);
                for node = 1:8
                    dN_dlocal = [dN_dxi(node); dN_deta(node); dN_dzeta(node)];
                    dN_dglobal = invJ * dN_dlocal;
                    dN_dxyz(1, node) = dN_dglobal(1);
                    dN_dxyz(2, node) = dN_dglobal(2);
                    dN_dxyz(3, node) = dN_dglobal(3);
                    
                    col_start = (node-1)*3 + 1;
                    B(1, col_start)   = dN_dglobal(1);
                    B(2, col_start+1) = dN_dglobal(2);
                    B(3, col_start+2) = dN_dglobal(3);
                    B(4, col_start)   = dN_dglobal(2);
                    B(4, col_start+1) = dN_dglobal(1);
                    B(5, col_start+1) = dN_dglobal(3);
                    B(5, col_start+2) = dN_dglobal(2);
                    B(6, col_start)   = dN_dglobal(3);
                    B(6, col_start+2) = dN_dglobal(1);
                end
                
                B_cell{idx} = B;
                detJ_cell(idx) = detJ;
                gp_weights(idx) = w(i) * w(j) * w(k);
                idx = idx + 1;
            end
        end
    end
end

function [von_mises, principal] = calculateBrickStresses(u, D, B_cell, detJ_cell, gp_weights)
    % Calculate stresses at all Gauss points
    nGauss = length(B_cell);
    von_mises = zeros(1, nGauss);
    principal = zeros(1, nGauss);
    
    for i = 1:nGauss
        B = B_cell{i};
        strain = B * u;
        stress = D * strain;
        
        % Von Mises stress
        sxx = stress(1); syy = stress(2); szz = stress(3);
        sxy = stress(4); syz = stress(5); szx = stress(6);
        von_mises(i) = sqrt(0.5*((sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 + ...
                          6*(sxy^2 + syz^2 + szx^2)));
        
        % Principal stress (max)
        stress_matrix = [sxx, sxy, szx;
                        sxy, syy, syz;
                        szx, syz, szz];
        principal_stresses = eig(stress_matrix);
        principal(i) = max(principal_stresses);
    end
end

function [von_mises, principal, locations] = calculateBrickStressDistribution(u, D, nodes)
    % Calculate stress distribution throughout element
    gp = [-1/sqrt(3), 1/sqrt(3)];
    nGauss = length(gp)^3;
    
    von_mises = zeros(1, nGauss);
    principal = zeros(1, nGauss);
    locations = zeros(nGauss, 3);
    
    idx = 1;
    for i = 1:length(gp)
        for j = 1:length(gp)
            for k = 1:length(gp)
                xi = gp(i); eta = gp(j); zeta_val = gp(k);
                
                % Get shape functions for this Gauss point
                [N, dN_dxi, dN_deta, dN_dzeta] = shapeFunc8Node(xi, eta, zeta_val);
                
                % Compute location of this Gauss point
                x = sum(N .* nodes(:,1));
                y = sum(N .* nodes(:,2));
                z = sum(N .* nodes(:,3));
                locations(idx, :) = [x, y, z];
                
                % Compute Jacobian and B matrix
                J = zeros(3, 3);
                for node = 1:8
                    J(1,1) = J(1,1) + dN_dxi(node) * nodes(node, 1);
                    J(1,2) = J(1,2) + dN_dxi(node) * nodes(node, 2);
                    J(1,3) = J(1,3) + dN_dxi(node) * nodes(node, 3);
                    
                    J(2,1) = J(2,1) + dN_deta(node) * nodes(node, 1);
                    J(2,2) = J(2,2) + dN_deta(node) * nodes(node, 2);
                    J(2,3) = J(2,3) + dN_deta(node) * nodes(node, 3);
                    
                    J(3,1) = J(3,1) + dN_dzeta(node) * nodes(node, 1);
                    J(3,2) = J(3,2) + dN_dzeta(node) * nodes(node, 2);
                    J(3,3) = J(3,3) + dN_dzeta(node) * nodes(node, 3);
                end
                
                invJ = inv(J);
                
                % Form B matrix
                B = zeros(6, 24);
                for node = 1:8
                    dN_dlocal = [dN_dxi(node); dN_deta(node); dN_dzeta(node)];
                    dN_dglobal = invJ * dN_dlocal;
                    
                    col_start = (node-1)*3 + 1;
                    B(1, col_start)   = dN_dglobal(1);
                    B(2, col_start+1) = dN_dglobal(2);
                    B(3, col_start+2) = dN_dglobal(3);
                    B(4, col_start)   = dN_dglobal(2);
                    B(4, col_start+1) = dN_dglobal(1);
                    B(5, col_start+1) = dN_dglobal(3);
                    B(5, col_start+2) = dN_dglobal(2);
                    B(6, col_start)   = dN_dglobal(3);
                    B(6, col_start+2) = dN_dglobal(1);
                end
                
                % Calculate stress
                strain = B * u;
                stress = D * strain;
                
                % Von Mises stress
                sxx = stress(1); syy = stress(2); szz = stress(3);
                sxy = stress(4); syz = stress(5); szx = stress(6);
                von_mises(idx) = sqrt(0.5*((sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 + ...
                                  6*(sxy^2 + syz^2 + szx^2)));
                
                % Principal stress
                stress_matrix = [sxx, sxy, szx;
                                sxy, syy, syz;
                                szx, syz, szz];
                principal_stresses = eig(stress_matrix);
                principal(idx) = max(principal_stresses);
                
                idx = idx + 1;
            end
        end
    end
end

function [N, dN_dxi, dN_deta, dN_dzeta] = shapeFunc8Node(xi, eta, zeta)
    % Shape functions for 8-node brick element
    N = zeros(8,1);
    dN_dxi   = zeros(8,1);
    dN_deta  = zeros(8,1);
    dN_dzeta = zeros(8,1);
    
    N(1) = (1-xi)*(1-eta)*(1-zeta)/8;
    N(2) = (1+xi)*(1-eta)*(1-zeta)/8;
    N(3) = (1+xi)*(1+eta)*(1-zeta)/8;
    N(4) = (1-xi)*(1+eta)*(1-zeta)/8;
    N(5) = (1-xi)*(1-eta)*(1+zeta)/8;
    N(6) = (1+xi)*(1-eta)*(1+zeta)/8;
    N(7) = (1+xi)*(1+eta)*(1+zeta)/8;
    N(8) = (1-xi)*(1+eta)*(1+zeta)/8;
    
    dN_dxi(1) = -(1-eta)*(1-zeta)/8;
    dN_dxi(2) =  (1-eta)*(1-zeta)/8;
    dN_dxi(3) =  (1+eta)*(1-zeta)/8;
    dN_dxi(4) = -(1+eta)*(1-zeta)/8;
    dN_dxi(5) = -(1-eta)*(1+zeta)/8;
    dN_dxi(6) =  (1-eta)*(1+zeta)/8;
    dN_dxi(7) =  (1+eta)*(1+zeta)/8;
    dN_dxi(8) = -(1+eta)*(1+zeta)/8;
    
    dN_deta(1) = -(1-xi)*(1-zeta)/8;
    dN_deta(2) = -(1+xi)*(1-zeta)/8;
    dN_deta(3) =  (1+xi)*(1-zeta)/8;
    dN_deta(4) =  (1-xi)*(1-zeta)/8;
    dN_deta(5) = -(1-xi)*(1+zeta)/8;
    dN_deta(6) = -(1+xi)*(1+zeta)/8;
    dN_deta(7) =  (1+xi)*(1+zeta)/8;
    dN_deta(8) =  (1-xi)*(1+zeta)/8;
    
    dN_dzeta(1) = -(1-xi)*(1-eta)/8;
    dN_dzeta(2) = -(1+xi)*(1-eta)/8;
    dN_dzeta(3) = -(1+xi)*(1+eta)/8;
    dN_dzeta(4) = -(1-xi)*(1+eta)/8;
    dN_dzeta(5) =  (1-xi)*(1-eta)/8;
    dN_dzeta(6) =  (1+xi)*(1-eta)/8;
    dN_dzeta(7) =  (1+xi)*(1+eta)/8;
    dN_dzeta(8) =  (1-xi)*(1+eta)/8;
end
