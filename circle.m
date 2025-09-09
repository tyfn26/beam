clear all;clc;
% Define two connected 3D lines
% Line 1: from point A to point B
A = [0, 0, 0];
B = [1, 0, 0];

% Line 2: from point B to point C (B is the common point)
C = [2, 0, 0];

%% Calculate direction vectors
v1 = B - A;  % Direction of first line
v2 = C - B;  % Direction of second line

% Normalize direction vectors
v1_unit = v1 / norm(v1);
v2_unit = v2 / norm(v2);

%% Calculate the angle between the lines at point B
dot_product = dot(v1_unit, v2_unit);
angle_rad = acos(dot_product);
angle_deg = rad2deg(angle_rad);

fprintf('Angle at junction point: %.2f degrees\n', angle_deg);

%% Find the bisector vector
bisector = v1_unit + v2_unit;
bisector = bisector / norm(bisector);

%% Draw circle normal to the bisector
radius = 0.8;
theta = linspace(0, 2*pi, 10);
junction_point = B;

% Find two orthogonal vectors in the plane perpendicular to the bisector
% We need to find two vectors that are both perpendicular to the bisector

% Method 1: Find first perpendicular vector
if abs(bisector(1)) > abs(bisector(2))
    temp = [0, 1, 0];  % Use y-axis if x-component is large
else
    temp = [1, 0, 0];  % Use x-axis otherwise
end
u = cross(bisector, temp);
u = u / norm(u);

% Second perpendicular vector (orthogonal to both bisector and u)
w = cross(bisector, u);
w = w / norm(w);

% Generate circle points in the plane perpendicular to bisector
circle_points = zeros(3, length(theta));
for i = 1:length(theta)
    circle_points(:, i) = A' + radius * (cos(theta(i)) * u' + sin(theta(i)) * w');
end
##
##%% Plot everything
##figure;
##hold on;
##grid on;
##axis equal;
##
##% Plot the two lines
##plot3([A(1), B(1)], [A(2), B(2)], [A(3), B(3)], 'b-', 'LineWidth', 3);
##plot3([B(1), C(1)], [B(2), C(2)], [B(3), C(3)], 'r-', 'LineWidth', 3);
##
##% Plot the points
##plot3(A(1), A(2), A(3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
##plot3(B(1), B(2), B(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
##plot3(C(1), C(2), C(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
##
##% Plot the circle (normal to bisector)
plot3(circle_points(1,:), circle_points(2,:), circle_points(3,:), ...
      'm-', 'LineWidth', 3);
##
##% Plot direction vectors at junction point
##quiver3(B(1), B(2), B(3), v1_unit(1), v1_unit(2), v1_unit(3), ...
##        'b', 'LineWidth', 2, 'MaxHeadSize', 0.3);
##quiver3(B(1), B(2), B(3), v2_unit(1), v2_unit(2), v2_unit(3), ...
##        'r', 'LineWidth', 2, 'MaxHeadSize', 0.3);
##quiver3(B(1), B(2), B(3), bisector(1), bisector(2), bisector(3), ...
##        'g', 'LineWidth', 2, 'MaxHeadSize', 0.3);
##
##% Plot the plane normal vectors for visualization
##quiver3(B(1), B(2), B(3), u(1), u(2), u(3), ...
##        'c', 'LineWidth', 1.5, 'MaxHeadSize', 0.2);
##quiver3(B(1), B(2), B(3), w(1), w(2), w(3), ...
##        'y', 'LineWidth', 1.5, 'MaxHeadSize', 0.2);
##
##% Labels and title
##xlabel('X');
##ylabel('Y');
##zlabel('Z');
##title(sprintf('Circle Normal to Bisector (Angle: %.1f°)', angle_deg));
##legend('Line 1 (A→B)', 'Line 2 (B→C)', 'Point A', 'Junction Point B', 'Point C', ...
##       'Circle (⊥ to bisector)', 'v1 direction', 'v2 direction', 'Bisector', ...
##       'Circle axis u', 'Circle axis w');
##
##% Set a nice viewing angle
##view(45, 30);
##hold off;
##
##%% Display results
##fprintf('\nLine 1: A = [%.1f, %.1f, %.1f] to B = [%.1f, %.1f, %.1f]\n', A, B);
##fprintf('Line 2: B = [%.1f, %.1f, %.1f] to C = [%.1f, %.1f, %.1f]\n', B, C);
##fprintf('Direction v1: [%.3f, %.3f, %.3f]\n', v1_unit);
##fprintf('Direction v2: [%.3f, %.3f, %.3f]\n', v2_unit);
##fprintf('Bisector vector: [%.3f, %.3f, %.3f]\n', bisector);
##fprintf('Circle plane normal (u): [%.3f, %.3f, %.3f]\n', u);
##fprintf('Circle plane normal (w): [%.3f, %.3f, %.3f]\n', w);
##
##%% Verify orthogonality (should be very small numbers)
##fprintf('\nOrthogonality checks:\n');
##fprintf('u • bisector = %.6f (should be 0)\n', dot(u, bisector));
##fprintf('w • bisector = %.6f (should be 0)\n', dot(w, bisector));
##fprintf('u • w = %.6f (should be 0)\n', dot(u, w));
