function [T1,T2,R] = radius(line1_start,line1_end,line2_start,line2_end,radius,N)
##line1_start = [0, 0, 0];
##line1_end = [15, 2, 2];
##line2_start = [15, 2, 2]; % Çizgi 1 ile birleşiyor
##line2_end = [50, 2, 2];
##
##% Radyus değeri
##radius = 30;

% Çizgilerin bağlantı noktasını kontrol et
if isequal(line1_end, line2_start)
    connection_point = line1_end;
    fprintf('Çizgiler bağlantı noktasında birleşiyor.\n');
else
    error('Çizgiler bağlantı noktasında birleşmiyor.');
end

% Çizgilerin yön vektörlerini bul
dir1 = (line1_end - line1_start) / norm(line1_end - line1_start);
dir2 = (line2_end - line2_start) / norm(line2_end - line2_start);

% Düzlem normal vektörünü hesapla (yayın bulunacağı düzlem)
normal = cross(dir1, dir2);
normal = normal / norm(normal);

% Açıyı hesapla (-dir1 ve dir2 arasındaki açı)
theta = acos(dot(-dir1, dir2));

% Açıortay vektörünü bul
angle_bisector = (-dir1 + dir2);
angle_bisector = angle_bisector / norm(angle_bisector);

% Radyus merkezini hesapla
center = connection_point + angle_bisector * (radius / sin(theta/2));

% Teğet noktalarını projeksiyon ile bul
% Line 1 için teğet noktası
t1 = dot(center - line1_start, dir1);
tangent_point1 = line1_start + t1 * dir1;
T1=tangent_point1;
% Line 2 için teğet noktası
t2 = dot(center - line2_start, dir2);
tangent_point2 = line2_start + t2 * dir2;
T2=tangent_point2;
% Yay parametrelerini hesapla
v1 = tangent_point1 - center;
v2 = tangent_point2 - center;
theta_start = atan2(norm(cross(v1, v2)), dot(v1, v2));
theta_end = 0;

% 3D yayı oluştur
t = linspace(0, theta_start, N);
arc_points = zeros(length(t), 3);
for i = 1:length(t)
    arc_points(i,:) = center + radius * (v1/norm(v1)) * cos(t(i)) + radius * (cross(normal, v1)/norm(v1)) * sin(t(i));
end
R=arc_points;
R(1,:)=[];
R(size(R,1),:)=[];
##% Grafik çizim
##figure;
##hold on;
##plot3([line1_start(1), line1_end(1)], [line1_start(2), line1_end(2)], [line1_start(3), line1_end(3)], 'b', 'LineWidth', 2);
##plot3([line2_start(1), line2_end(1)], [line2_start(2), line2_end(2)], [line2_start(3), line2_end(3)], 'r', 'LineWidth', 2);
##plot3(arc_points(:,1), arc_points(:,2), arc_points(:,3), 'g', 'LineWidth', 2);
##plot3(center(1), center(2), center(3), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
##plot3(tangent_point1(1), tangent_point1(2), tangent_point1(3), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
##plot3(tangent_point2(1), tangent_point2(2), tangent_point2(3), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 8);
##
##axis equal;
##grid on;
##xlabel('X');
##ylabel('Y');
##zlabel('Z');
##legend('Çizgi 1', 'Çizgi 2', 'Radyus', 'Merkez', 'Teğet 1', 'Teğet 2');
##hold off;
