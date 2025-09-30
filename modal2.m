%% Çelik Prizma Free-Free Modal Analizi - DÜZELTİLMİŞ
% 2x1x3 mm boyutları, 1 mm hex elemanlar
clear all; close all; clc;

%% Malzeme Özellikleri (Çelik)
E = 200e9;       % Pa, Elastisite modülü
nu = 0.3;        % Poisson oranı
rho = 7850;      % kg/m^3, Yoğunluk

%% Geometri ve Mesh
Lx = 2e-3; Ly = 1e-3; Lz = 3e-3;  % m, Boyutlar
elem_size = 1e-3;                  % m, Eleman boyutu

% Mesh parametreleri
nx = Lx/elem_size; ny = Ly/elem_size; nz = Lz/elem_size;
num_nodes_x = nx + 1;
num_nodes_y = ny + 1; 
num_nodes_z = nz + 1;
num_nodes = num_nodes_x * num_nodes_y * num_nodes_z;
num_elements = nx * ny * nz;

fprintf('Mesh Bilgisi:\n');
fprintf('  Düğüm sayısı: %d\n', num_nodes);
fprintf('  Eleman sayısı: %d\n', num_elements);
fprintf('  Serbestlik derecesi: %d\n', num_nodes*3);

%% Düğüm Koordinatlarını Oluştur
nodes = zeros(num_nodes, 3);
node_id = 1;

for k = 1:num_nodes_z
    for j = 1:num_nodes_y
        for i = 1:num_nodes_x
            nodes(node_id, :) = [(i-1)*elem_size, (j-1)*elem_size, (k-1)*elem_size];
            node_id = node_id + 1;
        end
    end
end

%% Eleman Bağlantılarını Oluştur
elements = zeros(num_elements, 8);
elem_id = 1;

for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            n0 = (i-1) + (j-1)*num_nodes_x + (k-1)*num_nodes_x*num_nodes_y;
            n1 = n0 + 1;
            n2 = n0 + num_nodes_x + 1;
            n3 = n0 + num_nodes_x;
            n4 = n0 + num_nodes_x*num_nodes_y;
            n5 = n1 + num_nodes_x*num_nodes_y;
            n6 = n2 + num_nodes_x*num_nodes_y;
            n7 = n3 + num_nodes_x*num_nodes_y;
            
            elements(elem_id, :) = [n0+1, n1+1, n2+1, n3+1, n4+1, n5+1, n6+1, n7+1];
            elem_id = elem_id + 1;
        end
    end
end

%% DÜZELTİLMİŞ Hex8 Eleman için Sertlik Matrisi Fonksiyonu
function Ke = hex8_stiffness_corrected(E, nu, nodes)
    % Düzeltilmiş D matrisi
    G = E / (2 * (1 + nu));
    lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    
    D = [lambda + 2*G, lambda, lambda, 0, 0, 0;
         lambda, lambda + 2*G, lambda, 0, 0, 0;
         lambda, lambda, lambda + 2*G, 0, 0, 0;
         0, 0, 0, G, 0, 0;
         0, 0, 0, 0, G, 0;
         0, 0, 0, 0, 0, G];
    
    % 8 noktalı tam integrasyon
    gp = [-0.577350269189626, 0.577350269189626];
    w = [1, 1];
    
    Ke = zeros(24, 24);
    
    for i = 1:2
        for j = 1:2
            for k = 1:2
                xi = gp(i); eta = gp(j); zeta = gp(k);
                
                % DÜZELTİLMİŞ şekil fonksiyonları türevleri
                dN_dxi = 0.125 * [
                    -(1-eta)*(1-zeta),  (1-eta)*(1-zeta),  (1+eta)*(1-zeta), -(1+eta)*(1-zeta), ...
                    -(1-eta)*(1+zeta),  (1-eta)*(1+zeta),  (1+eta)*(1+zeta), -(1+eta)*(1+zeta);
                    
                    -(1-xi)*(1-zeta), -(1+xi)*(1-zeta),  (1+xi)*(1-zeta),  (1-xi)*(1-zeta), ...
                    -(1-xi)*(1+zeta), -(1+xi)*(1+zeta),  (1+xi)*(1+zeta),  (1-xi)*(1+zeta);
                    
                    -(1-xi)*(1-eta), -(1+xi)*(1-eta), -(1+xi)*(1+eta), -(1-xi)*(1+eta), ...
                     (1-xi)*(1-eta),  (1+xi)*(1-eta),  (1+xi)*(1+eta),  (1-xi)*(1+eta)
                ];
                
                % Jacobian matrisi
                J = dN_dxi * nodes;
                detJ = det(J);
                
                if detJ <= 0
                    error('Negatif Jacobian determinantı!');
                end
                
                J_inv = inv(J);
                
                % Şekil fonksiyonlarının global türevleri
                dN_dx = J_inv * dN_dxi;
                
                % B matrisi (strain-displacement)
                B = zeros(6, 24);
                for m = 1:8
                    % Strain-displacement ilişkisi
                    B(1, 3*m-2) = dN_dx(1, m);  % ε_xx
                    B(2, 3*m-1) = dN_dx(2, m);  % ε_yy  
                    B(3, 3*m)   = dN_dx(3, m);  % ε_zz
                    
                    B(4, 3*m-2) = dN_dx(2, m);  % γ_xy
                    B(4, 3*m-1) = dN_dx(1, m);
                    
                    B(5, 3*m-1) = dN_dx(3, m);  % γ_yz
                    B(5, 3*m)   = dN_dx(2, m);
                    
                    B(6, 3*m-2) = dN_dx(3, m);  % γ_zx
                    B(6, 3*m)   = dN_dx(1, m);
                end
                
                % Sertlik matrisine katkı
                Ke = Ke + (B' * D * B) * detJ * w(i) * w(j) * w(k);
            end
        end
    end
    
    % Simetri kontrolü ve küçük asimetrileri düzelt
    Ke = 0.5 * (Ke + Ke');
end

%% Basitleştirilmiş Kütle Matrisi (Konsantre Kütle)
function Me = hex8_lumped_mass(rho, nodes)
    % Hacim hesapla (8 noktalı integrasyonla)
    gp = [-0.577350269189626, 0.577350269189626];
    w = [1, 1];
    
    volume = 0;
    for i = 1:2
        for j = 1:2
            for k = 1:2
                xi = gp(i); eta = gp(j); zeta = gp(k);
                
                dN_dxi = 0.125 * [
                    -(1-eta)*(1-zeta),  (1-eta)*(1-zeta),  (1+eta)*(1-zeta), -(1+eta)*(1-zeta), ...
                    -(1-eta)*(1+zeta),  (1-eta)*(1+zeta),  (1+eta)*(1+zeta), -(1+eta)*(1+zeta);
                    
                    -(1-xi)*(1-zeta), -(1+xi)*(1-zeta),  (1+xi)*(1-zeta),  (1-xi)*(1-zeta), ...
                    -(1-xi)*(1+zeta), -(1+xi)*(1+zeta),  (1+xi)*(1+zeta),  (1-xi)*(1+zeta);
                    
                    -(1-xi)*(1-eta), -(1+xi)*(1-eta), -(1+xi)*(1+eta), -(1-xi)*(1+eta), ...
                     (1-xi)*(1-eta),  (1+xi)*(1-eta),  (1+xi)*(1+eta),  (1-xi)*(1+eta)
                ];
                
                J = dN_dxi * nodes;
                detJ = det(J);
                volume = volume + detJ * w(i) * w(j) * w(k);
            end
        end
    end
    
    % Toplam kütle
    total_mass = rho * volume;
    
    % Konsantre kütle matrisi (her düğüme eşit kütle dağıt)
    node_mass = total_mass / 8;
    Me = zeros(24, 24);
    for i = 1:8
        Me(3*i-2, 3*i-2) = node_mass;
        Me(3*i-1, 3*i-1) = node_mass;
        Me(3*i, 3*i) = node_mass;
    end
end

%% Global Matrisleri Oluştur
fprintf('Global matrisler oluşturuluyor...\n');

K_global = zeros(num_nodes*3, num_nodes*3);
M_global = zeros(num_nodes*3, num_nodes*3);

% Her eleman için matrisleri hesapla ve assemble et
for elem_idx = 1:num_elements
    if mod(elem_idx, 10) == 0
        fprintf('  Eleman %d/%d işleniyor...\n', elem_idx, num_elements);
    end
    
    % Eleman düğüm koordinatları
    elem_nodes = nodes(elements(elem_idx, :), :);
    
    % Eleman matrisleri
    Ke = hex8_stiffness_corrected(E, nu, elem_nodes);
    Me = hex8_lumped_mass(rho, elem_nodes);
    
    % Global matrise assemble et
    for i = 1:8
        for j = 1:8
            node_i = elements(elem_idx, i);
            node_j = elements(elem_idx, j);
            
            dof_i = (node_i-1)*3 + (1:3);
            dof_j = (node_j-1)*3 + (1:3);
            
            K_global(dof_i, dof_j) = K_global(dof_i, dof_j) + Ke(3*i-2:3*i, 3*j-2:3*j);
            M_global(dof_i, dof_j) = M_global(dof_i, dof_j) + Me(3*i-2:3*i, 3*j-2:3*j);
        end
    end
end

fprintf('Global matrisler oluşturuldu.\n');

%% Modal Analiz - Geliştirilmiş
fprintf('Modal analiz yapılıyor...\n');

% Küçük bir regularization ekle (rigid body modları için)
K_reg = K_global + 1e-9 * max(diag(K_global)) * eye(size(K_global));

% Özdeğer problemi çöz (sparse solver kullan)
[V, D] = eigs(K_reg, M_global, 20, 'sm');

% Özdeğerleri sırala
lambda = diag(D);
[lambda, idx] = sort(lambda);
frequencies = sqrt(abs(lambda)) / (2*pi);
V = V(:, idx);

%% ANSYS SONUÇLARI İLE KARŞILAŞTIRMA
ansys_freqs = [3.8e5, 4.4e5, 6.1e5, 7.4e5, 8.0e5, 8.16e5];
num_ansys_modes = length(ansys_freqs);

fprintf('\nANSYS ile Karşılaştırma:\n');
fprintf('Mod |  MATLAB (Hz)  |  ANSYS (Hz)   |  Fark (%%)\n');
fprintf('----|---------------|---------------|----------\n');

for i = 1:min(12, length(frequencies))
    if i <= num_ansys_modes
        freq_matlab = frequencies(6+i); % İlk 6 rigid body modu atla
        freq_ansys = ansys_freqs(i);
        error_percent = abs(freq_matlab - freq_ansys) / freq_ansys * 100;
        fprintf('%3d | %12.0f | %12.0f | %8.1f\n', i, freq_matlab, freq_ansys, error_percent);
    else
        fprintf('%3d | %12.0f |              |\n', i, frequencies(6+i));
    end
end

%% Hassasiyet Analizi
fprintf('\nHassasiyet Analizi:\n');

% Malzeme parametrelerini değiştirerek test et
E_test = E * [0.95, 1.0, 1.05]; % ±%5 değişim

for e_idx = 1:length(E_test)
    K_test = K_global * (E_test(e_idx)/E);
    [~, D_test] = eigs(K_test, M_global, 1, 'sm');
    freq_test = sqrt(abs(diag(D_test))) / (2*pi);
    fprintf('E = %.2fE: f1 = %.0f Hz\n', E_test(e_idx)/E, freq_test(1));
end

%% Görselleştirme
figure('Position', [100, 100, 1400, 900]);

% 1. Mesh yapısı
subplot(2,4,1);
plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
hold on;
for i = 1:size(elements, 1)
    elem_nodes = nodes(elements(i, :), :);
    plot3(elem_nodes([1,2,3,4,1],1), elem_nodes([1,2,3,4,1],2), elem_nodes([1,2,3,4,1],3), 'b-');
    plot3(elem_nodes([5,6,7,8,5],1), elem_nodes([5,6,7,8,5],2), elem_nodes([5,6,7,8,5],3), 'b-');
    for j = 1:4
        plot3([elem_nodes(j,1), elem_nodes(j+4,1)], ...
              [elem_nodes(j,2), elem_nodes(j+4,2)], ...
              [elem_nodes(j,3), elem_nodes(j+4,3)], 'b-');
    end
end
title('Mesh Yapısı');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
axis equal; grid on;

% 2-4. İlk 3 esnek mod
for i = 1:3
    subplot(2,4,i+1);
    mode_num = 6 + i;
    mode_shape = V(:, mode_num);
    
    % Deforme edilmiş şekil
    scale = 0.2 / max(abs(mode_shape));
    deformed_nodes = nodes + scale * reshape(mode_shape, 3, [])';
    
    % Orijinal ve deforme edilmiş yapıyı çiz
    plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'bo', 'MarkerSize', 2);
    hold on;
    plot3(deformed_nodes(:,1), deformed_nodes(:,2), deformed_nodes(:,3), 'ro', 'MarkerSize', 3);
    
    for elem_idx = 1:size(elements, 1)
        orig_elem = nodes(elements(elem_idx, :), :);
        def_elem = deformed_nodes(elements(elem_idx, :), :);
        
        % Orijinal eleman
        for edge = [1,2,3,4,1; 5,6,7,8,5; 1,5; 2,6; 3,7; 4,8]'
            plot3(orig_elem(edge,1), orig_elem(edge,2), orig_elem(edge,3), 'b-', 'LineWidth', 0.5);
        end
        
        % Deforme edilmiş eleman
        for edge = [1,2,3,4,1; 5,6,7,8,5; 1,5; 2,6; 3,7; 4,8]'
            plot3(def_elem(edge,1), def_elem(edge,2), def_elem(edge,3), 'r-', 'LineWidth', 2);
        end
    end
    
    title(sprintf('Mod %d: %.0f Hz', i, frequencies(mode_num)));
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    axis equal; grid on;
end

% 5. Frekans karşılaştırması
subplot(2,4,5);
mod_numbers = 1:min(6, num_ansys_modes);
matlab_freqs = frequencies(7:7+length(mod_numbers)-1);
bar_data = [ansys_freqs(mod_numbers)'/1000, matlab_freqs/1000];
bar(mod_numbers, bar_data);
legend('ANSYS', 'MATLAB', 'Location', 'northwest');
xlabel('Mod Numarası');
ylabel('Frekans (kHz)');
title('Frekans Karşılaştırması');
grid on;

% 6. Hata analizi
subplot(2,4,6);
error_percent = abs(matlab_freqs - ansys_freqs(mod_numbers)') ./ ansys_freqs(mod_numbers)' * 100;
bar(mod_numbers, error_percent);
xlabel('Mod Numarası');
ylabel('Hata (%)');
title('ANSYS-MATLAB Farkı');
grid on;

% 7. Matris yoğunluğu
subplot(2,4,7);
spy(K_global);
title('Sertlik Matrisi Yoğunluğu');

% 8. Özdeğer dağılımı
subplot(2,4,8);
semilogy(1:20, abs(lambda(1:20)), 'bo-', 'LineWidth', 2);
xlabel('Mod Numarası');
ylabel('|Özdeğer|');
title('Özdeğer Dağılımı');
grid on;

%% Sonuçları Kaydet
fprintf('\nSonuçlar kaydediliyor...\n');

results = table((1:20)', frequencies(1:20), 'VariableNames', {'Mod_No', 'Frekans_Hz'});
writetable(results, 'modal_analiz_sonuclari_duzeltilmis.csv');

comparison_results = table((1:num_ansys_modes)', ansys_freqs', frequencies(7:7+num_ansys_modes-1), ...
    'VariableNames', {'Mod_No', 'ANSYS_Hz', 'MATLAB_Hz'});
writetable(comparison_results, 'ansys_matlab_karsilastirma.csv');

fprintf('Analiz tamamlandı.\n');

% Toplam kütle ve enerji kontrolü
total_mass = trace(M_global) / 3;
theoretical_mass = rho * Lx * Ly * Lz;
fprintf('\nKütle Kontrolü:\n');
fprintf('Hesaplanan: %.6e kg\n', total_mass);
fprintf('Teorik:     %.6e kg\n', theoretical_mass);
fprintf('Fark: %.2f%%\n', abs(total_mass - theoretical_mass)/theoretical_mass * 100);
