%% ÇELİK PRİZMA FREE-FREE MODAL ANALİZİ
% 2x1x3 mm boyutları, 1 mm hex elemanlar
% Tutarlı kütle matrisi ile
clear all; close all; clc;

%% MALZEME ÖZELLİKLERİ (Çelik)
E = 200e9;       % Pa, Elastisite modülü
nu = 0.3;        % Poisson oranı
rho = 7850;      % kg/m^3, Yoğunluk

%% GEOMETRİ VE MESH
Lx = 2e-3; Ly = 1e-3; Lz = 3e-3;  % m, Boyutlar
elem_size = 1e-3;                  % m, Eleman boyutu

% Mesh parametreleri
nx = Lx/elem_size; ny = Ly/elem_size; nz = Lz/elem_size;
num_nodes_x = nx + 1;
num_nodes_y = ny + 1; 
num_nodes_z = nz + 1;
num_nodes = num_nodes_x * num_nodes_y * num_nodes_z;
num_elements = nx * ny * nz;

fprintf('=== MESH BİLGİSİ ===\n');
fprintf('Boyutlar: %.1fx%.1fx%.1f mm\n', Lx*1000, Ly*1000, Lz*1000);
fprintf('Düğüm sayısı: %d\n', num_nodes);
fprintf('Eleman sayısı: %d\n', num_elements);
fprintf('Serbestlik derecesi: %d\n', num_nodes*3);
fprintf('\n');

%% DÜĞÜM KOORDİNATLARINI OLUŞTUR
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

%% ELEMAN BAĞLANTILARINI OLUŞTUR
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

%% HEX8 ELEMAN İÇİN SERTLİK MATRİSİ FONKSİYONU
function Ke = hex8_stiffness(E, nu, nodes)
    % D matrisi (3D izotropik malzeme)
    G = E / (2 * (1 + nu));
    lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    
    D = [lambda + 2*G, lambda, lambda, 0, 0, 0;
         lambda, lambda + 2*G, lambda, 0, 0, 0;
         lambda, lambda, lambda + 2*G, 0, 0, 0;
         0, 0, 0, G, 0, 0;
         0, 0, 0, 0, G, 0;
         0, 0, 0, 0, 0, G];
    
    % Gauss noktaları (2x2x2)
    gp = [-0.577350269189626, 0.577350269189626];
    w = [1, 1];
    
    Ke = zeros(24, 24);
    
    for i = 1:2
        for j = 1:2
            for k = 1:2
                xi = gp(i); eta = gp(j); zeta = gp(k);
                
                % Şekil fonksiyonları türevleri
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
                J_inv = inv(J);
                
                % Şekil fonksiyonlarının global türevleri
                dN_dx = J_inv * dN_dxi;
                
                % B matrisi
                B = zeros(6, 24);
                for m = 1:8
                    B(1, 3*m-2) = dN_dx(1, m);
                    B(2, 3*m-1) = dN_dx(2, m);
                    B(3, 3*m)   = dN_dx(3, m);
                    
                    B(4, 3*m-2) = dN_dx(2, m);
                    B(4, 3*m-1) = dN_dx(1, m);
                    
                    B(5, 3*m-1) = dN_dx(3, m);
                    B(5, 3*m)   = dN_dx(2, m);
                    
                    B(6, 3*m-2) = dN_dx(3, m);
                    B(6, 3*m)   = dN_dx(1, m);
                end
                
                % Sertlik matrisine katkı
                Ke = Ke + B' * D * B * detJ * w(i) * w(j) * w(k);
            end
        end
    end
    
    % Simetri kontrolü
    Ke = 0.5 * (Ke + Ke');
end

%% HEX8 ELEMAN İÇİN TUTARLI KÜTLE MATRİSİ FONKSİYONU
function Me = hex8_consistent_mass(rho, nodes)
    % Gauss noktaları
    gp = [-0.577350269189626, 0.577350269189626];
    w = [1, 1];
    
    Me = zeros(24, 24);
    
    for i = 1:2
        for j = 1:2
            for k = 1:2
                xi = gp(i); eta = gp(j); zeta = gp(k);
                
                % Şekil fonksiyonları
                N = zeros(8, 1);
                N(1) = 0.125 * (1-xi)*(1-eta)*(1-zeta);
                N(2) = 0.125 * (1+xi)*(1-eta)*(1-zeta);
                N(3) = 0.125 * (1+xi)*(1+eta)*(1-zeta);
                N(4) = 0.125 * (1-xi)*(1+eta)*(1-zeta);
                N(5) = 0.125 * (1-xi)*(1-eta)*(1+zeta);
                N(6) = 0.125 * (1+xi)*(1-eta)*(1+zeta);
                N(7) = 0.125 * (1+xi)*(1+eta)*(1+zeta);
                N(8) = 0.125 * (1-xi)*(1+eta)*(1+zeta);
                
                % Şekil fonksiyonları türevleri (Jacobian için)
                dN_dxi = 0.125 * [
                    -(1-eta)*(1-zeta),  (1-eta)*(1-zeta),  (1+eta)*(1-zeta), -(1+eta)*(1-zeta), ...
                    -(1-eta)*(1+zeta),  (1-eta)*(1+zeta),  (1+eta)*(1+zeta), -(1+eta)*(1+zeta);
                    
                    -(1-xi)*(1-zeta), -(1+xi)*(1-zeta),  (1+xi)*(1-zeta),  (1-xi)*(1-zeta), ...
                    -(1-xi)*(1+zeta), -(1+xi)*(1+zeta),  (1+xi)*(1+zeta),  (1-xi)*(1+zeta);
                    
                    -(1-xi)*(1-eta), -(1+xi)*(1-eta), -(1+xi)*(1+eta), -(1-xi)*(1+eta), ...
                     (1-xi)*(1-eta),  (1+xi)*(1-eta),  (1+xi)*(1+eta),  (1-xi)*(1+eta)
                ];
                
                % Jacobian
                J = dN_dxi * nodes;
                detJ = det(J);
                
                % N matrisi (3x24)
                N_matrix = zeros(3, 24);
                for m = 1:8
                    N_matrix(1, 3*m-2) = N(m);
                    N_matrix(2, 3*m-1) = N(m);
                    N_matrix(3, 3*m)   = N(m);
                end
                
                % Tutarlı kütle matrisine katkı
                Me = Me + rho * (N_matrix' * N_matrix) * detJ * w(i) * w(j) * w(k);
            end
        end
    end
end

%% GLOBAL MATRİSLERİ OLUŞTUR
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
    Ke = hex8_stiffness(E, nu, elem_nodes);
    Me = hex8_consistent_mass(rho, elem_nodes);
    
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

fprintf('Global matrisler oluşturuldu.\n\n');

%% MASS MATRİSİ ANALİZİ
fprintf('=== MASS MATRİSİ ANALİZİ ===\n');

% Temel özellikler
diag_vals = diag(M_global);
off_diag = M_global - diag(diag_vals);

fprintf('Boyut: %dx%d\n', size(M_global, 1), size(M_global, 2));
fprintf('Toplam kütle: %.6e kg\n', trace(M_global)/3);
fprintf('Teorik kütle: %.6e kg\n', rho * Lx * Ly * Lz);
fprintf('Kütle korunumu: %.2e kg\n', abs(trace(M_global)/3 - rho * Lx * Ly * Lz));
fprintf('Simetrik: %s\n', string(isequal(M_global, M_global')));
fprintf('Pozitif tanımlı: %s\n', string(all(eig(M_global) > 1e-20)));
fprintf('Köşegen ortalama: %.2e kg\n', mean(diag_vals));
fprintf('Köşegen std: %.2e kg\n', std(diag_vals));
fprintf('Maks köşegen-dışı: %.2e kg\n', max(abs(off_diag(:))));
fprintf('Doluluk oranı: %.1f%%\n', nnz(off_diag)/numel(off_diag)*100);
fprintf('\n');

%% MODAL ANALİZ
fprintf('Modal analiz yapılıyor...\n');

% ANSYS referans değerleri
ansys_freqs = [3.8e5, 4.4e5, 6.1e5, 7.4e5, 8.0e5, 8.16e5];

% Regularization (rigid body modlar için)
K_reg = K_global + 1e-9 * max(diag(K_global)) * eye(size(K_global));

% Özdeğer problemi çöz
[V, D] = eigs(K_reg, M_global, 20, 'sm');

% Özdeğerleri sırala
lambda = diag(D);
[lambda, idx] = sort(lambda);
frequencies = sqrt(abs(lambda)) / (2*pi);
V = V(:, idx);

%% SONUÇLARI GÖSTER
fprintf('\n=== MODAL ANALİZ SONUÇLARI ===\n');
fprintf('Mod | Frekans (Hz) | ANSYS (Hz) | Fark (%%) | Tip\n');
fprintf('----|--------------|------------|----------|------------\n');

for i = 1:15
    if i <= 6
        freq_type = 'Rigid Body';
        ansys_val = 0;
        error_pct = 0;
    else
        freq_type = 'Esnek Mod';
        ansys_idx = i - 6;
        if ansys_idx <= length(ansys_freqs)
            ansys_val = ansys_freqs(ansys_idx);
            error_pct = abs(frequencies(i) - ansys_val) / ansys_val * 100;
        else
            ansys_val = NaN;
            error_pct = NaN;
        end
    end
    
    if i <= 6
        fprintf('%3d | %12.0f | %10s | %8s | %s\n', i, frequencies(i), '-', '-', freq_type);
    else
        fprintf('%3d | %12.0f | %10.0f | %8.2f | %s\n', i, frequencies(i), ansys_val, error_pct, freq_type);
    end
end

fprintf('\nOrtalama hata (esnek modlar): %.2f%%\n', mean(abs(frequencies(7:12) - ansys_freqs) ./ ansys_freqs * 100));

%% MASS MATRİSİ DETAYLI ÇIKTI
fprintf('\n=== MASS MATRİSİ DETAYLI ÇIKTI ===\n');

fprintf('\nİlk 12x12 Blok (kg):\n');
fprintf('     ');
for j = 1:12
    fprintf('%8d ', j);
end
fprintf('\n     ');
for j = 1:12
    fprintf('---------');
end
fprintf('\n');

for i = 1:12
    fprintf('%3d |', i);
    for j = 1:12
        if abs(M_global(i,j)) > 1e-20
            if i == j
                fprintf('%8.2e*', M_global(i,j));
            else
                fprintf('%8.2e ', M_global(i,j));
            end
        else
            fprintf('%8s ', '0');
        end
    end
    fprintf('\n');
end
fprintf('* = köşegen eleman\n');

fprintf('\nDüğüm Kütle Dağılımı (ilk 8 düğüm):\n');
fprintf('Düğüm |   M_xx (kg)   |   M_yy (kg)   |   M_zz (kg)   |  Toplam (kg)\n');
fprintf('------|---------------|---------------|---------------|-------------\n');

for node = 1:8
    dof_x = (node-1)*3 + 1;
    dof_y = (node-1)*3 + 2;
    dof_z = (node-1)*3 + 3;
    
    m_xx = M_global(dof_x, dof_x);
    m_yy = M_global(dof_y, dof_y);
    m_zz = M_global(dof_z, dof_z);
    total = m_xx + m_yy + m_zz;
    
    fprintf('%5d | %13.2e | %13.2e | %13.2e | %13.2e\n', node, m_xx, m_yy, m_zz, total);
end

%% GÖRSELLEŞTİRME
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

% 2. Mass matrisi yoğunluğu
subplot(2,4,2);
spy(M_global);
title('Mass Matrisi Yoğunluğu');
xlabel('DOF'); ylabel('DOF');

% 3. Köşegen elemanlar
subplot(2,4,3);
plot(diag_vals, 'bo-', 'MarkerSize', 2);
title('Köşegen Elemanlar');
xlabel('DOF'); ylabel('Kütle (kg)');
grid on;

% 4. Frekans karşılaştırması
subplot(2,4,4);
mod_nums = 1:6;
bar_data = [frequencies(7:12)', ansys_freqs(1:6)'];
bar(mod_nums, bar_data);
legend('MATLAB', 'ANSYS', 'Location', 'northwest');
xlabel('Mod No'); ylabel('Frekans (Hz)');
title('Frekans Karşılaştırması');
grid on;

% 5. İlk mod şekli
subplot(2,4,5);
mode_shape = V(:, 7); % İlk esnek mod
scale = 0.1 / max(abs(mode_shape));
deformed_nodes = nodes + scale * reshape(mode_shape, 3, [])';

plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'bo', 'MarkerSize', 2);
hold on;
plot3(deformed_nodes(:,1), deformed_nodes(:,2), deformed_nodes(:,3), 'ro', 'MarkerSize', 3);

for elem_idx = 1:size(elements, 1)
    orig_elem = nodes(elements(elem_idx, :), :);
    def_elem = deformed_nodes(elements(elem_idx, :), :);
    
    for edge = [1,2,3,4,1; 5,6,7,8,5; 1,5; 2,6; 3,7; 4,8]'
        plot3(orig_elem(edge,1), orig_elem(edge,2), orig_elem(edge,3), 'b-', 'LineWidth', 0.5);
        plot3(def_elem(edge,1), def_elem(edge,2), def_elem(edge,3), 'r-', 'LineWidth', 2);
    end
end

title(sprintf('Mod 1: %.0f Hz', frequencies(7)));
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
axis equal; grid on;

% 6. Hata analizi
subplot(2,4,6);
errors = abs(frequencies(7:12) - ansys_freqs) ./ ansys_freqs * 100;
bar(mod_nums, errors);
xlabel('Mod No'); ylabel('Hata (%)');
title('ANSYS-MATLAB Farkı');
grid on;

% 7. Mass matrisi büyüklük haritası
subplot(2,4,7);
imagesc(log10(abs(M_global) + 1e-20));
colorbar;
title('Mass Matrisi (log scale)');
xlabel('DOF'); ylabel('DOF');

% 8. Özdeğer dağılımı
subplot(2,4,8);
semilogy(1:15, abs(lambda(1:15)), 'bo-', 'LineWidth', 2);
xlabel('Mod No'); ylabel('|Özdeğer|');
title('Özdeğer Dağılımı');
grid on;

%% DOSYALARA KAYDET
fprintf('\nSonuçlar kaydediliyor...\n');

% Mass matrisi
csvwrite('mass_matrix_consistent.csv', M_global);

% Frekans sonuçları
freq_results = table((1:20)', frequencies(1:20), 'VariableNames', {'Mod_No', 'Frekans_Hz'});
writetable(freq_results, 'modal_frekanslari.csv');

% Karşılaştırma sonuçları
comparison_results = table((1:6)', frequencies(7:12), ansys_freqs', ...
    abs(frequencies(7:12) - ansys_freqs') ./ ansys_freqs' * 100, ...
    'VariableNames', {'Mod_No', 'MATLAB_Hz', 'ANSYS_Hz', 'Hata_Yuzde'});
writetable(comparison_results, 'ansys_matlab_karsilastirma.csv');

fprintf('Kaydedilen dosyalar:\n');
fprintf('- mass_matrix_consistent.csv\n');
fprintf('- modal_frekanslari.csv\n');
fprintf('- ansys_matlab_karsilastirma.csv\n');

fprintf('\n=== ANALİZ TAMAMLANDI ===\n');
