%% REFERANSA GÖRE DÜZELTİLMİŞ K MATRİSİ
fprintf('=== REFERANSA GÖRE K MATRİSİ DÜZELTİLİYOR ===\n\n');

%% Düzeltilmiş Hex8 Sertlik Matrisi Fonksiyonu
function [Ke, Be, D] = hex8_stiffness_corrected(E, nu, nodes)
    % Malzeme matrisi D (referanstaki gibi)
    C = zeros(6,6);
    
    % İzotropik malzeme için elastisite matrisi
    lambda = E*nu/((1+nu)*(1-2*nu));
    mu = E/(2*(1+nu));
    
    C(1:3,1:3) = lambda;
    for i=1:3
        C(i,i) = C(i,i) + 2*mu;
    end
    for i=4:6
        C(i,i) = mu;
    end
    
    D = C;
    
    % Gauss integrasyon noktaları ve ağırlıkları
    a = 1/sqrt(3);
    gp = [-a, a];
    w = [1, 1];
    
    Ke = zeros(24,24);
    Be = zeros(6,24);
    
    for i1=1:2
        for i2=1:2
            for i3=1:2
                xi = gp(i1);
                eta = gp(i2);
                zeta = gp(i3);
                
                % Şekil fonksiyonları ve türevleri
                N = zeros(1,8);
                dNdxi = zeros(3,8);
                
                % Şekil fonksiyonları
                N(1) = (1-xi)*(1-eta)*(1-zeta)/8;
                N(2) = (1+xi)*(1-eta)*(1-zeta)/8;
                N(3) = (1+xi)*(1+eta)*(1-zeta)/8;
                N(4) = (1-xi)*(1+eta)*(1-zeta)/8;
                N(5) = (1-xi)*(1-eta)*(1+zeta)/8;
                N(6) = (1+xi)*(1-eta)*(1+zeta)/8;
                N(7) = (1+xi)*(1+eta)*(1+zeta)/8;
                N(8) = (1-xi)*(1+eta)*(1+zeta)/8;
                
                % Doğal koordinatlarda türevler
                dNdxi(1,1) = -(1-eta)*(1-zeta)/8;
                dNdxi(1,2) =  (1-eta)*(1-zeta)/8;
                dNdxi(1,3) =  (1+eta)*(1-zeta)/8;
                dNdxi(1,4) = -(1+eta)*(1-zeta)/8;
                dNdxi(1,5) = -(1-eta)*(1+zeta)/8;
                dNdxi(1,6) =  (1-eta)*(1+zeta)/8;
                dNdxi(1,7) =  (1+eta)*(1+zeta)/8;
                dNdxi(1,8) = -(1+eta)*(1+zeta)/8;
                
                dNdxi(2,1) = -(1-xi)*(1-zeta)/8;
                dNdxi(2,2) = -(1+xi)*(1-zeta)/8;
                dNdxi(2,3) =  (1+xi)*(1-zeta)/8;
                dNdxi(2,4) =  (1-xi)*(1-zeta)/8;
                dNdxi(2,5) = -(1-xi)*(1+zeta)/8;
                dNdxi(2,6) = -(1+xi)*(1+zeta)/8;
                dNdxi(2,7) =  (1+xi)*(1+zeta)/8;
                dNdxi(2,8) =  (1-xi)*(1+zeta)/8;
                
                dNdxi(3,1) = -(1-xi)*(1-eta)/8;
                dNdxi(3,2) = -(1+xi)*(1-eta)/8;
                dNdxi(3,3) = -(1+xi)*(1+eta)/8;
                dNdxi(3,4) = -(1-xi)*(1+eta)/8;
                dNdxi(3,5) =  (1-xi)*(1-eta)/8;
                dNdxi(3,6) =  (1+xi)*(1-eta)/8;
                dNdxi(3,7) =  (1+xi)*(1+eta)/8;
                dNdxi(3,8) =  (1-xi)*(1+eta)/8;
                
                % Jacobian matrisi
                J = zeros(3,3);
                for i=1:8
                    J(1,1) = J(1,1) + dNdxi(1,i)*nodes(i,1);
                    J(1,2) = J(1,2) + dNdxi(1,i)*nodes(i,2);
                    J(1,3) = J(1,3) + dNdxi(1,i)*nodes(i,3);
                    J(2,1) = J(2,1) + dNdxi(2,i)*nodes(i,1);
                    J(2,2) = J(2,2) + dNdxi(2,i)*nodes(i,2);
                    J(2,3) = J(2,3) + dNdxi(2,i)*nodes(i,3);
                    J(3,1) = J(3,1) + dNdxi(3,i)*nodes(i,1);
                    J(3,2) = J(3,2) + dNdxi(3,i)*nodes(i,2);
                    J(3,3) = J(3,3) + dNdxi(3,i)*nodes(i,3);
                end
                
                detJ = det(J);
                Jinv = inv(J);
                
                % Global türevler
                dNdx = zeros(3,8);
                for i=1:8
                    dNdx(:,i) = Jinv * dNdxi(:,i);
                end
                
                % B matrisi (strain-displacement)
                B = zeros(6,24);
                for i=1:8
                    B(1,3*i-2) = dNdx(1,i);
                    B(2,3*i-1) = dNdx(2,i);
                    B(3,3*i)   = dNdx(3,i);
                    B(4,3*i-2) = dNdx(2,i);
                    B(4,3*i-1) = dNdx(1,i);
                    B(5,3*i-1) = dNdx(3,i);
                    B(5,3*i)   = dNdx(2,i);
                    B(6,3*i-2) = dNdx(3,i);
                    B(6,3*i)   = dNdx(1,i);
                end
                
                % Sertlik matrisine katkı
                Ke = Ke + B' * D * B * detJ * w(i1) * w(i2) * w(i3);
                Be = B; % Son B matrisini sakla
            end
        end
    end
end

%% Düzeltilmiş Global K Matrisini Oluştur
fprintf('Düzeltilmiş K matrisi oluşturuluyor...\n');

K_global_corrected = zeros(num_nodes*3, num_nodes*3);

for elem_idx = 1:num_elements
    if mod(elem_idx, 10) == 0
        fprintf('  Eleman %d/%d işleniyor...\n', elem_idx, num_elements);
    end
    
    % Eleman düğüm koordinatları
    elem_nodes = nodes(elements(elem_idx, :), :);
    
    % Düzeltilmiş sertlik matrisi
    [Ke, Be, De] = hex8_stiffness_corrected(E, nu, elem_nodes);
    
    % Global matrise assemble et
    for i = 1:8
        for j = 1:8
            node_i = elements(elem_idx, i);
            node_j = elements(elem_idx, j);
            
            dof_i = (node_i-1)*3 + (1:3);
            dof_j = (node_j-1)*3 + (1:3);
            
            K_global_corrected(dof_i, dof_j) = K_global_corrected(dof_i, dof_j) + Ke(3*i-2:3*i, 3*j-2:3*j);
        end
    end
end

fprintf('Düzeltilmiş K matrisi oluşturuldu.\n');

%% K MATRİSİ DOĞRULAMA
fprintf('\n=== DÜZELTİLMİŞ K MATRİSİ DOĞRULAMA ===\n');

% 1. Simetri kontrolü
sym_error = norm(K_global_corrected - K_global_corrected', 'fro') / norm(K_global_corrected, 'fro');
fprintf('Simetri hatası: %.2e\n', sym_error);

% 2. Özdeğer analizi
eig_vals = eig(K_global_corrected);
num_zero_eig = sum(abs(eig_vals) < 1e-10);
num_neg_eig = sum(eig_vals < -1e-10);

fprintf('Sıfır özdeğer sayısı (rigid body): %d/6\n', num_zero_eig);
fprintf('Negatif özdeğer sayısı: %d\n', num_neg_eig);

% 3. Rigid body mod testi
fprintf('Rigid body mod testi...\n');
rigid_errors = zeros(6,1);
for i = 1:6
    rigid_errors(i) = norm(K_global_corrected * rigid_modes(:,i));
end
fprintf('Maksimum rigid body hatası: %.2e\n', max(rigid_errors));

%% 10N KUVVET İLE DEPLASMAN TESTİ (DÜZELTİLMİŞ K MATRİSİ İLE)
fprintf('\n=== 10N KUVVET TESTİ (DÜZELTİLMİŞ) ===\n');

% Kuvvet vektörü - X yönünde çekme
F_test = zeros(num_nodes*3, 1);
for i = 1:num_nodes
    if nodes(i,1) == Lx  % x=Lx düzlemindeki düğümler
        F_test(3*i-2) = 10 / (num_nodes_y * num_nodes_z);
    end
end

% Sınır şartları - x=0 düzlemi sabit
constrained_dofs = [];
for i = 1:num_nodes
    if nodes(i,1) == 0
        constrained_dofs = [constrained_dofs, 3*i-2, 3*i-1, 3*i];
    end
end
constrained_dofs = unique(constrained_dofs);

% Kondansasyon
free_dofs = setdiff(1:num_nodes*3, constrained_dofs);
K_cond = K_global_corrected(free_dofs, free_dofs);
F_cond = F_test(free_dofs);

% Deplasmanları çöz
u_corrected = zeros(num_nodes*3, 1);
u_cond = K_cond \ F_cond;
u_corrected(free_dofs) = u_cond;

% Maksimum deplasman
max_u_corrected = max(abs(u_corrected));

% Teorik deplasman
A_cross = Ly * Lz;
delta_theoretical = (10 * Lx) / (A_cross * E);

fprintf('Teorik deplasman: %.2e m\n', delta_theoretical);
fprintf('Hesaplanan deplasman: %.2e m\n', max_u_corrected);
fprintf('Hata: %.1f%%\n', abs(max_u_corrected - delta_theoretical)/delta_theoretical * 100);

%% ESKİ VE YENİ K MATRİSLERİNİ KARŞILAŞTIR
fprintf('\n=== ESKİ vs YENİ K MATRİSİ KARŞILAŞTIRMASI ===\n');

% Norm farkı
K_diff_norm = norm(K_global_corrected - K_global, 'fro') / norm(K_global_corrected, 'fro');
fprintf('K matrisi norm farkı: %.2f%%\n', K_diff_norm * 100);

% Deplasman farkı
u_old = zeros(num_nodes*3, 1);
if exist('u_x', 'var')
    u_old = u_x;
end

if norm(u_old) > 0
    u_diff = norm(u_corrected - u_old) / norm(u_corrected) * 100;
    fprintf('Deplasman farkı: %.1f%%\n', u_diff);
end

%% MODAL ANALİZ TEKRARI (DÜZELTİLMİŞ K MATRİSİ İLE)
fprintf('\n=== DÜZELTİLMİŞ MODAL ANALİZ ===\n');

K_reg = K_global_corrected + 1e-12 * eye(size(K_global_corrected));
[V_corr, D_corr] = eigs(K_reg, M_global, 20, 'sm');
lambda_corr = diag(D_corr);
[lambda_corr, idx_corr] = sort(lambda_corr);
frequencies_corr = sqrt(abs(lambda_corr)) / (2*pi);

fprintf('\nDÜZELTİLMİŞ FREKANSLAR:\n');
fprintf('Mod | Önceki (Hz) | Düzeltilmiş (Hz) | ANSYS (Hz) | Fark (%%)\n');
fprintf('----|--------------|-------------------|------------|----------\n');

for i = 7:12
    ansys_idx = i - 6;
    if ansys_idx <= length(ansys_freqs)
        error_old = abs(frequencies(i) - ansys_freqs(ansys_idx)) / ansys_freqs(ansys_idx) * 100;
        error_new = abs(frequencies_corr(i) - ansys_freqs(ansys_idx)) / ansys_freqs(ansys_idx) * 100;
        
        fprintf('%3d | %11.0f | %16.0f | %9.0f | %7.2f\n', ...
                i-6, frequencies(i), frequencies_corr(i), ansys_freqs(ansys_idx), error_new);
    end
end

%% GÖRSELLEŞTİRME
figure('Position', [100, 100, 1200, 800]);

% 1. Deplasman karşılaştırması
subplot(2,3,1);
if exist('u_x', 'var')
    plot(1:num_nodes, abs(u_x(1:3:end)), 'ro-', 1:num_nodes, abs(u_corrected(1:3:end)), 'bx-');
    legend('Eski K', 'Düzeltilmiş K', 'Location', 'northwest');
else
    plot(1:num_nodes, abs(u_corrected(1:3:end)), 'bx-');
    legend('Düzeltilmiş K', 'Location', 'northwest');
end
xlabel('Düğüm Numarası');
ylabel('X Deplasmanı (m)');
title('Deplasman Karşılaştırması');
grid on;

% 2. Frekans karşılaştırması
subplot(2,3,2);
mod_nums = 1:6;
if exist('frequencies', 'var')
    bar_data = [frequencies(7:12)', frequencies_corr(7:12)', ansys_freqs(1:6)'];
    bar(mod_nums, bar_data);
    legend('Eski', 'Düzeltilmiş', 'ANSYS', 'Location', 'northwest');
else
    bar_data = [frequencies_corr(7:12)', ansys_freqs(1:6)'];
    bar(mod_nums, bar_data);
    legend('Düzeltilmiş', 'ANSYS', 'Location', 'northwest');
end
xlabel('Mod No');
ylabel('Frekans (Hz)');
title('Frekans Karşılaştırması');
grid on;

% 3. K matrisi farkı
subplot(2,3,3);
imagesc(log10(abs(K_global_corrected - K_global) + 1e-20));
colorbar;
title('K Matrisi Farkı (log scale)');
xlabel('DOF'); ylabel('DOF');

% 4. Düzeltilmiş deplasman dağılımı
subplot(2,3,4);
displacement_mag = sqrt(u_corrected(1:3:end).^2 + u_corrected(2:3:end).^2 + u_corrected(3:3:end).^2);
scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 50, displacement_mag, 'filled');
colorbar;
title('Düzeltilmiş Deplasman (m)');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
axis equal;

% 5. Hata analizi
subplot(2,3,5);
if exist('frequencies', 'var')
    errors_old = abs(frequencies(7:12) - ansys_freqs) ./ ansys_freqs * 100;
    errors_new = abs(frequencies_corr(7:12) - ansys_freqs) ./ ansys_freqs * 100;
    bar([errors_old; errors_new]');
    legend('Eski', 'Düzeltilmiş', 'Location', 'northwest');
else
    errors_new = abs(frequencies_corr(7:12) - ansys_freqs) ./ ansys_freqs * 100;
    bar(errors_new);
    legend('Düzeltilmiş', 'Location', 'northwest');
end
xlabel('Mod No');
ylabel('Hata (%)');
title('ANSYS Hata Karşılaştırması');
grid on;

% 6. Teorik vs Hesaplanan
subplot(2,3,6);
theoretical_disp = (10 * Lx) / (Ly * Lz * E) * ones(num_nodes,1);
calculated_disp = abs(u_corrected(1:3:end));
plot(1:num_nodes, theoretical_disp, 'g-', 1:num_nodes, calculated_disp, 'b-');
legend('Teorik', 'Hesaplanan', 'Location', 'northwest');
xlabel('Düğüm No');
ylabel('Deplasman (m)');
title('Teorik vs Hesaplanan Deplasman');
grid on;

%% SONUÇLARI KAYDET
fprintf('\nDüzeltilmiş sonuçlar kaydediliyor...\n');

csvwrite('K_matrix_corrected.csv', K_global_corrected);

% Deplasman sonuçları
corr_disp_data = [(1:num_nodes)', nodes, reshape(u_corrected, 3, [])'];
header = {'Node', 'X', 'Y', 'Z', 'U_x', 'U_y', 'U_z'};
fid = fopen('displacement_corrected.csv', 'w');
fprintf(fid, '%s,', header{1:end-1});
fprintf(fid, '%s\n', header{end});
fclose(fid);
dlmwrite('displacement_corrected.csv', corr_disp_data, '-append');

% Frekans sonuçları
freq_corr_results = table((1:20)', frequencies_corr(1:20), ...
                         'VariableNames', {'Mod_No', 'Frekans_Hz_Corrected'});
writetable(freq_corr_results, 'frequencies_corrected.csv');

fprintf('Kaydedilen dosyalar:\n');
fprintf('- K_matrix_corrected.csv\n');
fprintf('- displacement_corrected.csv\n');
fprintf('- frequencies_corrected.csv\n');

fprintf('\n=== K MATRİSİ DÜZELTME TAMAMLANDI ===\n');
