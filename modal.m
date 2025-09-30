%% TUTARLI KÜTLE MATRİSİ (Consistent Mass Matrix) Hesaplama
fprintf('Tutarlı kütle matrisi hesaplanıyor...\n');

M_global_consistent = zeros(num_nodes*3, num_nodes*3);

%% Tutarlı Kütle Matrisi Fonksiyonu
function Me = hex8_consistent_mass(rho, nodes)
    % 8 noktalı Gauss integrasyonu
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
                
                % Jacobian matrisi
                J = dN_dxi * nodes;
                detJ = det(J);
                
                % N matrisi (3x24) - her şekil fonksiyonu 3 DOF'a katkı yapar
                N_matrix = zeros(3, 24);
                for m = 1:8
                    N_matrix(1, 3*m-2) = N(m);  % x DOF
                    N_matrix(2, 3*m-1) = N(m);  % y DOF
                    N_matrix(3, 3*m)   = N(m);  % z DOF
                end
                
                % Tutarlı kütle matrisine katkı: ρ * Nᵀ * N * detJ
                Me = Me + rho * (N_matrix' * N_matrix) * detJ * w(i) * w(j) * w(k);
            end
        end
    end
end

%% Global Tutarlı Kütle Matrisini Oluştur
for elem_idx = 1:num_elements
    if mod(elem_idx, 10) == 0
        fprintf('  Eleman %d/%d işleniyor...\n', elem_idx, num_elements);
    end
    
    % Eleman düğüm koordinatları
    elem_nodes = nodes(elements(elem_idx, :), :);
    
    % Tutarlı kütle matrisi
    Me = hex8_consistent_mass(rho, elem_nodes);
    
    % Global matrise assemble et
    for i = 1:8
        for j = 1:8
            node_i = elements(elem_idx, i);
            node_j = elements(elem_idx, j);
            
            dof_i = (node_i-1)*3 + (1:3);
            dof_j = (node_j-1)*3 + (1:3);
            
            M_global_consistent(dof_i, dof_j) = M_global_consistent(dof_i, dof_j) + Me(3*i-2:3*i, 3*j-2:3*j);
        end
    end
end

fprintf('Tutarlı kütle matrisi oluşturuldu.\n');

%% Tutarlı Kütle Matrisini Analiz Et
fprintf('\n=== TUTARLI KÜTLE MATRİSİ ANALİZİ ===\n');

fprintf('Tutarlı Kütle Matrisi Özellikleri:\n');
fprintf('  Boyut: %dx%d\n', size(M_global_consistent, 1), size(M_global_consistent, 2));
fprintf('  Toplam kütle: %.6e kg\n', trace(M_global_consistent)/3);
fprintf('  Teorik kütle: %.6e kg\n', rho * Lx * Ly * Lz);
fprintf('  Maksimum değer: %.6e\n', max(M_global_consistent(:)));
fprintf('  Minimum değer: %.6e\n', min(M_global_consistent(M_global_consistent > 0)));
fprintf('  Simetrik mi?: %s\n', string(isequal(M_global_consistent, M_global_consistent')));
fprintf('  Pozitif tanımlı mi?: %s\n', string(all(eig(M_global_consistent) > 1e-20)));

%% İki Kütle Matrisini Karşılaştır
fprintf('\n=== KONSANTRE vs TUTARLI KÜTLE MATRİSİ KARŞILAŞTIRMASI ===\n');

fprintf('Toplam Kütle Karşılaştırması:\n');
fprintf('  Konsantre kütle: %.6e kg\n', trace(M_global)/3);
fprintf('  Tutarlı kütle:   %.6e kg\n', trace(M_global_consistent)/3);
fprintf('  Teorik kütle:    %.6e kg\n', rho * Lx * Ly * Lz);

fprintf('\nKöşegen Eleman Karşılaştırması:\n');
diag_lumped = diag(M_global);
diag_consistent = diag(M_global_consistent);

fprintf('  Konsantre - Ortalama: %.2e, Std: %.2e\n', mean(diag_lumped), std(diag_lumped));
fprintf('  Tutarlı   - Ortalama: %.2e, Std: %.2e\n', mean(diag_consistent), std(diag_consistent));

%% Tutarlı Kütle Matrisinin İlk 12x12 Blokunu Göster
fprintf('\nTutarlı Kütle Matrisi - İlk 12x12 Blok (kg):\n');
fprintf('DOF\\DOF');
for j = 1:12
    fprintf('%10d', j);
end
fprintf('\n');
fprintf('------');
for j = 1:12
    fprintf('----------');
end
fprintf('\n');

for i = 1:12
    fprintf('%4d  |', i);
    for j = 1:12
        if abs(M_global_consistent(i,j)) > 1e-20
            fprintf('%10.2e', M_global_consistent(i,j));
        else
            fprintf('%10s', '0');
        end
    end
    fprintf('\n');
end

%% Modal Analizi Tutarlı Kütle Matrisi ile Tekrarla
fprintf('\nModal analiz (tutarlı kütle matrisi ile) yapılıyor...\n');

K_reg = K_global + 1e-9 * max(diag(K_global)) * eye(size(K_global));

% Özdeğer problemi çöz
[V_consistent, D_consistent] = eigs(K_reg, M_global_consistent, 20, 'sm');

% Özdeğerleri sırala
lambda_consistent = diag(D_consistent);
[lambda_consistent, idx_consistent] = sort(lambda_consistent);
frequencies_consistent = sqrt(abs(lambda_consistent)) / (2*pi);
V_consistent = V_consistent(:, idx_consistent);

%% Frekans Karşılaştırması
fprintf('\n=== FREKANS KARŞILAŞTIRMASI ===\n');
fprintf('Mod | Konsantre(Hz) | Tutarlı(Hz)  | Fark(%%) | ANSYS(Hz)\n');
fprintf('----|---------------|--------------|---------|-----------\n');

for i = 1:min(6, length(frequencies_consistent)-6)
    freq_lumped = frequencies(6+i);
    freq_consistent = frequencies_consistent(6+i);
    freq_ansys = ansys_freqs(i);
    
    error_percent = abs(freq_consistent - freq_lumped) / freq_lumped * 100;
    
    fprintf('%3d | %12.0f | %11.0f | %7.1f | %9.0f\n', ...
            i, freq_lumped, freq_consistent, error_percent, freq_ansys);
end

%% Tutarlı Kütle Matrisini Kaydet
fprintf('\nTutarlı kütle matrisi kaydediliyor...\n');

csvwrite('mass_matrix_consistent_full.csv', M_global_consistent);
writematrix(M_global_consistent, 'mass_matrix_consistent_formatted.csv');

% Köşegen ve köşegen dışı eleman analizi
fprintf('\nTutarlı Kütle Matrisi Köşegen-Dışı Eleman Analizi:\n');
off_diag_elements = M_global_consistent - diag(diag(M_global_consistent));
max_off_diag = max(abs(off_diag_elements(:)));
fprintf('  Maksimum köşegen-dışı eleman: %.2e\n', max_off_diag);
fprintf('  Köşegen-dışı eleman oranı: %.2f%%\n', ...
    nnz(off_diag_elements) / numel(M_global_consistent) * 100);

%% Görselleştirme
figure('Position', [100, 100, 1500, 1000]);

% 1. Tutarlı kütle matrisi yoğunluğu
subplot(2,4,1);
spy(M_global_consistent);
title('Tutarlı Kütle Matrisi Yoğunluğu');
xlabel('DOF'); ylabel('DOF');

% 2. Köşegen eleman karşılaştırması
subplot(2,4,2);
plot(1:72, diag_lumped, 'ro-', 1:72, diag_consistent, 'bx-', 'MarkerSize', 2);
legend('Konsantre', 'Tutarlı', 'Location', 'best');
title('Köşegen Eleman Karşılaştırması');
xlabel('DOF'); ylabel('Kütle (kg)');
grid on;

% 3. Frekans karşılaştırması
subplot(2,4,3);
mod_numbers = 1:6;
bar_data = [frequencies(7:12)', frequencies_consistent(7:12)', ansys_freqs(1:6)'];
bar(mod_numbers, bar_data);
legend('Konsantre', 'Tutarlı', 'ANSYS', 'Location', 'northwest');
xlabel('Mod Numarası'); ylabel('Frekans (Hz)');
title('Frekans Karşılaştırması');
grid on;

% 4. Tutarlı kütle matrisi büyüklük haritası
subplot(2,4,4);
imagesc(log10(abs(M_global_consistent) + 1e-20));
colorbar;
title('Tutarlı Kütle Matrisi (log scale)');
xlabel('DOF'); ylabel('DOF');

% 5. Köşegen-dışı eleman dağılımı
subplot(2,4,5);
[nz_i, nz_j, nz_val] = find(off_diag_elements);
scatter(nz_i, nz_val, 10, 'filled');
title('Köşegen-Dışı Elemanlar');
xlabel('DOF Numarası'); ylabel('Kütle Değeri (kg)');
grid on;

% 6. Her düğümün toplam kütlesi (tutarlı)
subplot(2,4,6);
node_masses_consistent = zeros(num_nodes, 1);
for node = 1:num_nodes
    dof_start = (node-1)*3 + 1;
    node_masses_consistent(node) = sum(diag_consistent(dof_start:dof_start+2));
end
bar(node_masses_consistent);
title('Tutarlı - Düğüm Kütleleri');
xlabel('Düğüm Numarası'); ylabel('Kütle (kg)');
grid on;

% 7. Kütle korunumu hatası
subplot(2,4,7);
mass_errors = abs(node_masses_consistent - node_mass_theoretical);
bar(mass_errors);
title('Düğüm Kütlesi Hataları');
xlabel('Düğüm Numarası'); ylabel('Hata (kg)');
grid on;

% 8. Modal etkin kütle
subplot(2,4,8);
effective_masses = zeros(6,1);
for i = 1:6
    mode = V_consistent(:, 6+i);
    effective_masses(i) = mode' * M_global_consistent * mode;
end
bar(effective_masses);
title('Modal Etkin Kütleler');
xlabel('Mod Numarası'); ylabel('Etkin Kütle (kg)');
grid on;

fprintf('\nTutarlı kütle matrisi analizi tamamlandı.\n');
