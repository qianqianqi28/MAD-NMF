clear; close all; clc;
cd('D:\MAD-NMF\AnalySize-master\AnalySize-master');
restoredefaultpath
addpath(genpath('D:\MAD-NMF\AnalySize-master\AnalySize-master'))
rng(42);
colnames = [0.24 0.29 0.35 0.41 0.49 0.58 0.69 0.82 0.98 1.2 1.4 1.6 2 2.3 2.8 ...
     3.3 3.9 4.7 5.5 6.6 7.8 9.3 11.1 13.1 15.6 18.6 22.1 26.3 31.3 37.2 ...
     44.2 52.6 62.5 74.3 88.4 105 125 149 177 210 250 297 354 420 500 ...
     595 707 841 1000 1189 1414 1681 2000];
basis_temp = readtable('Appendix B test datasets.xlsx', 'Sheet', 'Coversand EMs');
true_basis = table2array(basis_temp(1:end, 2:5));
true_basis = true_basis ./ sum(true_basis, 1);
true_basis = true_basis';
nBins = size(true_basis,2);
nEM   = 4;

n_specimens = 200;
% ---- Generate the coversand dataset (noise‑free) --------------------
min_abundance = 0.08;
max_abundance = 0.49;

% Generate random matrix in [min_abundance, max_abundance]
abund_highmix = constrained_composition(n_specimens, nEM, min_abundance, max_abundance);


% sum(abund_highmix, 2)
highly_mixed_samples = abund_highmix * true_basis;   % size n_specimens x nBins
fprintf('Coversand dataset (200 samples, noise‑free):\n');
fprintf('mean abundances: %.3f %.3f %.3f %.3f\n', mean(abund_highmix));
fprintf('abundance ranges: [%.3f, %.3f]%%\n', ...
        min(abund_highmix(:)), max(abund_highmix(:)));
fprintf('sample sum check: sum = %f (should be 1)\n', mean(sum(highly_mixed_samples,2)));


X = highly_mixed_samples;
collabel =  string(colnames);
colnames = 1:length(log(colnames));

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); 
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(colnames, true_basis(1, :), 'b', 'LineWidth', 3)  
hold on
plot(colnames, true_basis(2, :), 'r', 'LineWidth', 3)  
plot(colnames, true_basis(3, :), 'g', 'LineWidth', 3)  
plot(colnames, true_basis(4, :), 'c', 'LineWidth', 3)  
legend('True EM 1', 'True EM 2', 'True EM 3', 'True EM 4', 'FontSize', 16)
xticks(colnames(1:4:end));
xticklabels(collabel(1:4:end))
xlabel('Grain size', 'FontSize', 26, 'FontWeight', 'bold');
ylabel('Volume content', 'FontSize', 26, 'FontWeight', 'bold');
ylim([0 0.35])
exportgraphics(gcf,'plots/coversanddatatrueems.pdf','ContentType','vector');

%% HALS
MaxIter = 1000;
Reps = 1;
r = 4; 

%% min dis
lambdamin = 3.5;
Regs = [100, lambdamin];
disp('Running AnalySize:');
[S_min, A_min, ~] = HALS_NMF(X, r, MaxIter, Reps, Regs);
S_min = S_min';
A_min = A_min';
[S_min, A_min] = matchCol(S_min, true_basis', A_min);
S_min = S_min';
A_min = A_min';

%% max Dis
% lambdamax = 2.5;
lambdamax = 3.5;
Regs = [100, lambdamax];
disp('Running max-dis NMF:');
[S_max, A_max, ~] = MAD_NMF_SimplexProj(X, r, MaxIter, Reps, Regs);
S_max = S_max';
A_max = A_max';
[S_max, A_max] = matchCol(S_max, true_basis', A_max);
S_max = S_max';
A_max = A_max';

%% Figures min
figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); 
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
h1 = plot(colnames, true_basis, 'k-+', 'LineWidth', 3);
hold on
h2 = plot(colnames, S_min(1, :), 'b', 'LineWidth', 3)  
h3 = plot(colnames, S_min(2, :), 'r', 'LineWidth', 3)  
h4 = plot(colnames, S_min(3, :), 'g', 'LineWidth', 3)  
h5 = plot(colnames, S_min(4, :), 'c', 'LineWidth', 3) 
lgd = legend([h1(1), h2(1), h3(1), h4(1), h5(1)], {'True EMs', "Esti. EM 1", "Esti. EM 2", "Esti. EM 3", "Esti. EM 4"}, 'Location','best', FontSize=16)
xticks(colnames(1:4:end));
xticklabels(collabel(1:4:end));
xlabel('Grain size', 'FontSize', 26, 'FontWeight', 'bold');
ylabel('Volume content', 'FontSize', 26, 'FontWeight', 'bold');
ylim([0 0.35])
exportgraphics(gcf,'plots/coversanddataminems.pdf','ContentType','vector');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); 
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
h1 = plot(abund_highmix(:, 1), abund_highmix(:, 1),'Color', 'black', 'LineWidth', 2) 
axis equal
hold on
% plot(abund_highmix(:, 2), abund_highmix(:, 2),'Color', 'black', 'LineWidth', 2) 
h2 = plot(abund_highmix(:, 1), A_min(:, 1), 'bo', 'LineWidth', 2)
h3 = plot(abund_highmix(:, 2), A_min(:, 2),'ro', 'LineWidth', 2)  
h4 = plot(abund_highmix(:, 3), A_min(:, 3),'go', 'LineWidth', 2) 
h5 = plot(abund_highmix(:, 4), A_min(:, 4),'co', 'LineWidth', 2) 

% Add the x-axis label with the mu symbol
xlabel("True Abundances", 'FontSize', 26, 'FontWeight', 'bold')
% Add the y-axis label
ylabel("Estimated Abundances", 'FontSize', 26, 'FontWeight', 'bold')
legend([h1(1), h2(1), h3(1), h4(1), h5(1)], {"1:1", "Specimens: EM1", "Specimens: EM2", "Specimens: EM3", "Specimens: EM4"}, Location="northwest", FontSize=16)
lambdamin_rounded = round(lambdamin, 3);
xlim([0 1])
ylim([0 1])
exportgraphics(gcf, 'plots/coversanddataminabundances.pdf', 'ContentType', 'vector');


%% Figures max
figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); 
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);

h1 = plot(colnames, true_basis, 'k-+', 'LineWidth', 3);
hold on
h2 = plot(colnames, S_max(1, :), 'b', 'LineWidth', 3)  
h3 = plot(colnames, S_max(2, :), 'r', 'LineWidth', 3)  
h4 = plot(colnames, S_max(3, :), 'g', 'LineWidth', 3)  
h5 = plot(colnames, S_max(4, :), 'c', 'LineWidth', 3) 
lgd = legend([h1(1), h2(1), h3(1), h4(1), h5(1)], {'True EMs', "Esti. EM 1", "Esti. EM 2", "Esti. EM 3", "Esti. EM 4"}, 'Location','best', FontSize=16)
xticks(colnames(1:4:end));
xticklabels(collabel(1:4:end))
xlabel('Grain size', 'FontSize', 26, 'FontWeight', 'bold');
ylabel('Volume content', 'FontSize', 26, 'FontWeight', 'bold');
ylim([0 0.35])
exportgraphics(gcf,'plots/coversanddatamaxems.pdf','ContentType','vector');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); 
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
h1 = plot(abund_highmix(:, 1), abund_highmix(:, 1),'Color', 'black', 'LineWidth', 2) 
axis equal
hold on
% plot(abund_highmix(:, 2), abund_highmix(:, 2),'Color', 'black', 'LineWidth', 2) 
h2 = plot(abund_highmix(:, 1), A_max(:, 1), 'bo', 'LineWidth', 2)
h3 = plot(abund_highmix(:, 2), A_max(:, 2),'ro', 'LineWidth', 2)  
h4 = plot(abund_highmix(:, 3), A_max(:, 3),'go', 'LineWidth', 2) 
h5 = plot(abund_highmix(:, 4), A_max(:, 4),'co', 'LineWidth', 2) 

% Add the x-axis label with the mu symbol
xlabel("True Abundances", 'FontSize', 26, 'FontWeight', 'bold')
% Add the y-axis label
ylabel("Estimated Abundances", 'FontSize', 26, 'FontWeight', 'bold')
legend([h1(1), h2(1), h3(1), h4(1), h5(1)], {"1:1", "Specimens: EM1", "Specimens: EM2", "Specimens: EM3", "Specimens: EM4"}, Location="northwest", FontSize=16)
lgd.FontSize = 26;
lambdamax_rounded = round(lambdamax, 3);
xlim([0 1])
ylim([0 1])
exportgraphics(gcf,'plots/coversanddatamaxabundances.pdf','ContentType','vector');


%% angle
angles_sample_min = NaN(n_specimens, 1);
angles_sample_max = NaN(n_specimens, 1);

for i = 1:n_specimens
    angles_sample_min(i) = spectral_angle(abund_highmix(i, :)', A_min(i, :)');
    angles_sample_max(i) = spectral_angle(abund_highmix(i, :)', A_max(i, :)');

end

angles_class_min = NaN(r, 1);
angles_class_max = NaN(r, 1);
for i = 1:r
    angles_class_min(i) = spectral_angle(true_basis(i,:)', S_min(i,:)');
    angles_class_max(i) = spectral_angle(true_basis(i,:)', S_max(i,:)');

end
angles_class_mean = [mean(angles_class_min * 180/pi), mean(angles_class_max * 180/pi)]

angles_sample_mean = [mean(angles_sample_min * 180/pi), mean(angles_sample_max * 180/pi)]


[K, n] = size(S_min);

CK = eye(K) - ones(K)/K;
Cn = eye(n) - ones(n)/n;
dis = [norm(CK * S_min * Cn, 'fro')^2, norm(CK * S_max * Cn, 'fro')^2]

dis = [norm(CK * S_min, 'fro')^2, norm(CK * S_max, 'fro')^2]
