clear all; close all; clc;
cd('D:\MAD-NMF\AnalySize-master\AnalySize-master');
restoredefaultpath
addpath(genpath('D:\MAD-NMF\AnalySize-master\AnalySize-master'))
rng(42); % Set seed for reproducibility
%% Parameters
n_specimens = 99;          % Number of specimens
min_abundance = 0.13;  % Minimum abundance (13%)
max_abundance = 0.87;  % Maximum abundance (87%)
% Grain size vector
T = readtable('D:\MAD-NMF\AnalySize-master\AnalySize-master\Example_Data\Example_Data_2.xlsx','Sheet','True_Data');
grain_size = T{:,1};
collabel =  string(grain_size);
%% Define Lognormal Distributions for Sources
% Source 1
% Parameters for Source 1
sigma1 = 0.55;
mu1 = 2 + sigma1^2;  

% Parameters for Source 2
sigma2 = 1;
mu2 = 6+sigma2^2;

% Calculate lognormal probability density functions
source1 = (1 ./ (grain_size * sigma1 * sqrt(2 * pi))) .* ...
          exp(-((log(grain_size) - mu1).^2) / (2 * sigma1^2));
source1 = source1/sum(source1);
source2 = (1 ./ (grain_size * sigma2 * sqrt(2 * pi))) .* ...
          exp(-((log(grain_size) - mu2).^2) / (2 * sigma2^2));
source2 = source2/sum(source2);
endmember_matrix = NaN(2, length(source1));
endmember_matrix(1, :) = source1;
endmember_matrix(2, :) = source2;

%% Generate True Abundances
% Initialize abundance matrix
abundance_matrix = NaN(n_specimens, 2);
% Generate random abundances for Source 1 in range [0.13, 0.87]
source1_abundance = min_abundance + (max_abundance - min_abundance) * rand(n_specimens, 1);
% Set abundances for Source 2 as 1 - Source 1 abundance
source2_abundance = 1 - source1_abundance;

% Combine into abundance matrix
[~, sortIndex] = sort(source1_abundance, "descend");

abundance_matrix(:, 1) = source1_abundance(sortIndex);
abundance_matrix(:, 2) = source2_abundance(sortIndex);

%% Generate Mixed Distributions for Each Specimen
mixed_distributions = abundance_matrix * endmember_matrix;

% add noise
noise_factor = 1 + 0.01*randn(n_specimens, length(grain_size));
noisy_samples = mixed_distributions .* noise_factor;

% Avoid negative values
noisy_samples(noisy_samples < 0) = 0;
% Renormalise to 100% per sample
noisy_mixed_distributions = noisy_samples ./ sum(noisy_samples,2);

X = noisy_mixed_distributions;

MaxIter = 1000;
Reps = 1;

r = 2; 

disp('Running AnalySize:');
lambdamin = 1;
Regs = [5, lambdamin];
[S_min, A_min, ~] = HALS_NMF(X, r, MaxIter, Reps, Regs);
S_min = S_min';
A_min = A_min';
[S_min, A_min] = matchCol(S_min, endmember_matrix', A_min);
S_min = S_min';
A_min = A_min';

disp('Running max-dis NMF:');
lambdamax = 1;
Regs = [5, lambdamax];
[S_max, A_max, ~] = MAD_NMF_SimplexProj(X, r, MaxIter, Reps, Regs);

S_max = S_max';
A_max = A_max';
[S_max, A_max] = matchCol(S_max, endmember_matrix', A_max);
S_max = S_max';
A_max = A_max';

%% FIGURES
% figure('Position', [100, 100, 800, 800]);
% set(gcf, 'PaperSize', [21.0 21.0]); 
% set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
% plot(1:length(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 3)  
% hold on
% plot(1:length(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 3)  
% xticks(1:3:length(grain_size));
% xticklabels(collabel(1:3:end));
% xlabel('Grain size', 'FontSize', 26, 'FontWeight', 'bold');
% ylabel('Volume content', 'FontSize', 26, 'FontWeight', 'bold');
% legend('True EM 1', 'True EM 2', 'FontSize', 16)
% exportgraphics(gcf,'plots/noisytwodatatrueems.pdf','ContentType','vector');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); 
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
h1 = plot(1:length(grain_size), endmember_matrix, 'k-+', 'LineWidth', 3)  % First column in red
hold on
h2 = plot(1:length(grain_size), S_min(1,:), 'b', 'LineWidth', 3)  % First column in red
h3 = plot(1:length(grain_size), S_min(2,:), 'r', 'LineWidth', 3)  % Second column in blue
xticks(1:3:length(grain_size));
xticklabels(collabel(1:3:end));
% Add the x-axis label with the mu symbol
xlabel("Grain size", 'FontSize', 26, 'FontWeight', 'bold');
% Add the y-axis label
ylabel("Volume content", 'FontSize', 26, 'FontWeight', 'bold');
legend([h1(1), h2(1), h3(1)], {'True EMs', "Esti. EM 1", "Esti. EM 2"}, 'FontSize', 16)
exportgraphics(gcf,'plots/noisytwodataminems.pdf','ContentType','vector');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); 
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
h1 = plot(abundance_matrix(:, 1), abundance_matrix(:, 1),'Color', 'black', 'LineWidth', 3)  % Second column in blue
axis equal
hold on
% plot(abundance_matrix(:, 2), abundance_matrix(:, 2),'Color', 'black', 'LineWidth', 3)  % Second column in blue
h2 = plot(abundance_matrix(:, 1), A_min(:, 1), 'bo', 'LineWidth', 3)  % First column in red
h3 = plot(abundance_matrix(:, 2), A_min(:, 2),'ro', 'LineWidth', 3)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("True Abundances", 'FontSize', 26, 'FontWeight', 'bold');
% Add the y-axis label
ylabel("Estimated Abundances", 'FontSize', 26, 'FontWeight', 'bold');
legend([h1(1), h2(1), h3(1)], {"1:1", "Specimens: EM1", "Specimens: EM2"}, Location="northwest", FontSize=16)
lambdamin_rounded = round(lambdamin, 3);
xlim([0 1])
ylim([0 1])
exportgraphics(gcf,'plots/noisytwodataminabundances.pdf','ContentType','vector');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
h1 = plot(1:length(grain_size), endmember_matrix, 'k-+', 'LineWidth', 3)  % First column in red
hold on
h2 = plot(1:length(grain_size), S_max(1, :), 'b', 'LineWidth', 3)  % First column in red
h3 = plot(1:length(grain_size), S_max(2, :), 'r', 'LineWidth', 3)  % Second column in blue
xticks(1:3:length(grain_size));
xticklabels(collabel(1:3:end));
% Add the x-axis label with the mu symbol
xlabel("Grain size", 'FontSize', 26, 'FontWeight', 'bold');
% Add the y-axis label
ylabel("Volume content", 'FontSize', 26, 'FontWeight', 'bold');
legend([h1(1), h2(1), h3(1)], {'True EMs', "Esti. EM 1", "Esti. EM 2"}, 'FontSize', 16)
exportgraphics(gcf,'plots/noisytwodatamaxems.pdf','ContentType','vector');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]);
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
h1 = plot(abundance_matrix(:, 1), abundance_matrix(:, 1),'Color', 'black', 'LineWidth', 3)  % Second column in blue
axis equal
hold on
% plot(abundance_matrix(:, 2), abundance_matrix(:, 2),'Color', 'black', 'LineWidth', 3)  % Second column in blue
h2 = plot(abundance_matrix(:, 1), A_max(:, 1), 'bo', 'LineWidth', 3)  % First column in red
h3 = plot(abundance_matrix(:, 2), A_max(:, 2),'ro', 'LineWidth', 3)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("True Abundances", 'FontSize', 26, 'FontWeight', 'bold');
% Add the y-axis label
ylabel("Estimated Abundances", 'FontSize', 26, 'FontWeight', 'bold');
legend([h1(1), h2(1), h3(1)], {"1:1", "Specimens: EM1", "Specimens: EM2"}, Location="northwest", FontSize=16)
xlim([0 1])
ylim([0 1])
exportgraphics(gcf,'plots/noisytwodatamaxabundances.pdf','ContentType','vector');

%% angle
angles_sample_min = NaN(n_specimens, 1);
angles_sample_max = NaN(n_specimens, 1);

for i = 1:n_specimens
    angles_sample_min(i) = spectral_angle(abundance_matrix(i, :)', A_min(i, :)');
    angles_sample_max(i) = spectral_angle(abundance_matrix(i, :)', A_max(i, :)');

end

angles_class_min = NaN(r, 1);
angles_class_max = NaN(r, 1);
for i = 1:r
    angles_class_min(i) = spectral_angle(endmember_matrix(i,:)', S_min(i,:)');
    angles_class_max(i) = spectral_angle(endmember_matrix(i,:)', S_max(i,:)');
end
angles_class_mean = [mean(angles_class_min * 180/pi), mean(angles_class_max * 180/pi)]
angles_sample_mean = [mean(angles_sample_min * 180/pi), mean(angles_sample_max * 180/pi)]


[K, n] = size(S_min);

CK = eye(K) - ones(K)/K;
Cn = eye(n) - ones(n)/n;
dis = [norm(CK * S_min * Cn, 'fro')^2, norm(CK * S_max * Cn, 'fro')^2]

dis = [norm(CK * S_min, 'fro')^2, norm(CK * S_max, 'fro')^2]
