clear;clc;

%% 1. Load data, surface, color
load_demo
load_surface
load_color

%% 2. Surface plot
%  2.0 Background brain
bkg_brain = ones(32492,1);
plot_hemispheres(bkg_brain, surf_lh, 'views', 'lm'); colormap(hot); 
set(findall(gcf,'type','axes'),'CLim',[1 2]);

%  2.1 Individual participant Voxel_MidPET_hmz
load('/host/verges/tank/data/mo/demo/01Voxel_data_hmz.mat')
indi_MidPET_hmz = Voxel_all_MidPET_hmz(1,1:32492)';
plot_hemispheres(zscore(indi_MidPET_hmz)', surf_lh, 'views', 'l'); colormap(flipud(blue2red))
set(findall(gcf,'type','axes'),'CLim',[-3 3]);

%  2.2 Individual participant Voxel_MidPET_hmz_flip_asym_wPDB
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')
indi_MidPET_hmz_flip_asym_wPDB = Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(1,1:32492);
plot_hemispheres(zscore(indi_MidPET_hmz_flip_asym_wPDB)', surf_lh, 'views', 'l'); colormap(flipud(blue2red))
set(findall(gcf,'type','axes'),'CLim',[-3 3]);

%% 3. creat data for hmz 
%  3.1 unhamonized data
n_per_group = 15;
n_groups = 5;
x1 = []; x2 = [];
y1 = []; y2 = [];
group_labels = [];

x_centers = [2, 4, 6, 8, 10];
for i = 1:n_groups
    
    x1_group = randn(n_per_group, 1) * 0.8 + x_centers(i);
    y1_group = x1_group * 0.5 + randn(n_per_group, 1) * 0.8 - 3; 

    x2_group = randn(n_per_group, 1) * 0.8 + x_centers(i);
    y2_group = x2_group * 0.5 + randn(n_per_group, 1) * 0.3 - 3; 
    
    x1 = [x1; x1_group];
    y1 = [y1; y1_group];

    x2 = [x2; x2_group];
    y2 = [y2; y2_group];

    group_labels = [group_labels; ones(n_per_group, 1) * i];
end

load_paper
% subplot(2,2,1)
subplot(10, 7, 3);
hold on;
for i = 1:n_groups
    idx = group_labels == i;
    scatter(x1(idx), y1(idx), 10, colors_cyto_nor(i,:), 'filled', 'MarkerFaceAlpha', 0.6, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.1);
end
p = polyfit(x1, y1, 1);
x_fit = linspace(min(x1), max(x1), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 0.5);
box off;
xlim([-1 11]);
ylim([-4 4]); 

%  3.1 harmonized data
% subplot(2,2,2)
subplot(10, 7, 4);
y_norm = (y2 - min(y2)) / (max(y2) - min(y2));
y_norm = 1 - y_norm;
color_idx = round(y_norm * 255) + 1;
color_idx(color_idx > 256) = 256;
color_idx(color_idx < 1) = 1;
hold on;
for i = 1:length(x2)
    scatter(x2(i), y2(i), 10, blue(color_idx(i),:), 'filled', 'MarkerFaceAlpha', 0.6, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.1);
end
plot(x_fit, y_fit, 'k-', 'LineWidth', 0.5);
box off;
xlim([-1 11]);
ylim([-4 4]); 

%% 4. wscore
%  4.1 scatter
% subplot(2,2,3): 
subplot(10, 7, [10,11]);
hold on;
n_controls = 100;
age = 20 + rand(n_controls, 1) * 60;
age_squared = age.^2;
brain_measure = 0.5 + age * 0.008 - (age - 50).^2 * 0.00008 + randn(n_controls, 1) * 0.08;
scatter(age, brain_measure, 10, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.6, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.1);
tbl = table(age, age_squared, brain_measure, 'VariableNames', {'age', 'age_squared', 'brain_measure'});
mdl = fitlm(tbl, 'brain_measure ~ age + age_squared');
age_fit = linspace(min(age), max(age), 200)';
age_fit_squared = age_fit.^2;
tbl_predict = table(age_fit, age_fit_squared, 'VariableNames', {'age', 'age_squared'});
brain_fit = predict(mdl, tbl_predict);
plot(age_fit, brain_fit, 'Color', [0.4 0.5 0.7], 'LineWidth', 0.5);
residuals = mdl.Residuals.Raw;
SD = std(residuals);
plot(age_fit, brain_fit + 2*SD, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
plot(age_fit, brain_fit - 2*SD, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);

patient_age = 70;
patient_age_squared = patient_age^2;
tbl_patient = table(patient_age, patient_age_squared, 'VariableNames', {'age', 'age_squared'});
patient_expected = predict(mdl, tbl_patient);
patient_observed = patient_expected - 0.35;
plot(patient_age, patient_observed, 'o', 'MarkerSize', 2, ...
     'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'k', 'LineWidth', 0.1);
plot([patient_age patient_age], [patient_observed patient_expected], ...
     'r--', 'LineWidth', 1);
xlim([15 90]);
ylim([0 1.6]);
box off;

%  4.2 distribution
% subplot(2,2,4):
subplot(10, 7, 12);
hold on;
residuals = randn(1000, 1);
[f, xi] = ksdensity(residuals, 'Function', 'pdf', 'NumPoints', 300);
fill(xi, f, [0.7 0.7 0.7], 'EdgeColor', 'k', 'FaceAlpha', 0.6, 'LineWidth', 0.1);
x_dist = linspace(-4, 4, 300);
y_dist = normpdf(x_dist, 0, 1);
plot(x_dist, y_dist, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.5);
patient_residual = -2.5;
y_lim = get(gca, 'YLim');
plot([patient_residual patient_residual], [0 y_lim(2)], ...
     'r-', 'LineWidth', 1);
xlim([-4 4]);
box off;
axis off;