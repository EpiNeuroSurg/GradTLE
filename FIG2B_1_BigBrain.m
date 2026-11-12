clear;clc;

%% 1. load data, surface, color
load_color
load_surface
load_brh

%% 2. load BigBrain data
profileData = load('/host/verges/tank/data/mo/software/BBW_BigData/spaces/tpl-fs_LR/tpl-fs_LR_den-32k_desc-profiles.txt');
Hist_G1   = -1*zscore(gifti('/host/verges/tank/data/mo/software/BBW_BigData/spaces/tpl-fs_LR/tpl-fs_LR_hemi-L_den-32k_desc-Hist_G1.shape.gii').cdata);
Hist_G2   = -1*zscore(gifti('/host/verges/tank/data/mo/software/BBW_BigData/spaces/tpl-fs_LR/tpl-fs_LR_hemi-L_den-32k_desc-Hist_G2.shape.gii').cdata);
Hist_mean = -1*zscore(mean(profileData(:, 1:end/2), 1)');
Hist_skew = -1*zscore(skewness(profileData(:, 1:end/2), 0, 1)');

Hist_G1(economo_parc(1:end/2,1) == 0) = min(Hist_G1);
plot_hemispheres(Hist_G1, surf_lh, 'views', 'lm'); colormap(cmap_G1)

Hist_G2(economo_parc(1:end/2,1) == 0) = min(Hist_G2);
plot_hemispheres(Hist_G2, surf_lh, 'views', 'lm'); colormap(cmap_G2)

Hist_mean(economo_parc(1:end/2,1) == 0) = min(Hist_mean);
plot_hemispheres(Hist_mean, surf_lh, 'views', 'lm'); colormap(cmap_mean)
set(findall(gcf,'type','axes'),'CLim',[-1 1]);

Hist_skew(economo_parc(1:end/2,1) == 0) = min(Hist_skew);
plot_hemispheres(Hist_skew, surf_lh, 'views', 'lm'); colormap(cmap_skew)
set(findall(gcf,'type','axes'),'CLim',[-1 1]);

%% 3. compute mean ± 95% CI per contour band for each BigBrain variable
vars      = [Hist_G1, Hist_G2, Hist_mean, Hist_skew];
var_names = {'Hist\_G1', 'Hist\_G2', 'Hist\_mean', 'Hist\_skew'};
cmaps     = {cmap_G1,   cmap_G2,   cmap_mean,   cmap_skew};

levels = 0:contour_step_mm:max_dist_mm;
nL     = numel(levels);
mean_all = nan(nL, 4);
ci_all   = nan(nL, 4);
n_all    = zeros(nL, 4);

for v = 1:4
    val = vars(:,v);
    for k = 1:nL
        band = is_TL(1:end/2) & abs(dist_roi(1:end/2) - levels(k)) <= contour_tol_mm;
        band = band & ~isnan(val);
        n_all(k,v) = nnz(band);
        if n_all(k,v) > 1
            s = val(band);
            mean_all(k,v) = mean(s, 'omitnan');
            ci_all(k,v)   = 1.96 * std(s, 'omitnan') / sqrt(n_all(k,v));
        end
    end
end

%% 4. plot the gradient
load_paper

n_colors = 192;
norm_idx = @(n) max(1, min(n_colors, round((n - min(n)) / max(max(n)-min(n),1) * (n_colors-1) + 1)));

% Plot 1: Hist_G1 overlayed with Hist_G2
ax1 = subplot(10, 7, [1, 2, 3, 8, 9, 10]);
hold(ax1, 'on');

for v = [1, 2]  % G1 and G2
    cmap_v = cmaps{v};
    cmap_v = cmap_v(65:256, :);
    idx_v  = norm_idx(n_all(:,v));
    valid = ~isnan(mean_all(:,v)) & ~isnan(ci_all(:,v));
    x     = levels(valid)';
    m     = mean_all(valid, v);
    ci    = ci_all(valid, v);
    idx_s = idx_v(valid);
    
    % CI fill
    fill(ax1, [x; flipud(x)], [m+ci; flipud(m-ci)], ...
        cmap_v(1,:), 'EdgeColor','none', 'FaceAlpha', 0.2);
    
    % Color-varying line
    for i = 1:numel(x)-1
        plot(ax1, x(i:i+1), m(i:i+1), '-', ...
            'Color', cmap_v(idx_s(i),:), 'LineWidth', 2);
    end
end

box(ax1, 'off');
xlabel(ax1, 'Geodesic distance (mm)');
ylabel(ax1, 'Hist\_G1 & Hist\_G2');
xlim(ax1, [-5 130]);
ylim(ax1, [-2 2]);
legend(ax1, {'G1', 'G2'}, 'Location', 'best');

% Colorbar for plot 1 (using G1 colormap and n_all)
colormap(ax1, cmaps{1}(65:256, :));
cb1 = colorbar(ax1, 'Location', 'eastoutside');
clim(ax1, [min(n_all(:,1)), max(n_all(:,1))]);
cb1.Label.String = 'n (G1)';
cb1.Ticks = [min(n_all(:,1)), max(n_all(:,1))];
cb1.TickLabels = {num2str(min(n_all(:,1))), num2str(max(n_all(:,1)))};

% Get position of ax1 to anchor the new colorbar
pos1 = get(ax1, 'Position');

% Create a new invisible axes for G2 colorbar
ax1b = axes('Position', [pos1(1)+pos1(3)+0.07, pos1(2), 0.02, pos1(4)]);
colormap(ax1b, cmaps{2}(65:256, :));
cb1b = colorbar(ax1b, 'Location', 'eastoutside');
clim(ax1b, [min(n_all(:,2)), max(n_all(:,2))]);
cb1b.Label.String = 'n (G2)';
cb1b.Ticks = [min(n_all(:,2)), max(n_all(:,2))];
cb1b.TickLabels = {num2str(min(n_all(:,2))), num2str(max(n_all(:,2)))};
axis(ax1b, 'off');  % hide the dummy axes

% Plot 2: Hist_mean overlayed with Hist_skew
ax2 = subplot(10, 7, [29, 30, 31, 36, 37, 38]);
hold(ax2, 'on');

for v = [3, 4]  % mean and skew
    cmap_v = cmaps{v};
    cmap_v = cmap_v(65:256, :);
    idx_v  = norm_idx(n_all(:,v));
    valid = ~isnan(mean_all(:,v)) & ~isnan(ci_all(:,v));
    x     = levels(valid)';
    m     = mean_all(valid, v);
    ci    = ci_all(valid, v);
    idx_s = idx_v(valid);
    
    % CI fill
    fill(ax2, [x; flipud(x)], [m+ci; flipud(m-ci)], ...
        cmap_v(1,:), 'EdgeColor','none', 'FaceAlpha', 0.2);
    
    % Color-varying line
    for i = 1:numel(x)-1
        plot(ax2, x(i:i+1), m(i:i+1), '-', ...
            'Color', cmap_v(idx_s(i),:), 'LineWidth', 2);
    end
end

box(ax2, 'off');
xlabel(ax2, 'Geodesic distance (mm)');
ylabel(ax2, 'Hist\_mean & Hist\_skew');
xlim(ax2, [-5 130]);
ylim(ax2, [-2,1.5])
legend(ax2, {'mean', 'skew'}, 'Location', 'best');

% Colorbar for plot 2 (using mean colormap and n_all)
colormap(ax2, cmaps{3}(65:256, :));
cb2 = colorbar(ax2, 'Location', 'eastoutside');
clim(ax2, [min(n_all(:,3)), max(n_all(:,3))]);
cb2.Label.String = 'n (mean)';
cb2.Ticks = [min(n_all(:,3)), max(n_all(:,3))];
cb2.TickLabels = {num2str(min(n_all(:,3))), num2str(max(n_all(:,3)))};

pos2 = get(ax2, 'Position');

ax2b = axes('Position', [pos2(1)+pos2(3)+0.07, pos2(2), 0.02, pos2(4)]);
colormap(ax2b, cmaps{4}(65:256, :));
cb2b = colorbar(ax2b, 'Location', 'eastoutside');
clim(ax2b, [min(n_all(:,4)), max(n_all(:,4))]);
cb2b.Label.String = 'n (skew)';
cb2b.Ticks = [min(n_all(:,4)), max(n_all(:,4))];
cb2b.TickLabels = {num2str(min(n_all(:,4))), num2str(max(n_all(:,4)))};
axis(ax2b, 'off');

hold on;
%% 5. plot the kde intensity 
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')
pet_val_pt = mean(Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:32492), 1)';
% 
% subplot(10,7, [4 5 11 12]);
% scatter_kde(Hist_G1, pet_val_pt, 'filled', 'MarkerSize', 10);
% colormap(gca, cmap_G1);   % <-- apply to current axes, not to h1
% colorbar;
% 
% subplot(10,7, [6 7 13 14]);
% scatter_kde(Hist_G2, pet_val_pt, 'filled', 'MarkerSize', 10);
% colormap(gca, cmap_G2);   % <-- apply to current axes, not to h1
% colorbar;
% 
% subplot(10,7, [32 33 39 40]);
% scatter_kde(Hist_mean, pet_val_pt, 'filled', 'MarkerSize', 10);
% colormap(gca, cmap_mean);   % <-- apply to current axes, not to h1
% colorbar;
% 
% subplot(10,7, [34 35 41 42]);
% scatter_kde(Hist_skew, pet_val_pt, 'filled', 'MarkerSize', 10);
% colormap(gca, cmap_skew);   % <-- apply to current axes, not to h1
% colorbar;
% 
% %% 6. permutation
% hold on;
% subplot(10,7,[18 19 25 26]);
% spin_permutation_test(Hist_G1, pet_val_pt,'n_perm', 10000, 'color',cmap_G1);
% 
% subplot(10,7,[20 21 27 28]);
% spin_permutation_test(Hist_G2, pet_val_pt,'n_perm', 10000, 'color',cmap_G2);
% 
% subplot(10,7,[46 47 53 54]);
% spin_permutation_test(Hist_mean, pet_val_pt,'n_perm', 10000, 'color',cmap_G1);
% 
% subplot(10,7,[48 49 55 56]);
% spin_permutation_test(Hist_skew, pet_val_pt,'n_perm', 10000, 'color',cmap_G1);

% --- Helper: remove outliers from a pair of vectors ---
clean_pair = @(a, b) deal( ...
    a(~isoutlier(a) & ~isoutlier(b)), ...
    b(~isoutlier(a) & ~isoutlier(b)) );

subplot(10,7, [4 5 11 12]);
[x1, y1] = clean_pair(Hist_G1, pet_val_pt);
scatter_kde(x1, y1, 'filled', 'MarkerSize', 10);
colormap(gca, cmap_G1);
colorbar off;
xlim([-2.5 2.5]);

subplot(10,7, [6 7 13 14]);
[x2, y2] = clean_pair(Hist_G2, pet_val_pt);
scatter_kde(x2, y2, 'filled', 'MarkerSize', 10);
colormap(gca, cmap_G2);
colorbar off;
xlim([-2.5 2.5]);

subplot(10,7, [32 33 39 40]);
[x3, y3] = clean_pair(Hist_mean, pet_val_pt);
scatter_kde(x3, y3, 'filled', 'MarkerSize', 10);
colormap(gca, cmap_mean);
colorbar off;
xlim([-2.5 2.5]);

subplot(10,7, [34 35 41 42]);
[x4, y4] = clean_pair(Hist_skew, pet_val_pt);
scatter_kde(x4, y4, 'filled', 'MarkerSize', 10);
colormap(gca, cmap_skew);
colorbar off;
xlim([-2.5 2.5]);

%% 6. permutation
% writematrix(pet_val_pt, '/host/verges/tank/data/mo/demo/BrainSMASH/pet_val_pt.txt');
% writematrix(Hist_G1, '/host/verges/tank/data/mo/demo/BrainSMASH/Hist_G1.txt');
% writematrix(Hist_G2, '/host/verges/tank/data/mo/demo/BrainSMASH/Hist_G2.txt');
% writematrix(Hist_mean, '/host/verges/tank/data/mo/demo/BrainSMASH/Hist_mean.txt');
% writematrix(Hist_skew, '/host/verges/tank/data/mo/demo/BrainSMASH/Hist_skew.txt');


hold on;
subplot(10,7,[18 19 25 26]);
spin_permutation_test(x1, y1, 'n_perm', 10000, 'color', cmap_G1);

subplot(10,7,[20 21 27 28]);
spin_permutation_test(x2, y2, 'n_perm', 10000, 'color', cmap_G2);

subplot(10,7,[46 47 53 54]);
spin_permutation_test(x3, y3, 'n_perm', 10000, 'color', cmap_mean);

subplot(10,7,[48 49 55 56]);
spin_permutation_test(x4, y4, 'n_perm', 10000, 'color', cmap_skew);

% %% 7. output
% fig = gcf;
% fig.Renderer = 'painters';
% exportgraphics(fig, '/host/verges/tank/data/mo/demo/FIG2B_BBW.pdf', 'ContentType', 'vector');