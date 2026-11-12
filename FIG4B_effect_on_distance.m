clear; clc;

%% 1. Load
load_color
load_brh
load('/host/verges/tank/data/mo/demo/42Voxel_Dvalue.mat')

%% 2. Define variables
% D_hippoFeats: 5 x 32492 (rows: HippoThick, HippoCurv, HippoFLAIR, HippoMyelin, HippoPET)
% D_Madal:      1 x 32492
vars      = {zscore(D_hippoFeats(1,:)'), -zscore(D_hippoFeats(2,:)'), ...
             -zscore(D_hippoFeats(3,:)'), zscore(D_hippoFeats(4,:)'), ...
             zscore(D_hippoFeats(5,:)'), -zscore(D_Madal')};
var_names = {'HippoThick','HippoCurv','HippoFLAIR','HippoMyelin','HippoPET','Mahalanobis'};
cmaps     = {cmap_purple, cmap_blue, cmap_green, cmap_orange, cmap_yellow, cmap_red};

%% 3. Compute mean ± 95% CI per contour band
levels = 0:contour_step_mm:max_dist_mm;
nL     = numel(levels);
nV     = 4;  % use left hemi only (1:end/2 assumes 64984 total; adjust if needed)
half   = numel(dist_roi)/2;
mean_all = nan(nL, 6);
ci_all   = nan(nL, 6);
n_all    = zeros(nL, 6);
for v = 1:6
    val = vars{v};
    for k = 1:nL
        band = is_TL(1:half) & abs(dist_roi(1:half) - levels(k)) <= contour_tol_mm;
        band = band & ~isnan(val(1:half));
        n_all(k,v) = nnz(band);
        if n_all(k,v) > 1
            s = val(band);
            mean_all(k,v) = mean(s, 'omitnan');
            ci_all(k,v)   = 1.96 * std(s, 'omitnan') / sqrt(n_all(k,v));
        end
    end
end
load_paper
n_colors = 192;
norm_idx = @(n) max(1, min(n_colors, round((n - min(n)) / max(max(n)-min(n), 1) * (n_colors-1) + 1)));

% --- Plot v = 1:5 ---
ax1 = subplot(10, 7, [23, 24, 25, 26, 30, 31, 32, 33, 37, 38, 39, 40]);
hold(ax1, 'on');
for v = 1:5
    cmap_v = cmaps{v}(65:256, :);
    idx_v  = norm_idx(n_all(:,v));
    valid  = ~isnan(mean_all(:,v)) & ~isnan(ci_all(:,v));
    x      = levels(valid)';
    m      = mean_all(valid, v);
    ci     = ci_all(valid, v);
    idx_s  = idx_v(valid);
    clr    = cmap_v(end, :);  % solid color for legend/fill
    % CI fill
    fill(ax1, [x; flipud(x)], [m+ci; flipud(m-ci)], clr, ...
        'EdgeColor','none', 'FaceAlpha', 0.2);
    % Color-varying line
    for i = 1:numel(x)-1
        plot(ax1, x(i:i+1), m(i:i+1), '-', ...
            'Color', cmap_v(idx_s(i),:), 'LineWidth', 2);
    end
end
box(ax1, 'off');
xlabel(ax1, 'Geodesic distance (mm)');
ylabel(ax1, 'CohensD');
xlim(ax1, [-5 130]); ylim(ax1, [-4 6]);

% --- Plot v = 6 ---
ax2 = subplot(10, 7, [44, 45, 46, 47, 51, 52, 53, 54, 58, 59, 60, 61]);
hold(ax2, 'on');
v = 6;
cmap_v = cmaps{v}(65:256, :);
idx_v  = norm_idx(n_all(:,v));
valid  = ~isnan(mean_all(:,v)) & ~isnan(ci_all(:,v));
x      = levels(valid)';
m      = mean_all(valid, v);
ci     = ci_all(valid, v);
idx_s  = idx_v(valid);
clr    = cmap_v(end, :);  % solid color for legend/fill
% CI fill
fill(ax2, [x; flipud(x)], [m+ci; flipud(m-ci)], clr, ...
    'EdgeColor','none', 'FaceAlpha', 0.2);
% Color-varying line
for i = 1:numel(x)-1
    plot(ax2, x(i:i+1), m(i:i+1), '-', ...
        'Color', cmap_v(idx_s(i),:), 'LineWidth', 2);
end
box(ax2, 'off');
xlabel(ax2, 'Geodesic distance (mm)');
ylabel(ax2, 'CohensD');
xlim(ax2, [-5 130]);

figure
spin_permutation_test(x, m,'n_perm', 10000);



%% 4. Permutation tests — dist_roi vs each variable
var_names_plot = {'HippoThick','HippoCurv','HippoFLAIR','HippoMyelin','HippoPET','Mahalanobis'};
cmaps_perm     = {cmap_purple, cmap_blue, cmap_green, cmap_orange, cmap_yellow, cmap_red};

% subplot positions for 6 panels (adjust indices to fit your layout)
subpos = {[18 19 25 26], [20 21 27 28], [22 23 29 30], ...
          [32 33 39 40], [34 35 41 42], [36 37 43 44]};

half = numel(dist_roi)/2;

for v = 1:6
    val  = vars{v};
    x_v  = dist_roi(1:half);
    y_v  = val(1:half);

    % keep only valid (non-NaN) TL vertices
    mask = is_TL(1:half) & ~isnan(x_v) & ~isnan(y_v);
    x_v  = x_v(mask);
    y_v  = y_v(mask);

    subplot(10, 7, subpos{v});
    hold on;
    spin_permutation_test(x_v, y_v, 'n_perm', 10000, 'color', cmaps_perm{v});
    title(var_names_plot{v}, 'Interpreter', 'none');
end