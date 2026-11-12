clear;clc;

%% 1. load_data
load_brh

%% 2. draw the gradient
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')

% ---- Extract PET values ----
pet_val_pt = mean(Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:32492), 1)';
pet_val_hc = mean(Voxel_HCs_neofeats_hmz_flip_asym_wPDB.MidPET(:,1:32492), 1)';

econ_lab = economo_parc(1:end/2, 1);

if numel(econ_lab) ~= size(F_V,1) || numel(pet_val_pt) ~= size(F_V,1)
    error("Label/PET length mismatch with fsLR vertices. labels=%d PET=%d V=%d", ...
        numel(econ_lab), numel(pet_val_pt), size(F_V,1));
end

% =================== 1) Contours + mean±95%CI PET per contour ===================
if isempty(max_dist_mm)
    max_dist_mm = nanmax(dist_roi);
end
levels = 0:contour_step_mm:max_dist_mm;
nL = numel(levels);

% Preallocate
mean_pt = nan(nL,1); ci_pt = nan(nL,1);
mean_hc = nan(nL,1); ci_hc = nan(nL,1);
n_pt    = zeros(nL,1); n_hc = zeros(nL,1);

for k = 1:nL
    L = levels(k);
    band = is_TL(1:end/2,1) & abs(dist_roi(1:end/2,1) - L) <= contour_tol_mm;

    % --- Patients ---
    band_pt = band & ~isnan(pet_val_pt);
    n_pt(k) = nnz(band_pt);
    if n_pt(k) > 1
        vals = pet_val_pt(band_pt);
        mean_pt(k) = mean(vals, 'omitnan');
        se = std(vals, 'omitnan') / sqrt(n_pt(k));
        ci_pt(k) = 1.96 * se;
    end

    % --- HC Controls ---
    band_hc = band & ~isnan(pet_val_hc);
    n_hc(k) = nnz(band_hc);
    if n_hc(k) > 1
        vals = pet_val_hc(band_hc);
        mean_hc(k) = mean(vals, 'omitnan');
        se = std(vals, 'omitnan') / sqrt(n_hc(k));
        ci_hc(k) = 1.96 * se;
    end
end

% Build summary table
T = table(levels(:), mean_pt, ci_pt, n_pt, mean_hc, ci_hc, n_hc, ...
    'VariableNames', {'distance_mm','mean_PET_PT','CI95_PT','n_PT', ...
                      'mean_PET_HC','CI95_HC','n_HC'});
disp(T(1:min(15,height(T)),:));

% =================== C) Mean ± 95% CI PET vs distance curve ===================
load_paper
ax_main = subplot(10, 7, [44, 45, 46, 47, 51, 52, 53, 54, 58, 59, 60, 61]);
hold(ax_main, 'on');

% --- CI fill colors (fixed) ---
pt_fill = blue(64,:);
hc_fill = gray(64,:);

% --- Gradient colormaps ---
cmap_pt  = blue(65:256, :);
cmap_hc  = gray(65:256, :);
n_colors = size(cmap_pt, 1);

% Normalize n_PT / n_HC to [1, n_colors]
norm_idx = @(n) max(1, min(n_colors, ...
    round( (n - min(n)) / max(max(n)-min(n), 1) * (n_colors-1) + 1 )));
idx_pt = norm_idx(n_pt);
idx_hc = norm_idx(n_hc);

% --- Valid mask ---
valid_hc = ~isnan(mean_hc) & ~isnan(ci_hc);
valid_pt = ~isnan(mean_pt) & ~isnan(ci_pt);

% --- CI fills ---
x_hc = levels(valid_hc)';
x_pt = levels(valid_pt)';

fill(ax_main, [x_hc; flipud(x_hc)], ...
    [mean_hc(valid_hc)+ci_hc(valid_hc); flipud(mean_hc(valid_hc)-ci_hc(valid_hc))], ...
    hc_fill, 'EdgeColor','none','FaceAlpha',0.3);

fill(ax_main, [x_pt; flipud(x_pt)], ...
    [mean_pt(valid_pt)+ci_pt(valid_pt); flipud(mean_pt(valid_pt)-ci_pt(valid_pt))], ...
    pt_fill, 'EdgeColor','none','FaceAlpha',0.3);

% --- Color-varying line for HC ---
lev_hc   = levels(valid_hc);
m_hc     = mean_hc(valid_hc);
idx_hc_v = idx_hc(valid_hc);
for i = 1:numel(lev_hc)-1
    plot(ax_main, lev_hc(i:i+1), m_hc(i:i+1), '-', ...
        'Color', cmap_hc(idx_hc_v(i),:), 'LineWidth', 2);
end

% --- Color-varying line for PT ---
lev_pt   = levels(valid_pt);
m_pt     = mean_pt(valid_pt);
idx_pt_v = idx_pt(valid_pt);
for i = 1:numel(lev_pt)-1
    plot(ax_main, lev_pt(i:i+1), m_pt(i:i+1), '-', ...
        'Color', cmap_pt(idx_pt_v(i),:), 'LineWidth', 2);
end

% --- Finalize main axes ---
box(ax_main, 'off');
xlabel(ax_main, 'Geodesic distance to bridgehead');
ylabel(ax_main, 'Mean PET within contour band');
ylim(ax_main, [-4.5 0.5]);
xlim(ax_main, [-5 130]);

% --- Two colorbars via dummy axes (one colormap each) ---
pos = ax_main.Position;  % [left bottom width height]

% PT colorbar (blue) - first slot to the right
ax_cb_pt = axes('Position', [pos(1)+pos(3)+0.01, pos(2), 0.015, pos(4)]);
colormap(ax_cb_pt, cmap_pt);
cb_pt = colorbar(ax_cb_pt, 'Location', 'eastoutside');
ax_cb_pt.Visible = 'off';
cb_pt.Label.String = 'n_{PT}';
cb_pt.Ticks = [0 1];
cb_pt.TickLabels = {num2str(min(n_pt)), num2str(max(n_pt))};

% HC colorbar (gray) - second slot, further right
ax_cb_hc = axes('Position', [pos(1)+pos(3)+0.01, pos(2), 0.015, pos(4)]);
colormap(ax_cb_hc, cmap_hc);
cb_hc = colorbar(ax_cb_hc, 'Location', 'eastoutside');
ax_cb_hc.Visible = 'off';
cb_hc.Position(1) = cb_pt.Position(1) + cb_pt.Position(3) + 0.02;  % shift right of PT bar
cb_hc.Label.String = 'n_{HC}';
cb_hc.Ticks = [0 1];
cb_hc.TickLabels = {num2str(min(n_hc)), num2str(max(n_hc))};

%% 4. abnormal ranges
% dist_plot = dist_roi;            % Nx1, concatenated L+R
% dist_plot(~is_TL) = 0;
% 
% range5 = dist_plot;
% range34 = dist_plot;
% range85 = dist_plot;
% 
% range5(dist_plot < 3 | dist_plot > 7) = 0;
% range34(dist_plot < 32 | dist_plot > 36) = 0;
% range85(dist_plot < 10 | dist_plot > 12) = 0;

% obj = plot_hemispheres(range5(1:end/2,1), surf_lh, 'views', 'lm');
% colormap(red);
% % set(findall(gcf,'type','axes'), 'CLim', [0 10]);
% % 
% obj = plot_hemispheres(range34(1:end/2,1), surf_lh, 'views', 'lm');
% colormap(red);
% % set(findall(gcf,'type','axes'), 'CLim', [30 40]);
% 
% obj = plot_hemispheres(range85(1:end/2,1), inflated_surf_lh, 'views', 'lm');
% colormap(red);
% % set(findall(gcf,'type','axes'), 'CLim', [10 20]);

% %% 5. intensity plot
% % x = distances (e.g., dist_roi within TL)
% % y = PET values (pet_val within TL)
% % mask = is_TL & ~isnan(dist_roi) & ~isnan(pet_val);
% % 
% % x = dist_roi(mask);
% % y = pet_val(mask);
% 
% x = levels';
% y = mean_pt;
% NN = [T.n_PT];
% % 2D histogram grid
% nx = 80; ny = 80;
% xedges = linspace(min(x), max(x), nx);
% yedges = linspace(min(y), max(y), ny);
% 
% N = histcounts2(x, y, xedges, yedges);
% 
% % Light smoothing (optional)
% N = imgaussfilt(N, 1.0);
% 
% % Plot filled contours
% [Xc, Yc] = meshgrid( ...
%     (xedges(1:end-1)+xedges(2:end))/2, ...
%     (yedges(1:end-1)+yedges(2:end))/2 );
% 
% figure('Color','w');
% contourf(Xc, Yc, N', 8, 'LineColor','none');  % transpose for orientation
% xlabel('Geodesic distance (mm)');
% ylabel('PET');
% title('2D density: PET vs distance');
% colorbar;
% 

%% 3. permutation
x1 = levels; y1 = mean_hc; n1 = [T.n_HC];
x2 = levels; y2 = mean_pt; n2 = [T.n_PT];
spin_permutation_test(x1(1:123), y1(1:123), 'n_perm', 10000);
spin_permutation_test(x2(1:123), y2(1:123), 'n_perm', 10000);

% save('/host/verges/tank/data/mo/demo/sample', 'x1', 'y1', 'n1', 'x2', 'y2', 'n2')
