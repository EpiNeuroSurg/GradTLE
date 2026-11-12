clear;clc;

%% 1. load data, surface, color
load_demo
load_surface
load_color

%% 2. Calculation the layer proportion of each brain region
FJK = economo_parc_5layers(economo_parc == 12); tabulate(FJK)
IA = economo_parc_5layers(economo_parc == 17); tabulate(IA)
IB = economo_parc_5layers(economo_parc == 18); tabulate(IB)

HA = economo_parc_5layers(economo_parc == 14); tabulate(HA)
HB = economo_parc_5layers(economo_parc == 15); tabulate(HB)
HC = economo_parc_5layers(economo_parc == 16); tabulate(HC)

TA = economo_parc_5layers(economo_parc == 37); tabulate(TA)
TB = economo_parc_5layers(economo_parc == 38); tabulate(TB)
TC = economo_parc_5layers(economo_parc == 39); tabulate(TC)
TD = economo_parc_5layers(economo_parc == 40); tabulate(TD)
TE = economo_parc_5layers(economo_parc == 41); tabulate(TE)
TF = economo_parc_5layers(economo_parc == 42); tabulate(TF)
TG = economo_parc_5layers(economo_parc == 43); tabulate(TG)

Region_prop = [...    
    20.00,  20.00,  20.00,  20.00,  20.00; ...
        0,  84.09,  0,      15.91,  0; ...     
     1.82,  97.57,  0,      0.61,   0; ...        
        0,  95.14,  2.87,   0,      1.99; ... 
    100.0,  0,      0,      0,      0; ...    
    18.03,  21.31,  60.66,  0,      0; ...     
     3.47,  3.82,   92.71,  0,      0; ...        
        0,  7.84,   90.44,  0,      1.72; ...        
        0,  5.47,   83.59,  0,      10.94; ...        
        0,  0,      3.62,   0,      96.38; ...
        0,  1.56,   9.38,   0,      89.06; ...        
        0,  97.69,  2.31,   0,      0; ...        
        0,  90.46,  9.54,   0,      0; ...     
     0.67,  94.64,  4.69,   0,      0];

figure; clf; hold on
for r = 1:size(Region_prop,1)
    bh = barh(r, Region_prop(r,:), 'stacked', 'BarWidth', 0.8);
    for c = 1:numel(bh)
        bh(c).FaceColor = colors_cyto_nor(c,:);
        bh(c).EdgeColor = 'k';
    end
end
ylim([0.5 size(Region_prop,1)+0.5]); set(gca,'YDir','reverse');
xlim([-10 110]);
set(gcf,'Renderer','painters');

%% 3. Plot cyto surface
% Create a custom colormap that maps each label to its specific color
max_label = max(allTempLbls);
custom_cmap = repmat([1 1 1], max_label + 1, 1);  
custom_cmap(1, :) = [1 1 1];  % White for label 0 (medial wall)

% Map labels to colors (index = label + 1 because MATLAB indexing starts at 1)
% Label 14 -> colors_cyto(1,:) [purple]
custom_cmap(14+1, :) = colors_cyto_nor(1, :);
% Label 12, 17, 18, 41, 42, 43 -> colors_cyto(2,:) [blue]
custom_cmap(12+1, :) = colors_cyto_nor(2, :);
custom_cmap(17+1, :) = colors_cyto_nor(2, :);
custom_cmap(18+1, :) = colors_cyto_nor(2, :);
custom_cmap(41+1, :) = colors_cyto_nor(2, :);
custom_cmap(42+1, :) = colors_cyto_nor(2, :);
custom_cmap(43+1, :) = colors_cyto_nor(2, :);
% Label 15, 16, 37, 38 -> colors_cyto(3,:) [green]
custom_cmap(15+1, :) = colors_cyto_nor(3, :);
custom_cmap(16+1, :) = colors_cyto_nor(3, :);
custom_cmap(37+1, :) = colors_cyto_nor(3, :);
custom_cmap(38+1, :) = colors_cyto_nor(3, :);
% Label 39, 40 -> colors_cyto(5,:) [light yellow]
custom_cmap(39+1, :) = colors_cyto_nor(5, :);
custom_cmap(40+1, :) = colors_cyto_nor(5, :);

%bkg inflated for AI edit 
bkg_inflated_brain = ones(32492,1);
plot_hemispheres(bkg_inflated_brain, inflated_surf_lh, 'views', 'lm'); colormap(hot); 
set(findall(gcf,'type','axes'),'CLim',[1 2]);

obj = plot_hemispheres(economo_parc(1:end/2), {inflated_surf_lh}, 'views', 'lm');
colormap(custom_cmap);
set(findall(gcf, 'type', 'axes'), 'CLim', [0 max_label]);

%  3.1 calculate the boundary
% Extract LEFT hemisphere labels & surface
nvert = size(inflated_surf_lh.coord,2);
label_L = double(economo_parc(1:nvert));   
tri_L   = inflated_surf_lh.tri;
coord_L = inflated_surf_lh.coord';

% Build all mesh edges
edges = [tri_L(:,[1 2]);
         tri_L(:,[2 3]);
         tri_L(:,[3 1])];
lab1 = label_L(edges(:,1));
lab2 = label_L(edges(:,2));
roi1 = ismember(lab1, allTempLbls);
roi2 = ismember(lab2, allTempLbls);

is_internal = roi1 & roi2 & (lab1 ~= lab2);
is_outer = (roi1 & ~roi2) | (~roi1 & roi2);
is_boundary = is_internal | is_outer;
boundary_edges = unique(sort(edges(is_boundary,:),2), 'rows');

v1 = coord_L(boundary_edges(:,1),:);
v2 = coord_L(boundary_edges(:,2),:);

X = [v1(:,1) v2(:,1) nan(size(v1,1),1)]';
Y = [v1(:,2) v2(:,2) nan(size(v1,1),1)]';
Z = [v1(:,3) v2(:,3) nan(size(v1,1),1)]';
X = X(:); Y = Y(:); Z = Z(:);

% Draw on ALL surface axes (BrainStat makes >1 axis; exclude colorbar)
axs = findobj(gcf,'Type','axes');
axs = axs(~strcmp(get(axs,'Tag'),'Colorbar'));
for k = 1:numel(axs)
    hold(axs(k),'on');
    plot3(axs(k), X, Y, Z, 'k-', 'LineWidth', 0.1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% economo cyto comparison
clear; clc;
rng('default')

%% load data
load_demo
load_color
load_surface

%% data preparation
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')

condition_names = {'1','2','3', '4', '5'};
num_cond = length(condition_names);

cytoData_HCs = zeros(length(Demo_HCs.ID), num_cond);
cytoData_PTs = zeros(length(Demo_disPTs.ID), num_cond);

for i = 1:num_cond
    idxTmp = ismember(economo_parc, allTempLbls) & (economo_parc_5layers == i);
    cytoData_HCs(:,i) = mean(Voxel_HCs_neofeats_hmz_flip_asym_wPDB.MidPET(:, idxTmp),2);  
    cytoData_PTs(:,i) = mean(Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(:, idxTmp),2);  
end

%% Cohen's d and p-value (two-sample t-test, HC vs PT per condition)
cohens_d = zeros(1, num_cond);
p_values  = zeros(1, num_cond);
t_values  = zeros(1, num_cond);

for i = 1:num_cond
    hc = cytoData_HCs(:, i);
    pt = cytoData_PTs(:, i);

    % Two-sample t-test
    [~, p, ~, stats] = ttest2(hc, pt);
    t_values(i)  = stats.tstat;
    p_values(i)  = p;

    % Pooled Cohen's d
    n_hc = length(hc);  n_pt = length(pt);
    pooled_std = sqrt(((n_hc-1)*std(hc)^2 + (n_pt-1)*std(pt)^2) / (n_hc + n_pt - 2));
    cohens_d(i) = (mean(pt) - mean(hc)) / pooled_std;
end
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values, 0.05, 'pdep', 'yes');

%% parameters
data_humans = cytoData_HCs;
data_dogs   = cytoData_PTs;

cond_spacing       = 1.0;    % condition 中心间距
violin_offset      = 0.22;   % violin 中心距 condition 中心
scatter_offset     = 0.12;   % scatter 中心距 condition 中心
violin_width_scale = 0.20;   % violin 最大宽度
human_alpha        = 0.20;   % Humans 所有元素不透明度（70% 透明）
dog_alpha          = 1.00;   % Dogs violin 不透明度
lw                 = 0.1;    % 统一线条粗细
box_half_w         = 0.04;   % 箱图半宽
scatter_sz         = 5;     % 散点大小
scatter_alpha      = 0.75;

cond_x = (0:num_cond-1) * cond_spacing;   % condition 中心 x 坐标

%% plot
figure('Units', 'centimeters', 'Position', [2 2 21 29.7]); 
subplot(10,7,[1 2 3 8 9 10]);
% figure
hold on;

for c = 1:num_cond
    cx  = cond_x(c);
    col = colors_cyto_nor(c,:);

    d_h = data_humans(:,c);
    d_d = data_dogs(:,c);

    % first group
    vx_h = cx - violin_offset;

    % --- violin
    [f_h, xi_h] = ksdensity(d_h, 'Bandwidth', 0.3);
    f_h = violin_width_scale * (f_h / max(f_h));
    fill([vx_h - f_h, vx_h*ones(1,length(xi_h))], ...
         [xi_h, fliplr(xi_h)], col, ...
         'FaceAlpha', human_alpha, 'EdgeColor', 'k', 'LineWidth', lw);

    % --- scatter（Humans）---
    sx_h = cx - scatter_offset;
    jitter_h = (box_half_w * 1.5) .* (0.5 - rand(size(d_h)));
    % scatter(sx_h + jitter_h, d_h, scatter_sz, ...
    %     'MarkerFaceColor', col, ...
    %     'MarkerEdgeColor', 'k', ...
    %     'LineWidth', lw, ...
    %     'MarkerFaceAlpha', human_alpha);

    % --- 箱图    
    bx_h = sx_h;
    % bx_h = vx_h + box_half_w * 0.5;
    q_h  = prctile(d_h, [25 50 75]);
    iqr_h = q_h(3) - q_h(1);
    wlo_h = min(d_h(d_h >= q_h(1) - 0.5*iqr_h));
    whi_h = max(d_h(d_h <= q_h(3) + 0.5*iqr_h));
    fill([bx_h-box_half_w, bx_h+box_half_w, bx_h+box_half_w, bx_h-box_half_w], ...
         [q_h(1), q_h(1), q_h(3), q_h(3)], col, ...
         'FaceAlpha', human_alpha, 'EdgeColor', 'k', 'LineWidth', lw);
    plot([bx_h-box_half_w, bx_h+box_half_w], [q_h(2), q_h(2)], '-', ...
         'Color', 'w', 'LineWidth', lw);
    plot([bx_h bx_h], [wlo_h q_h(1)], 'k-', 'LineWidth', lw);
    plot([bx_h bx_h], [q_h(3) whi_h], 'k-', 'LineWidth', lw);

    % second group
    vx_d = cx + violin_offset;

    % --- violin
    [f_d, xi_d] = ksdensity(d_d, 'Bandwidth', 0.3);
    f_d = violin_width_scale * (f_d / max(f_d));
    fill([vx_d + f_d, vx_d*ones(1,length(xi_d))], ...
         [xi_d, fliplr(xi_d)], col, ...
         'FaceAlpha', dog_alpha, 'EdgeColor', 'k', 'LineWidth', lw);

    % --- scatter（Dogs）---
    sx_d = cx + scatter_offset;
    jitter_d = (box_half_w * 1.5) .* (0.5 - rand(size(d_d)));
    % scatter(sx_d + jitter_d, d_d, scatter_sz, ...
    %     'MarkerFaceColor', col, ...
    %     'MarkerEdgeColor', 'k', ...
    %     'LineWidth', lw, ...
    %     'MarkerFaceAlpha', scatter_alpha);

    % --- box plot
    bx_d = sx_d;
    % bx_d = vx_d - box_half_w * 0.5;
    q_d  = prctile(d_d, [25 50 75]);
    iqr_d = q_d(3) - q_d(1);
    wlo_d = min(d_d(d_d >= q_d(1) - 0.5*iqr_d));
    whi_d = max(d_d(d_d <= q_d(3) + 0.5*iqr_d));
    fill([bx_d-box_half_w, bx_d+box_half_w, bx_d+box_half_w, bx_d-box_half_w], ...
         [q_d(1), q_d(1), q_d(3), q_d(3)], col, ...
         'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', lw);
    plot([bx_d-box_half_w, bx_d+box_half_w], [q_d(2), q_d(2)], '-', ...
         'Color', 'w', 'LineWidth', lw);
    plot([bx_d bx_d], [wlo_d q_d(1)], 'k-', 'LineWidth', lw);
    plot([bx_d bx_d], [q_d(3) whi_d], 'k-', 'LineWidth', lw);

end

%% ----（linkline）-------------------------------------------------
med_h = median(data_humans);
med_d = median(data_dogs);

for c = 1:num_cond-1
    col_avg = (colors_cyto_nor(c,:) + colors_cyto_nor(c+1,:)) / 2;
    xh1 = cond_x(c)   - violin_offset + box_half_w * 0.5;
    xh2 = cond_x(c+1) - violin_offset + box_half_w * 0.5;
    plot([xh1 xh2], [med_h(c) med_h(c+1)], '--', 'Color', col_avg, 'LineWidth', lw);
    xd1 = cond_x(c)   + violin_offset - box_half_w * 0.5;
    xd2 = cond_x(c+1) + violin_offset - box_half_w * 0.5;
    plot([xd1 xd2], [med_d(c) med_d(c+1)], '--', 'Color', col_avg, 'LineWidth', lw);
end

%% -others
set(gca, 'XTick', cond_x, 'XTickLabel', condition_names, ...
    'FontSize', 5, 'Box', 'off');
ylabel('Performance');
xlim([cond_x(1) - violin_offset - violin_width_scale - 0.08, ...
      cond_x(end) + violin_offset + violin_width_scale + 0.08]);

hp1 = patch(NaN, NaN, [0.5 0.5 0.5], 'FaceAlpha', human_alpha, 'EdgeColor', 'none');
hp2 = patch(NaN, NaN, [0.5 0.5 0.5], 'FaceAlpha', dog_alpha,   'EdgeColor', 'none');
% legend([hp1 hp2], {'Humans (70% transparent)', 'Dogs'}, 'Location', 'northeast');