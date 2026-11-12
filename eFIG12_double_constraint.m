clear; clc; close all;

%% 1. load
load('/host/verges/tank/data/mo/demo/56surgMetric_SZProb.mat');

datasets = {
    Xraw,   Yraw,   'a. Raw - Feature 1',   1;
    Xraw,   Yraw,   'b. Raw - Feature 2',   2;
    Xclean, Yclean, 'c. Clean - Feature 1', 1;
    Xclean, Yclean, 'd. Clean - Feature 2', 2;
};

feat_names = {'Feature1', 'Feature2'};

%% baseline

baseline_WY = 0.69;   

%% 3. parameters
UPPER_Q  = 0.999;   
LOWER_Q  = 0.001;   
N_GRID   = 300;     
BW_FRAC  = 0.10;    
SMOOTH_P = 0.85;    
N_BOOT   = 200;     

%% 4. cycle
figure('Color','w','Position',[50 50 1600 480]);

for p = 1:4
    X_in  = datasets{p,1};
    Y_in  = datasets{p,2};
    label = datasets{p,3};
    col   = datasets{p,4};

    %% --- (A)  Logistic  ---
    X_table = array2table(X_in, 'VariableNames', feat_names);
    mdl = fitclinear(X_table, Y_in, 'Learner', 'logistic');

    % SF = 1）
    posClass = 1;
    [~, posIdx] = ismember(posClass, mdl.ClassNames);

    %% --- (B) Manual PDP ---
    x_vals   = linspace(min(X_in(:,col)), max(X_in(:,col)), 100)';
    pdp_vals = zeros(size(x_vals));

    for i = 1:numel(x_vals)
        X_tmp = X_table;
        X_tmp.(feat_names{col})(:) = x_vals(i);
        [~, score_tmp] = predict(mdl, X_tmp);
        pdp_vals(i) = mean(score_tmp(:, posIdx));
    end

    % （optional）
    % [pd_off, x_off] = partialDependence(mdl, feat_names{col}, posClass, X_table, ...
    %     'QueryPoints', x_vals);
    % fprintf('Max PDP diff (panel %d): %.2e\n', p, max(abs(pdp_vals - pd_off)));

    %% --- (C) mapping PDP value bach to each sample ---
    TEM = X_in(:, col);
    WY  = interp1(x_vals, pdp_vals, TEM, 'linear', 'extrap');

    %% --- (D) clice window---
    x_fine = linspace(min(TEM), max(TEM), N_GRID)';
    bw     = (max(TEM) - min(TEM)) * BW_FRAC;
    q_up   = zeros(size(x_fine));
    q_lo   = zeros(size(x_fine));

    for k = 1:numel(x_fine)
        idx = abs(TEM - x_fine(k)) < bw;
        if sum(idx) < 10
            [~, ord] = sort(abs(TEM - x_fine(k)));
            idx = false(numel(TEM),1);
            idx(ord(1:min(50, numel(ord)))) = true;
        end
        q_up(k) = quantile(WY(idx), UPPER_Q);   % 99.9%
        q_lo(k) = quantile(WY(idx), LOWER_Q);   % 0.1%  ← 修正
    end

    %% --- (E) RCS ---
    sp_up = csaps(x_fine, q_up, SMOOTH_P);
    sp_lo = csaps(x_fine, q_lo, SMOOTH_P);
    y_up  = fnval(sp_up, x_fine);
    y_lo  = fnval(sp_lo, x_fine);

    %% --- (F) Bootstrap  ---
    n       = numel(TEM);
    boot_up = zeros(numel(x_fine), N_BOOT);
    boot_lo = zeros(numel(x_fine), N_BOOT);

    for b = 1:N_BOOT
        bidx = randsample(n, n, true);
        Tb = TEM(bidx);
        Wb = WY(bidx);

        qu_b = zeros(size(x_fine));
        ql_b = zeros(size(x_fine));

        for k = 1:numel(x_fine)
            idx2 = abs(Tb - x_fine(k)) < bw;
            if sum(idx2) < 5
                [~, ord2] = sort(abs(Tb - x_fine(k)));
                idx2 = false(numel(Tb),1);
                idx2(ord2(1:min(30, numel(ord2)))) = true;
            end
            qu_b(k) = quantile(Wb(idx2), UPPER_Q);
            ql_b(k) = quantile(Wb(idx2), LOWER_Q);   
        end

        boot_up(:,b) = fnval(csaps(x_fine, qu_b, SMOOTH_P), x_fine);
        boot_lo(:,b) = fnval(csaps(x_fine, ql_b, SMOOTH_P), x_fine);
    end

    % upper 95% CI
    ci_up_hi = prctile(boot_up, 97.5, 2);
    ci_up_lo = prctile(boot_up,  2.5, 2);
    % lower 95% CI  
    ci_lo_hi = prctile(boot_lo, 97.5, 2);
    ci_lo_lo = prctile(boot_lo,  2.5, 2);

    dy        = gradient(y_up, x_fine);
    d2y       = gradient(dy,   x_fine);
    curvature = abs(d2y) ./ (1 + dy.^2).^1.5;

    % find peak
    [WY_peak, idx_peak] = max(y_up);
    TEM_peak            = x_fine(idx_peak);

    if idx_peak > 1
        curv_left   = curvature(1:idx_peak);
        [~, idx_th] = max(curv_left);
    else
        idx_th = 1;
    end
    TEM_threshold = x_fine(idx_th);
    WY_threshold  = y_up(idx_th);

    %% --- (H)  ---
    [WY_lo_min, idx_lo_min] = min(y_lo);
    TEM_lo_min = x_fine(idx_lo_min);

    %% --- (I) ---
    below_baseline = y_lo < baseline_WY;

    %% --- (J) plot ---
    subplot(1, 4, p); hold on;

    %-- 底层散点（全部数据，灰色）
    scatter(TEM, WY, 5, [0.82 0.82 0.82], 'filled', 'MarkerFaceAlpha', 0.4);

    %-- 上约束以上的点（红色高亮）
    mask_up = WY > interp1(x_fine, y_up * 0.93, TEM, 'linear', 'extrap');
    scatter(TEM(mask_up), WY(mask_up), 18, [0.85 0.18 0.18], ...
        'filled', 'MarkerFaceAlpha', 0.8);

    %-- 下约束以下的点（蓝色高亮）
    mask_lo = WY < interp1(x_fine, y_lo * 1.07, TEM, 'linear', 'extrap');
    scatter(TEM(mask_lo), WY(mask_lo), 18, [0.18 0.40 0.75], ...
        'filled', 'MarkerFaceAlpha', 0.8);

    %-- 上约束 95% CI 阴影（红色，两侧对称）
    fill([x_fine; flipud(x_fine)], [ci_up_hi; flipud(ci_up_lo)], ...
         [0.85 0.18 0.18], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    %-- 下约束 95% CI 阴影（蓝色，两侧对称）← 修正：之前视觉上不对称
    fill([x_fine; flipud(x_fine)], [ci_lo_hi; flipud(ci_lo_lo)], ...
         [0.18 0.40 0.75], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    %-- 生态风险区（黄色背景）← 新增
    if any(below_baseline)
        x_risk = x_fine(below_baseline);
        y_ax   = ylim;
        fill([x_risk(1) x_risk(end) x_risk(end) x_risk(1)], ...
             [y_ax(1) y_ax(1) y_ax(2) y_ax(2)], ...
             [1 0.85 0.2], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    end

    %-- 上约束曲线（红）& 下约束曲线（蓝）
    plot(x_fine, y_up, 'r-', 'LineWidth', 2.5);
    plot(x_fine, y_lo, 'b-', 'LineWidth', 2.5);

    %-- 基线（绿色虚线）← 新增
    yline(baseline_WY, '--', 'Color', [0.10 0.60 0.10], ...
          'LineWidth', 1.4, 'Label', sprintf('Baseline=%.2f', baseline_WY), ...
          'LabelHorizontalAlignment', 'right');

    %-- 中位线
    yline(median(WY), 'k-', 'LineWidth', 1.2);

    %-- 上约束阈值（曲率变化点）← 修正：原代码用最大值点
    plot([TEM_threshold TEM_threshold], [min(WY) WY_threshold], ...
         'k--', 'LineWidth', 1.0);
    plot(TEM_threshold, WY_threshold, 'kv', ...
         'MarkerSize', 7, 'MarkerFaceColor', 'k');
    text(TEM_threshold, WY_threshold + 0.02, ...
         sprintf('Curvature\nThreshold\n(%.2f)', TEM_threshold), ...
         'FontSize', 7, 'HorizontalAlignment', 'center');

    %-- 下约束最低点
    plot(TEM_lo_min, WY_lo_min, 'b^', ...
         'MarkerSize', 7, 'MarkerFaceColor', 'b');

    xlabel(feat_names{col}, 'FontSize', 11);
    ylabel('SF Probability',  'FontSize', 11);
    title(label, 'FontSize', 11, 'FontWeight', 'bold');
    ylim([0.68 0.88]);
    box on; grid on;

    %-- 记录第4个子图的数据，供单独出图使用
    if p == 4
        panel4.TEM          = TEM;
        panel4.WY           = WY;
        panel4.x_fine       = x_fine;
        panel4.y_up         = y_up;
        panel4.y_lo         = y_lo;
        panel4.ci_up_hi     = ci_up_hi;
        panel4.ci_up_lo     = ci_up_lo;
        panel4.ci_lo_hi     = ci_lo_hi;
        panel4.ci_lo_lo     = ci_lo_lo;
        panel4.mask_up      = mask_up;
        panel4.mask_lo      = mask_lo;
        panel4.below_baseline = below_baseline;
        panel4.TEM_threshold  = TEM_threshold;
        panel4.WY_threshold   = WY_threshold;
        panel4.TEM_lo_min     = TEM_lo_min;
        panel4.WY_lo_min      = WY_lo_min;
        panel4.feat_col       = col;
        panel4.label          = label;
    end
end

%% 5. 图例（四合一图）
lg = legend({'All data', 'Upper pts', 'Lower pts', ...
             'Upper 95% CI', 'Lower 95% CI', ...
             'Ecological risk zone', ...
             'Upper constraint', 'Lower constraint', ...
             'Baseline', 'Median'}, ...
    'Orientation', 'horizontal', 'FontSize', 8, ...
    'Position', [0.18 0.01 0.64 0.04]);

%% 6. 单独绘制第4个子图（subplot布局，cbrewer配色，修正散点着色）
load_color
blue = max(0, min(1, blue));   
red  = max(0, min(1, red));

c_red  = red(200,:);
c_blue = blue(200,:);
c_gray = gray(56,:);

load_paper
ax = subplot(10, 7, [1 2 3 8 9 10 15 16 17 22 23 24]);
hold(ax, 'on');

hold on

d = panel4;

% ===== 下约束也用曲率点（替换最低点逻辑）=====
dy_lo        = gradient(d.y_lo, d.x_fine);
d2y_lo       = gradient(dy_lo, d.x_fine);
curvature_lo = abs(d2y_lo) ./ (1 + dy_lo.^2).^1.5;
[~, idx_lo_th] = max(curvature_lo);
d.TEM_lo_min = d.x_fine(idx_lo_th);
d.WY_lo_min  = d.y_lo(idx_lo_th);

pct_up  = prctile(d.WY, 95);   % 前5%高值 → 红
pct_lo  = prctile(d.WY, 5);    % 后5%低值 → 蓝
mask_above = d.WY >= pct_up;
mask_below = d.WY <= pct_lo;
mask_mid   = ~mask_above & ~mask_below;

% --- 1. 灰色散点 ---
scatter(ax, d.TEM(mask_mid), d.WY(mask_mid), 20, c_gray, ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7);

% --- 2. 红色散点（高值，上约束区）---
scatter(ax, d.TEM(mask_above), d.WY(mask_above), 20, c_red, ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.85);

% --- 3. 蓝色散点（低值，下约束区）---
scatter(ax, d.TEM(mask_below), d.WY(mask_below), 20, c_blue, ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.85);

% --- 4. 上约束 95% CI 阴影 ---
fill(ax, [d.x_fine; flipud(d.x_fine)], [d.ci_up_hi; flipud(d.ci_up_lo)], ...
     c_red, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% --- 5. 下约束 95% CI 阴影 ---
fill(ax, [d.x_fine; flipud(d.x_fine)], [d.ci_lo_hi; flipud(d.ci_lo_lo)], ...
     c_blue, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% --- 6. 生态风险区（黄色）---
if any(d.below_baseline)
    x_risk = d.x_fine(d.below_baseline);
    fill(ax, [x_risk(1) x_risk(end) x_risk(end) x_risk(1)], ...
         [0.68 0.68 0.88 0.88], ...
         [1 0.85 0.2], 'FaceAlpha', 0.13, 'EdgeColor', 'none');
end

% --- 7. 约束曲线 ---
plot(ax, d.x_fine, d.y_up, '-', 'Color', c_red,  'LineWidth', 2.5);
plot(ax, d.x_fine, d.y_lo, '-', 'Color', c_blue, 'LineWidth', 2.5);

% --- 8. 基线（黑色粗实线）---
yline(ax, baseline_WY, 'k-', 'LineWidth', 2.0);
text(ax, min(d.TEM) + 0.02*(max(d.TEM)-min(d.TEM)), baseline_WY + 0.003, ...
     sprintf('Baseline=%.3f', baseline_WY), ...
     'FontSize', 7.5, 'Color', 'k', 'FontWeight', 'bold');

% --- 9. 中位线 ---
yline(ax, median(d.WY), 'k--', 'LineWidth', 1.0);

% --- 10. 上约束阈值（曲率变化点）---
plot(ax, [d.TEM_threshold d.TEM_threshold], [0.68 d.WY_threshold], ...
     'k--', 'LineWidth', 1.2);
plot(ax, d.TEM_threshold, d.WY_threshold, 'o', ...
     'MarkerSize', 7, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.4);
text(ax, d.TEM_threshold, d.WY_threshold + 0.005, ...
     sprintf('(%.2f, %.4f)', d.TEM_threshold, d.WY_threshold), ...
     'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.2]);

% --- 11. 下约束最低点 ---
plot(ax, [d.TEM_lo_min d.TEM_lo_min], [0.68 d.WY_lo_min], ...
     'k--', 'LineWidth', 1.2);
plot(ax, d.TEM_lo_min, d.WY_lo_min, 'o', ...
     'MarkerSize', 7, 'MarkerEdgeColor', c_blue, 'MarkerFaceColor', 'w', 'LineWidth', 1.4);
text(ax, d.TEM_lo_min, d.WY_lo_min - 0.007, ...
     sprintf('%.4f', d.TEM_lo_min), ...
     'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', c_blue);

% --- 轴标签 & 样式 ---
xlabel(ax, feat_names{d.feat_col}, 'FontSize', 10);
ylabel(ax, 'SF Probability', 'FontSize', 10);
title(ax, d.label, 'FontSize', 11, 'FontWeight', 'bold');
ylim(ax, [0.68 0.88]);
xlim(ax, [0.1 1.1]);
set(ax, 'FontSize', 9, 'LineWidth', 0.8);
box(ax, 'on'); grid(ax, 'on');

legend(ax, {'Mid data', 'Upper pts', 'Lower pts', ...
             'Upper 95% CI', 'Lower 95% CI', ...
             'Ecological risk zone', ...
             'Upper constraint', 'Lower constraint', ...
             'Baseline', 'Median'}, ...
    'Location', 'southeast', 'FontSize', 7, 'Box', 'on');