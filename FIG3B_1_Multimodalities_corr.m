clear; clc;
%% 1. load data, surface, color
load_color

%% 2. Multimodalities correlation
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')
neofeats = {'CorThick', 'CorCurv', 'MidFLAIR', 'Myelin', 'MidPET'};
stat_val = zeros(6,4);
load_paper

for i = 1:length(neofeats)-1
    x_HCs_temp = mean(eval(['Voxel_HCs_neofeats_hmz_flip_asym_wPDB.' neofeats{i}]),1);
    y_HCs_temp = mean(eval(['Voxel_HCs_neofeats_hmz_flip_asym_wPDB.' neofeats{5}]),1);
    x_PTs_temp = mean(eval(['Voxel_PTs_neofeats_hmz_flip_asym_wPDB.' neofeats{i}]),1);
    y_PTs_temp = mean(eval(['Voxel_PTs_neofeats_hmz_flip_asym_wPDB.' neofeats{5}]),1);
    x_HCs_ips = x_HCs_temp(1,1:end/2)'; y_HCs_ips = y_HCs_temp(1,1:end/2)';
    x_PTs_ips = x_PTs_temp(1,1:end/2)'; y_PTs_ips = y_PTs_temp(1,1:end/2)';

    writematrix(y_HCs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/pet_val_hc_3B.txt');
    writematrix(y_PTs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/pet_val_pt_3B.txt');

    writematrix(x_HCs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/HCs_CorThick.txt');
    writematrix(x_PTs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/PTs_CorThick.txt');

    writematrix(x_HCs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/HCs_CorCurv.txt');
    writematrix(x_PTs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/PTs_CorCurv.txt');
        
    writematrix(x_HCs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/HCs_MidFLAIR.txt');
    writematrix(x_PTs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/PTs_MidFLAIR.txt');
        
    writematrix(x_HCs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/HCs_Myelin.txt');
    writematrix(x_PTs_ips, '/host/verges/tank/data/mo/demo/BrainSMASH/PTs_Myelin.txt');

    % --- Outlier Detection ---
    mask_HCs = ~isoutlier(x_HCs_ips, "gesd") & ~isoutlier(y_HCs_ips, "gesd");
    mask_PTs = ~isoutlier(x_PTs_ips, "gesd") & ~isoutlier(y_PTs_ips, "gesd");
    fprintf('Feature %s | HCs: %d/%d kept | PTs: %d/%d kept\n', ...
        neofeats{i}, sum(mask_HCs), numel(mask_HCs), ...
        sum(mask_PTs), numel(mask_PTs));
    x_HCs_clean = x_HCs_ips(mask_HCs); y_HCs_clean = y_HCs_ips(mask_HCs);
    x_PTs_clean = x_PTs_ips(mask_PTs); y_PTs_clean = y_PTs_ips(mask_PTs);

    % --- HCs scatter_kde with gray colormap ---
    ax_hcs = subplot(10, 7, i+10);
    scatter_kde(x_HCs_clean, y_HCs_clean, 'filled', 'MarkerSize', 10);
    colormap(ax_hcs, gray(16:end,:));
    colorbar(ax_hcs, 'off');
    xlim([-max(abs(x_HCs_clean)), max(abs(x_HCs_clean))])
    % ylim([min(y_HCs_clean), max(y_HCs_clean)])
    % xlim([-0.2 0.2]);
    ylim([-0.06 0.06]);
    
    % --- PTs scatter_kde with green colormap ---
    ax_pts = subplot(10, 7, i+3);
    scatter_kde(x_PTs_clean, y_PTs_clean, 'filled', 'MarkerSize', 10);
    colormap(ax_pts, cmap_skew(16:end,:));
    colorbar(ax_pts, 'off');
    xlim([-max(abs(x_PTs_clean)), max(abs(x_PTs_clean))])
    % ylim([min(y_PTs_clean), max(y_PTs_clean)])   
    % xlim([-3.5 3.5]);
    ylim([-5.1 2.5]);

    % --- Correlation stats ---
    [rho_HCs_ips, pval_HCs_ips] = corr(x_HCs_clean, y_HCs_clean);
    [rho_PTs_ips, pval_PTs_ips] = corr(x_PTs_clean, y_PTs_clean);
    stat_val(1,i) = rho_HCs_ips; stat_val(2,i) = pval_HCs_ips;
    stat_val(4,i) = rho_PTs_ips; stat_val(5,i) = pval_PTs_ips;

    % --- Spin Permutation Tests ---
    % ⚡ Speed tip: reduce n_perm; 1000 is sufficient for p<0.05,
    %    10000 only needed for very precise p-values
    % ax_spin_HCs = subplot(10, 7, i+24);
    % axes(ax_spin_HCs);
    % spin_permutation_test(x_HCs_clean, y_HCs_clean, ...
    %     'n_perm', 10000, 'color', 'gray');          % reduced from 2000
    % title(ax_spin_HCs, sprintf('Spin: %s vs MidPET (HCs)', neofeats{i}));
    % 
    % ax_spin_PTs = subplot(10, 7, i+31);
    % axes(ax_spin_PTs);
    % spin_permutation_test(x_PTs_clean, y_PTs_clean, ...
    %     'n_perm', 10000, 'color', 'blue');          % reduced from 10000
    % title(ax_spin_PTs, sprintf('Spin: %s vs MidPET (PTs)', neofeats{i}));
end

% %% 3. output
% fig = gcf;
% fig.Renderer = 'painters';
% exportgraphics(fig, '/host/verges/tank/data/mo/demo/FIG3B_withingroup_corr.pdf', 'ContentType', 'vector');

%% 4. Pyramid Chart for the pathological modulation


load_paper

subplot(10,7,[1 8])
hold on;
y = 1:4;
barh(y, stat_val(1,1:4), 'FaceColor',Gray_G2, 'EdgeColor','none');
barh(y, stat_val(4,1:4), 'FaceColor',green_skew, 'EdgeColor','none');
xline(0, 'k', 'LineWidth', 1);

ax = gca;
ax.YTick = y;
ax.YTickLabel = neofeats(1:4);
ax.YDir = 'reverse';

xlabel('xlabel'); ylabel('ylabel');
xlim([-0.5 0.5]); ylim([0.4 4+0.6]);