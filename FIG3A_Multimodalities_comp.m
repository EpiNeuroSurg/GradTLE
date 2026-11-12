clear;clc;

%% 1. load data
load_color
load_surface
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')

%% 2. plot the surface
featsName = {'CorThick','CorCurv','MidFLAIR','Myelin','MidPET'};

%  2.1 for good distribution
% for i = 1:length(featsName)
mesialIdx = economo_parc(1:32492) == 0;

% set the 'i'
% controls
indiFeats_HCs = mean(Voxel_HCs_neofeats_hmz_flip_asym_wPDB.(featsName{i})(:, 1:32492))';
outlierIdx_HCs = isoutlier(indiFeats_HCs(~mesialIdx), "gesd");
cleanData = indiFeats_HCs(~mesialIdx);
clim_vals = [min(cleanData(~outlierIdx_HCs)) max(cleanData(~outlierIdx_HCs))];
indiFeats_HCs(mesialIdx) = 0.1;
plot_hemispheres(indiFeats_HCs, surf_lh, 'views', 'lm'); colormap(flipud(blue))
set(findall(gcf,'type','axes'), 'CLim', [-0.1 0.1]);
cb = findall(gcf, 'type', 'colorbar');
set(cb, 'Ticks', linspace(clim_vals(1), clim_vals(2), 2), ...
    'TickLabels', round(linspace(clim_vals(1), clim_vals(2), 2), 2));

% patients
indiFeats_PTs = mean(Voxel_PTs_neofeats_hmz_flip_asym_wPDB.(featsName{i})(:, 1:32492))';
outlierIdx_PTs = isoutlier(indiFeats_PTs(~mesialIdx), "gesd");
cleanData_PTs = indiFeats_PTs(~mesialIdx);
clim_vals_PTs = [min(cleanData_PTs(~outlierIdx_PTs)) max(cleanData_PTs(~outlierIdx_PTs))];
indiFeats_PTs(mesialIdx) = clim_vals_PTs(2);
plot_hemispheres(indiFeats_PTs, surf_lh, 'views', 'lm'); colormap(flipud(blue))
set(findall(gcf,'type','axes'), 'CLim', clim_vals_PTs);
cb = findall(gcf, 'type', 'colorbar');
set(cb, 'Ticks', linspace(clim_vals_PTs(1), clim_vals_PTs(2), 2), ...
    'TickLabels', round(linspace(clim_vals_PTs(1), clim_vals_PTs(2), 2), 2));
% end

%  2.2 for bad distribution
mesialIdx = economo_parc(1:32492) == 0;
faces = surf_lh.tri;

% --- controls ---
indiFeats_HCs = mean(Voxel_HCs_neofeats_hmz_flip_asym_wPDB.(featsName{i})(:, 1:32492))';
nonMesialData = indiFeats_HCs(~mesialIdx);
p1  = prctile(nonMesialData, 0.2);
p99 = prctile(nonMesialData, 99.8);
spikeIdx = false(32492, 1);
spikeIdx(~mesialIdx) = nonMesialData < p1 | nonMesialData > p99;

indiFeats_HCs_plot = indiFeats_HCs;
for iter = 1:20
    remaining = find(spikeIdx);
    if isempty(remaining); fprintf('HCs converged at iter %d\n', iter); break; end
    for v = remaining'
        neighborIdx = unique(faces(any(faces == v, 2), :));
        neighborIdx = neighborIdx(neighborIdx ~= v);
        validNeighbors = neighborIdx(~spikeIdx(neighborIdx) & ~mesialIdx(neighborIdx));
        if ~isempty(validNeighbors)
            indiFeats_HCs_plot(v) = mean(indiFeats_HCs_plot(validNeighbors));
            spikeIdx(v) = false;
        end
    end
end
nonMesialClean = indiFeats_HCs_plot(~mesialIdx);
clim_vals_HCs = [min(nonMesialClean) max(nonMesialClean)];
indiFeats_HCs_plot(mesialIdx) = clim_vals_HCs(2);

plot_hemispheres(indiFeats_HCs_plot, surf_lh, 'views', 'lm'); colormap(flipud(blue))
set(findall(gcf,'type','axes'), 'CLim', clim_vals_HCs);
cb = findall(gcf, 'type', 'colorbar');
set(cb, 'Ticks', [clim_vals_HCs(1) clim_vals_HCs(2)], ...
    'TickLabels', round([clim_vals_HCs(1) clim_vals_HCs(2)], 2));

% --- patients ---
indiFeats_PTs = mean(Voxel_PTs_neofeats_hmz_flip_asym_wPDB.(featsName{i})(:, 1:32492))';
nonMesialData = indiFeats_PTs(~mesialIdx);
p1  = prctile(nonMesialData, 0.2);
p99 = prctile(nonMesialData, 99.8);
spikeIdx = false(32492, 1);
spikeIdx(~mesialIdx) = nonMesialData < p1 | nonMesialData > p99;

indiFeats_PTs_plot = indiFeats_PTs;
for iter = 1:20
    remaining = find(spikeIdx);
    if isempty(remaining); fprintf('PTs converged at iter %d\n', iter); break; end
    for v = remaining'
        neighborIdx = unique(faces(any(faces == v, 2), :));
        neighborIdx = neighborIdx(neighborIdx ~= v);
        validNeighbors = neighborIdx(~spikeIdx(neighborIdx) & ~mesialIdx(neighborIdx));
        if ~isempty(validNeighbors)
            indiFeats_PTs_plot(v) = mean(indiFeats_PTs_plot(validNeighbors));
            spikeIdx(v) = false;
        end
    end
end
nonMesialClean = indiFeats_PTs_plot(~mesialIdx);
clim_vals_PTs = [min(nonMesialClean) max(nonMesialClean)];
indiFeats_PTs_plot(mesialIdx) = clim_vals_PTs(2);

plot_hemispheres(indiFeats_PTs_plot, surf_lh, 'views', 'lm'); colormap(flipud(blue))
set(findall(gcf,'type','axes'), 'CLim', clim_vals_PTs);
cb = findall(gcf, 'type', 'colorbar');
set(cb, 'Ticks', [clim_vals_PTs(1) clim_vals_PTs(2)], ...
    'TickLabels', round([clim_vals_PTs(1) clim_vals_PTs(2)], 2));

%% 3. within-group comparison
left_blue  = blue(200,:);    % features 1,2,5
left_red   = red(200,:);     % features 3,4
right_col  = gray(200,:);

left_alpha  = 1.0;           % left: opaque
right_alpha = 0.30;          % right: transparent

% --- params
inMaskIdx = ismember(economo_parc, allTempLbls);
load_paper
feat_red   = [3 4];          % these two use red on the left
plot_order = [1 2 3 4 5];    % ensure 3rd/4th are #3 and #4

cx=1; v_off=0.2; s_off=0.2; box_w=0.08; v_w=0.15; lw=0.5;

% stats container
stats = struct('ttest_p', cell(1,numel(plot_order)), ...
               'cohen_d', cell(1,numel(plot_order)), ...
               'n1',      cell(1,numel(plot_order)), ...
               'n2',      cell(1,numel(plot_order)));

for k = 1:numel(plot_order)
    i = plot_order(k);
    left_col = left_blue; if ismember(i, feat_red), left_col = left_red; end

    subplot(10,7,7*k); hold on;

    % -------- data (KEEPING YOUR INPUTS AS-IS) --------
    d_h = mean(Voxel_PTs_neofeats_hmz_flip_asym_wPDB.(featsName{i})(:,1:32492))'; % left series
    d_d = mean(Voxel_HCs_neofeats_hmz_flip_asym_wPDB.(featsName{i})(:,1:32492))'; % right series
    msk = inMaskIdx(1:end/2,1);
    d_h = d_h(msk~=0); d_d = d_d(msk~=0);
    d_h = d_h(~isoutlier(d_h,'gesd')); 
    d_d = d_d(~isoutlier(d_d,'gesd'));

    % -------- stats: two-sample t-test + Cohen's d --------
    group1 = d_h; group2 = d_d;                             % respect your left/right data
    [~, p] = ttest2(group1, group2);
    m1 = mean(group1); m2 = mean(group2);
    s1 = std(group1);  s2 = std(group2);
    n1 = numel(group1); n2 = numel(group2);
    pooledSD = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2));
    d = (m1 - m2) / pooledSD;

    stats(i).ttest_p = p; stats(i).cohen_d = d; stats(i).n1 = n1; stats(i).n2 = n2;

    % -------- left violin+box (opaque) --------
    vx = cx - v_off; [f,xi] = ksdensity(d_h,'Bandwidth',0.3); f = v_w*(f/max(f));
    fill([vx - f, vx*ones(1,numel(xi))],[xi, fliplr(xi)], left_col, ...
        'FaceAlpha', left_alpha, 'EdgeColor','k','LineWidth',lw,'EdgeAlpha',1);

    bx = cx - s_off; q = prctile(d_h,[25 50 75]); iqr = q(3)-q(1);
    wlo = min(d_h(d_h >= q(1)-0.5*iqr)); whi = max(d_h(d_h <= q(3)+0.5*iqr));
    fill([bx-box_w bx+box_w bx+box_w bx-box_w],[q(1) q(1) q(3) q(3)], left_col, ...
        'FaceAlpha', left_alpha, 'EdgeColor','k','LineWidth',lw,'EdgeAlpha',1);
    plot([bx-box_w bx+box_w],[q(2) q(2)],'-','Color','w','LineWidth',lw);
    plot([bx bx],[wlo q(1)],'k-','LineWidth',lw); plot([bx bx],[q(3) whi],'k-','LineWidth',lw);

    % -------- right violin+box (transparent gray) --------
    vx = cx + v_off; [f,xi] = ksdensity(d_d,'Bandwidth',0.3); f = v_w*(f/max(f));
    fill([vx + f, vx*ones(1,numel(xi))],[xi, fliplr(xi)], right_col, ...
        'FaceAlpha', right_alpha, 'EdgeColor','k','LineWidth',lw,'EdgeAlpha',1);

    bx = cx + s_off; q = prctile(d_d,[25 50 75]); iqr = q(3)-q(1);
    wlo = min(d_d(d_d >= q(1)-0.5*iqr)); whi = max(d_d(d_d <= q(3)+0.5*iqr));
    fill([bx-box_w bx+box_w bx+box_w bx-box_w],[q(1) q(1) q(3) q(3)], right_col, ...
        'FaceAlpha', right_alpha, 'EdgeColor','k','LineWidth',lw,'EdgeAlpha',1);
    plot([bx-box_w bx+box_w],[q(2) q(2)],'-','Color','w','LineWidth',lw);
    plot([bx bx],[wlo q(1)],'k-','LineWidth',lw); plot([bx bx],[q(3) whi],'k-','LineWidth',lw);

    % -------- axes & on-plot annotation --------
    set(gca,'XTick',[cx-v_off cx+v_off],'XTickLabel',{'HCs','PTs'},'FontSize',5,'Box','off');
    xlim([cx - v_off - v_w - 0.08, cx + v_off + v_w + 0.08]);

    % show p and d (top-left, normalized coordinates)
    text(0.02, 0.95, sprintf('p=%.3e, d=%.3f', p, d), ...
        'Units','normalized','FontSize',6,'VerticalAlignment','top');
end