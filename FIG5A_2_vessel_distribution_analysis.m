clear;
clc;
close all;

%% ============================================================
% 1. Paths and colors
% =============================================================

DATA_FILE = ...
    '/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat';

HBT_FILE = [ ...
    '/host/verges/tank/data/mo/atlas/hippocampus/HBT/' ...
    'MNI152NLin2009cAsym/' ...
    'HippoHBT_space-HippUnfold_den-2mm_hemi-L.label.gii'];

ARTERIAL_DIR = ...
    '/host/verges/tank/data/mo/atlas/arterial/fsLR32k';

COLOR_FILE = '/host/verges/tank/data/mo/load_color.m';

%% Left arterial masks used for analysis

arterial_left_files = {
    fullfile(ARTERIAL_DIR,'ACTPL_hemi-L_fsLR32k.func.gii')
    fullfile(ARTERIAL_DIR,'MCATL_hemi-L_fsLR32k.func.gii')
    fullfile(ARTERIAL_DIR,'PCATL_hemi-L_fsLR32k.func.gii')
};

HBT_names = {'Head','Body','Tail'};
vessel_names = {'ACTP','MCAT','PCAT'};

files_to_check = [
    {
        DATA_FILE
        HBT_FILE
        COLOR_FILE
    }
    arterial_left_files(:)
];

for i = 1:numel(files_to_check)
    if ~isfile(files_to_check{i})
        error('File not found:\n%s',files_to_check{i});
    end
end

run(COLOR_FILE);

%% Colors

gray_color = [0.82 0.82 0.82];
highlight_red = red(end,:);

correlation_colormap = [
    flipud(blue)
    red
];

negative_color = blue(end,:);
positive_color = red(end,:);

%% ============================================================
% 2. Load PET features
% =============================================================

S = load(DATA_FILE);

MidPET = ...
    S.Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET;

HippoPET = ...
    S.Voxel_PTs_subfeats_hmz_flip_asym_wPDB.HippoPET;

fprintf('MidPET size:   %s\n',mat2str(size(MidPET)));
fprintf('HippoPET size: %s\n',mat2str(size(HippoPET)));

if size(MidPET,1) ~= size(HippoPET,1)
    error('MidPET and HippoPET have different patient numbers.');
end

N = size(MidPET,1);

%% ============================================================
% 3. Keep flipped ipsilateral/left PET features
% =============================================================

if size(MidPET,2) == 64984
    MidPET_ipsi = MidPET(:,1:32492);
elseif size(MidPET,2) == 32492
    MidPET_ipsi = MidPET;
else
    error('Unexpected MidPET width: %d.',size(MidPET,2));
end

if size(HippoPET,2) == 838
    HippoPET_ipsi = HippoPET(:,1:419);
elseif size(HippoPET,2) == 419
    HippoPET_ipsi = HippoPET;
else
    error('Unexpected HippoPET width: %d.',size(HippoPET,2));
end

%% ============================================================
% 4. HBT parcellation
% =============================================================

HBT_gii = gifti(HBT_FILE);
HBT_labels = double(HBT_gii.cdata(:));

if numel(HBT_labels) ~= 419
    error('HBT atlas must contain exactly 419 vertices.');
end

HBT_indices = {
    HBT_labels == 1
    HBT_labels == 2
    HBT_labels == 3
};

fprintf('\nHippocampal atlas:\n');

for i = 1:3
    fprintf('%s: %d vertices\n', ...
        HBT_names{i},nnz(HBT_indices{i}));
end

HBT_PET = nan(N,3);

for i = 1:3
    HBT_PET(:,i) = mean( ...
        HippoPET_ipsi(:,HBT_indices{i}), ...
        2,'omitnan');
end

HeadPET = HBT_PET(:,1);
BodyPET = HBT_PET(:,2);
TailPET = HBT_PET(:,3);

%% ============================================================
% 5. Load ipsilateral arterial masks
% =============================================================

arterial_left_masks = zeros(32492,3);

for i = 1:3

    left_gii = gifti(arterial_left_files{i});

    arterial_left_masks(:,i) = ...
        double(left_gii.cdata(:));

end

if size(arterial_left_masks,1) ~= 32492
    error('All arterial masks must contain 32,492 vertices.');
end

arterial_left_indices = ...
    isfinite(arterial_left_masks) & ...
    arterial_left_masks > 0;

fprintf('\nArterial territories:\n');

for i = 1:3
    fprintf( ...
        '%s: %d ipsilateral vertices\n', ...
        vessel_names{i}, ...
        nnz(arterial_left_indices(:,i)));
end

%% Mean ipsilateral PET in each arterial territory

VesselPET = nan(N,3);

for i = 1:3
    VesselPET(:,i) = mean( ...
        MidPET_ipsi(:,arterial_left_indices(:,i)), ...
        2,'omitnan');
end

ACTPPET = VesselPET(:,1);
MCATPET = VesselPET(:,2);
PCATPET = VesselPET(:,3);

%% ============================================================
% 6. HBT-to-territory correlations
% =============================================================

[rho_HBT_vessel,p_HBT_vessel] = corr( ...
    HBT_PET,VesselPET, ...
    'Type','Spearman', ...
    'Rows','pairwise');

%% BH-FDR across nine HBT-to-territory tests

p_vector = p_HBT_vessel(:);
valid_index = isfinite(p_vector);
valid_p = p_vector(valid_index);

[sorted_p,sort_order] = sort(valid_p);

number_tests = numel(sorted_p);

sorted_q = ...
    sorted_p .* number_tests ./ (1:number_tests)';

sorted_q = flipud(cummin(flipud(sorted_q)));
sorted_q(sorted_q > 1) = 1;

valid_q = nan(size(valid_p));
valid_q(sort_order) = sorted_q;

q_vector = nan(size(p_vector));
q_vector(valid_index) = valid_q;

q_HBT_vessel = reshape( ...
    q_vector,size(p_HBT_vessel));

fprintf('\nHBT-to-territory Spearman rho:\n');

disp(array2table( ...
    rho_HBT_vessel, ...
    'RowNames',HBT_names, ...
    'VariableNames',vessel_names));

fprintf('\nHBT-to-territory FDR q:\n');

disp(array2table( ...
    q_HBT_vessel, ...
    'RowNames',HBT_names, ...
    'VariableNames',vessel_names));

%% ============================================================
% 7. Correlations among arterial territories
% =============================================================

[rho_vessel,p_vessel] = corr( ...
    VesselPET, ...
    'Type','Spearman', ...
    'Rows','pairwise');

%% BH-FDR across three unique arterial correlations

upper_test_index = triu(true(3),1);
p_vector = p_vessel(upper_test_index);

valid_index = isfinite(p_vector);
valid_p = p_vector(valid_index);

[sorted_p,sort_order] = sort(valid_p);

number_tests = numel(sorted_p);

sorted_q = ...
    sorted_p .* number_tests ./ (1:number_tests)';

sorted_q = flipud(cummin(flipud(sorted_q)));
sorted_q(sorted_q > 1) = 1;

valid_q = nan(size(valid_p));
valid_q(sort_order) = sorted_q;

q_vector = nan(size(p_vector));
q_vector(valid_index) = valid_q;

q_vessel = zeros(3);
q_vessel(upper_test_index) = q_vector;
q_vessel = q_vessel + q_vessel';

fprintf('\nInter-territory Spearman correlations (3 decimals):\n');
fprintf('%10s %10s %10s %10s\n','',vessel_names{:});
for row = 1:3
    fprintf('%10s',vessel_names{row});
    fprintf(' %10.3f',rho_vessel(row,:));
    fprintf('\n');
end

%% ============================================================
% 8. Within-patient comparisons across HBT and arterial regions
% =============================================================

pair_index = [1 2; 1 3; 2 3];

complete_HBT = all(isfinite(HBT_PET),2);
HBT_complete = HBT_PET(complete_HBT,:);
[p_HBT_omnibus,HBT_friedman_table] = friedman(HBT_complete,1,'off');
HBT_chi_square = HBT_friedman_table{2,5};
HBT_kendall_W = HBT_chi_square / ...
    (size(HBT_complete,1)*(size(HBT_complete,2)-1));

HBT_pair_names = strings(3,1);
HBT_pair_N = zeros(3,1);
HBT_median_difference = nan(3,1);
HBT_pair_p = nan(3,1);

for pair = 1:3
    first_region = pair_index(pair,1);
    second_region = pair_index(pair,2);
    valid_pair = isfinite(HBT_PET(:,first_region)) & ...
                 isfinite(HBT_PET(:,second_region));
    first_values = HBT_PET(valid_pair,first_region);
    second_values = HBT_PET(valid_pair,second_region);
    HBT_pair_names(pair) = ...
        string(HBT_names{first_region}) + " vs " + ...
        string(HBT_names{second_region});
    HBT_pair_N(pair) = nnz(valid_pair);
    HBT_median_difference(pair) = median(first_values-second_values);
    HBT_pair_p(pair) = signrank(first_values,second_values);
end

[sorted_p,sort_order] = sort(HBT_pair_p);
sorted_q = sorted_p .* 3 ./ (1:3)';
sorted_q = flipud(cummin(flipud(sorted_q)));
sorted_q(sorted_q > 1) = 1;
HBT_pair_q = nan(3,1);
HBT_pair_q(sort_order) = sorted_q;

fprintf('\nHBT Friedman test: N=%d, chi-square=%.3f, p=%.3f, Kendall W=%.3f\n', ...
    size(HBT_complete,1),HBT_chi_square,p_HBT_omnibus,HBT_kendall_W);
HBT_stats_table = table( ...
    HBT_pair_names,HBT_pair_N,HBT_median_difference, ...
    HBT_pair_p,HBT_pair_q, ...
    'VariableNames',{'Comparison','N','MedianDifference','P','FDR_Q'});
disp(HBT_stats_table);

complete_vessel = all(isfinite(VesselPET),2);
Vessel_complete = VesselPET(complete_vessel,:);
[p_vessel_omnibus,vessel_friedman_table] = ...
    friedman(Vessel_complete,1,'off');
vessel_chi_square = vessel_friedman_table{2,5};
vessel_kendall_W = vessel_chi_square / ...
    (size(Vessel_complete,1)*(size(Vessel_complete,2)-1));

vessel_pair_names = strings(3,1);
vessel_pair_N = zeros(3,1);
vessel_median_difference = nan(3,1);
vessel_pair_p = nan(3,1);

for pair = 1:3
    first_region = pair_index(pair,1);
    second_region = pair_index(pair,2);
    valid_pair = isfinite(VesselPET(:,first_region)) & ...
                 isfinite(VesselPET(:,second_region));
    first_values = VesselPET(valid_pair,first_region);
    second_values = VesselPET(valid_pair,second_region);
    vessel_pair_names(pair) = ...
        string(vessel_names{first_region}) + " vs " + ...
        string(vessel_names{second_region});
    vessel_pair_N(pair) = nnz(valid_pair);
    vessel_median_difference(pair) = median(first_values-second_values);
    vessel_pair_p(pair) = signrank(first_values,second_values);
end

[sorted_p,sort_order] = sort(vessel_pair_p);
sorted_q = sorted_p .* 3 ./ (1:3)';
sorted_q = flipud(cummin(flipud(sorted_q)));
sorted_q(sorted_q > 1) = 1;
vessel_pair_q = nan(3,1);
vessel_pair_q(sort_order) = sorted_q;

fprintf(['\nArterial-territory Friedman test: N=%d, chi-square=%.3f, ' ...
    'p=%.3f, Kendall W=%.3f\n'], ...
    size(Vessel_complete,1),vessel_chi_square, ...
    p_vessel_omnibus,vessel_kendall_W);
vessel_stats_table = table( ...
    vessel_pair_names,vessel_pair_N,vessel_median_difference, ...
    vessel_pair_p,vessel_pair_q, ...
    'VariableNames',{'Comparison','N','MedianDifference','P','FDR_Q'});
disp(vessel_stats_table);

%% ============================================================
% 9. Main statistical figure
%    A+B: combined HBT associations and arterial matrix
%    C+D: horizontal mean bars with individual scatter
% =============================================================

load_paper;
set(gcf,'Color','w');

AB_position = [1:5 16:20 31:35 46:50 61:65];
C_position  = [6:10 21:25 36:40 51:55 66:70];
D_position  = [11:15 26:30 41:45 56:60 71:75];

%% ------------------------------------------------------------
% Panels A+B: network integrated with arterial correlation matrix
% -------------------------------------------------------------

ax1 = subplot(21,15,AB_position);
hold(ax1,'on');
axis(ax1,'off');

hippo_x = 0;
hippo_y = [3 2 1];

matrix_x = [2 3 4];
matrix_y = [3 2 1];

max_rho = max(abs(rho_HBT_vessel),[],'all','omitnan');
if isempty(max_rho) || isnan(max_rho) || max_rho == 0
    max_rho = 1;
end

% Association lines first, so matrix squares and circles remain on top.
for i = 1:3
    for j = 1:3

        current_rho = rho_HBT_vessel(i,j);
        current_q = q_HBT_vessel(i,j);

        if isnan(current_rho)
            continue
        end

        if current_q < 0.05
            line_style = '-';
        else
            line_style = '--';
        end

        normalized_r = abs(current_rho)/max_rho;
        line_color_index = 70+round(normalized_r*(230-70));
        line_color = blue(line_color_index,:);

        plot(ax1, ...
            [hippo_x+0.18 matrix_x(j)], ...
            [hippo_y(i) matrix_y(j)], ...
            'Color',line_color, ...
            'LineStyle',line_style, ...
            'LineWidth',1);
    end
end

% HBT circles.
scatter(ax1, ...
    repmat(hippo_x,3,1),hippo_y, ...
    300,[0.95 0.55 0.15],'filled', ...
    'MarkerEdgeColor',[0.2 0.2 0.2], ...
    'LineWidth',1);

for i = 1:3
    text(ax1,hippo_x-0.30,hippo_y(i),HBT_names{i}, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontSize',9);
end

% 3-by-3 arterial correlation matrix using the bar-plot blue range.
for row = 1:3
    for column = 1:3

        if row > column
            continue
        end

        current_rho = rho_vessel(row,column);
        color_fraction = max(0,min(1,current_rho));
        color_index = 70+round(color_fraction*(230-70));

        rectangle(ax1, ...
            'Position',[matrix_x(column)-0.5 matrix_y(row)-0.5 1 1], ...
            'FaceColor',blue(color_index,:), ...
            'EdgeColor',[0.85 0.85 0.85], ...
            'LineWidth',0.7);

        if row ~= column
            if q_vessel(row,column) < 0.05
                displayed_text = sprintf('%.3f*',current_rho);
            else
                displayed_text = sprintf('%.3f',current_rho);
            end

            if color_fraction > 0.55
                text_color = 'w';
            else
                text_color = 'k';
            end

            text(ax1,matrix_x(column),matrix_y(row),displayed_text, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'Color',text_color, ...
                'FontSize',9);
        end
    end
end

% Arterial circles occupy their corresponding diagonal matrix squares.
for j = 1:3
    scatter(ax1,matrix_x(j),matrix_y(j), ...
        300,green_skew,'filled', ...
        'MarkerEdgeColor',[0.2 0.2 0.2], ...
        'LineWidth',1);

    text(ax1,matrix_x(j),matrix_y(j),vessel_names{j}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',8, ...
        'FontWeight','bold');
end

xlim(ax1,[-0.9 4.6]);
ylim(ax1,[0.4 3.6]);
axis(ax1,'equal');

% No visible colorbar for the blue matrix.

%% ------------------------------------------------------------
% Panel C: HBT horizontal mean bars
% Head, Body, Tail from top to bottom; smallest x-value at left.
% -------------------------------------------------------------

ax3 = subplot(21,15,C_position);
hold(ax3,'on');

HBT_blue_indices = round(linspace(230,70,3));
HBT_bar_colors = blue(HBT_blue_indices,:);

rng(42);
scatter_size = 4;
jitter_width = 0.20;

for region = 1:3
    current_values = HBT_PET(:,region);
    current_values = current_values(isfinite(current_values));
    bar_mean = mean(current_values);
    bar_sem = std(current_values)/sqrt(numel(current_values));

    jitter = (rand(numel(current_values),1)-0.5)*jitter_width;

    scatter(ax3,current_values,region+jitter,scatter_size, ...
        HBT_bar_colors(region,:), ...
        'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeColor','none');

    barh(ax3,region,bar_mean,0.55, ...
        'FaceColor',HBT_bar_colors(region,:), ...
        'FaceAlpha',0.70, ...
        'EdgeColor','k', ...
        'LineWidth',0.5);

    % Mean +/- SEM at the bar endpoint.
    plot(ax3,[bar_mean-bar_sem bar_mean+bar_sem], ...
        [region region],'k-','LineWidth',1);
    plot(ax3,[bar_mean-bar_sem bar_mean-bar_sem], ...
        [region-0.08 region+0.08],'k-','LineWidth',1);
    plot(ax3,[bar_mean+bar_sem bar_mean+bar_sem], ...
        [region-0.08 region+0.08],'k-','LineWidth',1);
end

xline(ax3,0,'--', ...
    'Color',[0.5 0.5 0.5], ...
    'LineWidth',0.8);

ylim(ax3,[0.5 3.5]);
yticks(ax3,[]);
xlim(ax3,[-4 0]);
xticks(ax3,-4:1:0);

set(ax3, ...
    'YDir','reverse', ...
    'YAxisLocation','right', ...
    'XDir','normal');

xlabel(ax3,'Mean PET value');

ax3.FontSize = 10;
ax3.TickDir = 'out';
ax3.LineWidth = 1;
ax3.Box = 'off';

%% ------------------------------------------------------------
% Panel D: arterial-territory horizontal mean bars
% ACTP, MCAT, PCAT from top to bottom; smallest x-value at right.
% -------------------------------------------------------------

ax4 = subplot(21,15,D_position);
hold(ax4,'on');

vessel_blue_indices = round(linspace(230,70,3));
vessel_bar_colors = blue(vessel_blue_indices,:);

for region = 1:3
    current_values = VesselPET(:,region);
    current_values = current_values(isfinite(current_values));
    bar_mean = mean(current_values);
    bar_sem = std(current_values)/sqrt(numel(current_values));

    jitter = (rand(numel(current_values),1)-0.5)*jitter_width;

    scatter(ax4,current_values,region+jitter,scatter_size, ...
        vessel_bar_colors(region,:), ...
        'filled', ...
        'MarkerFaceAlpha',0.25, ...
        'MarkerEdgeColor','none');

    barh(ax4,region,bar_mean,0.55, ...
        'FaceColor',vessel_bar_colors(region,:), ...
        'FaceAlpha',0.70, ...
        'EdgeColor','k', ...
        'LineWidth',0.5);

    % Mean +/- SEM at the bar endpoint.
    plot(ax4,[bar_mean-bar_sem bar_mean+bar_sem], ...
        [region region],'k-','LineWidth',1);
    plot(ax4,[bar_mean-bar_sem bar_mean-bar_sem], ...
        [region-0.08 region+0.08],'k-','LineWidth',1);
    plot(ax4,[bar_mean+bar_sem bar_mean+bar_sem], ...
        [region-0.08 region+0.08],'k-','LineWidth',1);
end

xline(ax4,0,'--', ...
    'Color',[0.5 0.5 0.5], ...
    'LineWidth',0.8);

ylim(ax4,[0.5 3.5]);
yticks(ax4,[]);
xlim(ax4,[-4 0]);
xticks(ax4,-4:1:0);

set(ax4, ...
    'YDir','reverse', ...
    'YAxisLocation','left', ...
    'XDir','reverse');

xlabel(ax4,'Mean PET value');

ax4.FontSize = 10;
ax4.TickDir = 'out';
ax4.LineWidth = 1;
ax4.Box = 'off';
