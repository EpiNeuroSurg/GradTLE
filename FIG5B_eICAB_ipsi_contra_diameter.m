%% S11: eICAB ipsilateral-versus-contralateral vessel measurements
% Seven displayed figures; no internal functions; nothing is saved.
clear; clc;

BIDS = ['/host/verges/tank' ...
    '/data/mo/BIDS_20251201_TLE_BTTH1'];
DERIV = fullfile(BIDS, 'derivatives', 'eICAB');
CLINICAL_FILE = fullfile(BIDS, 'dis_BTTH1_BFTH.csv');

vessel_names = ["ICA", "ACA-A1", "MCA-M1", "PCom", ...
                "PCA-P1", "PCA-P2", "SCA", "AChA"];
left_labels  = ["LCAR", "LACA1", "LMCA1", "LPcom", ...
                "LPCA1", "LPCA2", "LSCA", "LAChA"];
right_labels = ["RCAR", "RACA1", "RMCA1", "RPcom", ...
                "RPCA1", "RPCA2", "RSCA", "RAChA"];

metric_rows = ["Diameter median", "Diameter mean", "Volume"];
metric_titles = ["Median diameter", "Mean diameter", "Volume"];
metric_units = ["mm", "mm", "mm^3"];

summary_files = dir(fullfile(DERIV, 'sub-*', 'ses-pre', 'anat', ...
    'metric_space', '*_CircleOfWillis', '*_SUMMARY.csv'));
if numel(summary_files) ~= 107
    error('Expected 107 eICAB summaries, found %d.', numel(summary_files));
end

subject_id = strings(107, 1);
summary_path = strings(107, 1);
for i = 1:107
    summary_path(i) = fullfile(summary_files(i).folder, summary_files(i).name);
    subject_id(i) = string(regexp(summary_path(i), 'sub-[^/\\]+', 'match', 'once'));
    if strlength(subject_id(i)) == 0
        error('Cannot identify subject from %s.', summary_path(i));
    end
end
[subject_id, order] = sort(subject_id);
summary_path = summary_path(order);
if numel(unique(subject_id)) ~= 107
    error('Duplicate eICAB subject summaries detected.');
end

clinical = readtable(CLINICAL_FILE, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
if ~ismember('Side', clinical.Properties.VariableNames)
    error('Column "Side" is missing from %s.', CLINICAL_FILE);
end

best_id_column = "";
best_matches = -1;
for c = 1:width(clinical)
    candidate = strip(string(clinical{:, c}));
    candidate = regexprep(candidate, '\.0$', '');
    add_prefix = ~startsWith(candidate, 'sub-') & ~ismissing(candidate);
    candidate(add_prefix) = "sub-" + candidate(add_prefix);
    current_matches = sum(ismember(subject_id, candidate));
    if current_matches > best_matches
        best_matches = current_matches;
        best_id_column = string(clinical.Properties.VariableNames{c});
    end
end
if best_matches ~= 107
    error('Best clinical ID column "%s" matched %d/107 subjects.', ...
        best_id_column, best_matches);
end

clinical_id = strip(string(clinical{:, char(best_id_column)}));
clinical_id = regexprep(clinical_id, '\.0$', '');
add_prefix = ~startsWith(clinical_id, 'sub-') & ~ismissing(clinical_id);
clinical_id(add_prefix) = "sub-" + clinical_id(add_prefix);
clinical_side_raw = lower(strip(string(clinical{:, 'Side'})));

n_subjects = 107;
n_vessels = numel(vessel_names);
n_metrics = numel(metric_rows);
clinical_side = strings(n_subjects, 1);
ipsilateral_values = nan(n_subjects, n_vessels, n_metrics);
contralateral_values = nan(n_subjects, n_vessels, n_metrics);

for i = 1:n_subjects
    clinical_row = find(clinical_id == subject_id(i));
    if numel(clinical_row) ~= 1
        error('%s has %d clinical matches.', subject_id(i), numel(clinical_row));
    end
    side_value = clinical_side_raw(clinical_row);
    if ismember(side_value, ["l", "left", "lt"]) || contains(side_value, "left")
        clinical_side(i) = "L";
    elseif ismember(side_value, ["r", "right", "rt"]) || contains(side_value, "right")
        clinical_side(i) = "R";
    else
        error('Unrecognized Side "%s" for %s.', side_value, subject_id(i));
    end

    raw = readcell(summary_path(i));
    column_labels = strip(string(raw(1, 2:end)));
    row_labels = strip(string(raw(2:end, 1)));

    for v = 1:n_vessels
        left_column = find(column_labels == left_labels(v));
        right_column = find(column_labels == right_labels(v));
        if numel(left_column) ~= 1 || numel(right_column) ~= 1
            error('Missing %s/%s in %s.', left_labels(v), right_labels(v), summary_path(i));
        end
        for m = 1:n_metrics
            metric_row = find(row_labels == metric_rows(m));
            if numel(metric_row) ~= 1
                error('Missing "%s" in %s.', metric_rows(m), summary_path(i));
            end
            left_value = str2double(string(raw{metric_row + 1, left_column + 1}));
            right_value = str2double(string(raw{metric_row + 1, right_column + 1}));
            if clinical_side(i) == "L"
                ipsilateral_values(i, v, m) = left_value;
                contralateral_values(i, v, m) = right_value;
            else
                ipsilateral_values(i, v, m) = right_value;
                contralateral_values(i, v, m) = left_value;
            end
        end
    end
end

fprintf('Subjects: %d; left-sided: %d; right-sided: %d\n', ...
    n_subjects, sum(clinical_side == "L"), sum(clinical_side == "R"));

n_analyses = n_metrics * 2;
n_result_rows = n_analyses * n_vessels;
result_metric = strings(n_result_rows, 1);
result_outliers = false(n_result_rows, 1);
result_vessel = strings(n_result_rows, 1);
result_n_before = zeros(n_result_rows, 1);
result_n_removed = zeros(n_result_rows, 1);
result_n = zeros(n_result_rows, 1);
result_ipsi_mean = nan(n_result_rows, 1);
result_contra_mean = nan(n_result_rows, 1);
result_difference = nan(n_result_rows, 1);
result_ci_low = nan(n_result_rows, 1);
result_ci_high = nan(n_result_rows, 1);
result_p = nan(n_result_rows, 1);
result_dz = nan(n_result_rows, 1);
result_row = 0;

for analysis_index = 1:n_analyses
    metric_index = ceil(analysis_index / 2);
    remove_outliers = mod(analysis_index, 2) == 0;

    if remove_outliers
        figure_title = metric_titles(metric_index) + ...
            ": paired-difference 3-SD outlier removal";
    else
        figure_title = metric_titles(metric_index) + ": all paired data";
    end

    figure('Color', 'w', 'Position', [60 70 1500 720]);
    tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    for v = 1:n_vessels
        ipsi_all = ipsilateral_values(:, v, metric_index);
        contra_all = contralateral_values(:, v, metric_index);
        valid = isfinite(ipsi_all) & isfinite(contra_all);
        ipsi_all = ipsi_all(valid);
        contra_all = contra_all(valid);
        difference_all = ipsi_all - contra_all;
        n_before = numel(difference_all);

        paired_outlier = false(n_before, 1);
        if remove_outliers && n_before >= 3
            difference_sd_all = std(difference_all);
            if difference_sd_all > 0
                paired_outlier = abs(difference_all - mean(difference_all)) > ...
                    3 * difference_sd_all;
            end
        end

        ipsi = ipsi_all(~paired_outlier);
        contra = contra_all(~paired_outlier);
        difference = ipsi - contra;
        n_removed = sum(paired_outlier);
        n = numel(difference);

        mean_difference = mean(difference);
        difference_sd = std(difference);
        difference_se = difference_sd / sqrt(n);
        critical_t = tinv(0.975, n - 1);
        ci_low = mean_difference - critical_t * difference_se;
        ci_high = mean_difference + critical_t * difference_se;
        if difference_sd > 0
            [~, p_value] = ttest(ipsi, contra);
            cohen_dz = mean_difference / difference_sd;
        elseif mean_difference == 0
            p_value = 1;
            cohen_dz = 0;
        else
            p_value = 0;
            cohen_dz = sign(mean_difference) * Inf;
        end

        nexttile;
        hold on;
        for k = 1:n
            plot([1 2], [ipsi(k) contra(k)], '-', ...
                'Color', [0.84 0.84 0.84], 'LineWidth', 0.5);
        end
        scatter(ones(n, 1), ipsi, 13, [0.85 0.15 0.15], 'filled', ...
            'MarkerFaceAlpha', 0.5);
        scatter(2 * ones(n, 1), contra, 13, [0.10 0.35 0.80], 'filled', ...
            'MarkerFaceAlpha', 0.5);
        boxchart(ones(n, 1), ipsi, 'BoxFaceColor', [0.85 0.15 0.15], ...
            'MarkerStyle', 'none');
        boxchart(2 * ones(n, 1), contra, 'BoxFaceColor', [0.10 0.35 0.80], ...
            'MarkerStyle', 'none');
        xlim([0.5 2.5]);
        xticks([1 2]);
        xticklabels({'Ipsilateral', 'Contralateral'});
        ylabel(sprintf('eICAB %s (%s)', lower(metric_titles(metric_index)), ...
            metric_units(metric_index)));
        if remove_outliers
            title(sprintf('%s: n=%d, removed=%d, p=%.3g, d_z=%.2f', ...
                vessel_names(v), n, n_removed, p_value, cohen_dz), ...
                'FontWeight', 'bold');
        else
            title(sprintf('%s: paired n=%d, p=%.3g, d_z=%.2f', ...
                vessel_names(v), n, p_value, cohen_dz), 'FontWeight', 'bold');
        end
        grid on; box off;

        result_row = result_row + 1;
        result_metric(result_row) = metric_titles(metric_index);
        result_outliers(result_row) = remove_outliers;
        result_vessel(result_row) = vessel_names(v);
        result_n_before(result_row) = n_before;
        result_n_removed(result_row) = n_removed;
        result_n(result_row) = n;
        result_ipsi_mean(result_row) = mean(ipsi);
        result_contra_mean(result_row) = mean(contra);
        result_difference(result_row) = mean_difference;
        result_ci_low(result_row) = ci_low;
        result_ci_high(result_row) = ci_high;
        result_p(result_row) = p_value;
        result_dz(result_row) = cohen_dz;
    end
    sgtitle("eICAB " + figure_title, 'FontWeight', 'bold');
end

Results = table(result_metric, result_outliers, result_vessel, ...
    result_n_before, result_n_removed, result_n, result_ipsi_mean, ...
    result_contra_mean, result_difference, result_ci_low, result_ci_high, ...
    result_p, result_dz, ...
    'VariableNames', {'Metric', 'Outliers_removed', 'Vessel', 'N_before', ...
    'N_removed', 'N_analyzed', 'Ipsi_mean', 'Contra_mean', ...
    'Mean_difference', 'CI95_low', 'CI95_high', 'P_paired_t', 'Cohen_dz'});

disp('===== eICAB SIX-ANALYSIS SUMMARY =====');
disp(Results);
disp('Nothing was saved.');

%% Figure 7: selected diameter x volume comparisons after outlier removal
% Step4A manual boxplot style; no title/sgtitle objects are used.

load_color;
load_paper;
set(gcf,'Color','w');

selected_vessels = [8 3 5]; % AChA, MCA-M1, PCA-P1
selected_names = vessel_names(selected_vessels);
diameter_metric_index = 2; % Diameter mean, mm
volume_metric_index = 3;   % Volume, mm^3

figure7_ipsi = cell(3,1);
figure7_contra = cell(3,1);
figure7_dz = nan(3,1);
figure7_p = nan(3,1);
figure7_n = zeros(3,1);
figure7_removed = zeros(3,1);
figure7_y_limits = nan(3,2);

for panel = 1:3
    v = selected_vessels(panel);
    ipsi_all = ipsilateral_values(:,v,diameter_metric_index) .* ...
        ipsilateral_values(:,v,volume_metric_index);
    contra_all = contralateral_values(:,v,diameter_metric_index) .* ...
        contralateral_values(:,v,volume_metric_index);
    valid = isfinite(ipsi_all) & isfinite(contra_all);
    ipsi_all = ipsi_all(valid);
    contra_all = contra_all(valid);

    difference_all = ipsi_all-contra_all;
    paired_outlier = false(numel(difference_all),1);
    if numel(difference_all)>=3
        difference_sd_all = std(difference_all);
        if difference_sd_all>0
            paired_outlier = abs(difference_all-mean(difference_all)) > ...
                3*difference_sd_all;
        end
    end

    ipsi = ipsi_all(~paired_outlier);
    contra = contra_all(~paired_outlier);
    difference = ipsi-contra;

    figure7_ipsi{panel} = ipsi;
    figure7_contra{panel} = contra;
    figure7_n(panel) = numel(difference);
    figure7_removed(panel) = sum(paired_outlier);
    if std(difference)>0
        figure7_dz(panel) = mean(difference)/std(difference);
        [~,figure7_p(panel)] = ttest(ipsi,contra);
    elseif mean(difference)==0
        figure7_dz(panel) = 0;
        figure7_p(panel) = 1;
    else
        figure7_dz(panel) = sign(mean(difference))*Inf;
        figure7_p(panel) = 0;
    end

    panel_values = [ipsi;contra];
    value_range = max(panel_values)-min(panel_values);
    if value_range<=0
        value_range = max(abs(panel_values));
    end
    value_padding = max(0.10*value_range,eps);
    provisional_limits = [max(0,min(panel_values)-value_padding), ...
                          max(panel_values)+value_padding];
    limit_span = diff(provisional_limits);
    rounding_step = 10^floor(log10(limit_span));
    if limit_span/rounding_step<2
        rounding_step = rounding_step/5;
    elseif limit_span/rounding_step<5
        rounding_step = rounding_step/2;
    end
    figure7_y_limits(panel,:) = [ ...
        floor(provisional_limits(1)/rounding_step)*rounding_step, ...
        ceil(provisional_limits(2)/rounding_step)*rounding_step];
end

% Fixed negative Cohen's dz range requested for this figure.
cohen_y_limits = [-0.35 -0.01];

% Exact 21 x 15 subplot blocks (subplot indices are numbered row-wise).
cohen_position = [112:120 127:135]; % rows 8:9, columns 7:15
box_positions = {
    [37:39 52:54 67:69]  % rows 3:5, columns 7:9
    [40:42 55:57 70:72] % rows 3:5, columns 10:12
    [43:45 58:60 73:75] % rows 3:5, columns 13:15
};

%% Cohen's dz bar-and-line panel
ax_d = subplot(21,15,cohen_position);
hold(ax_d,'on');
bar(ax_d,1:3,figure7_dz,0.25, ...
    'FaceColor',blue(192,:), ...
    'FaceAlpha',0.75, ...
    'EdgeColor','k', ...
    'LineWidth',0.5);
plot(ax_d,1:3,figure7_dz,'-o', ...
    'Color',blue(256,:), ...
    'MarkerFaceColor',blue(256,:), ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1, ...
    'MarkerSize',3);
xlim(ax_d,[0.7 3.3]);
ylim(ax_d,cohen_y_limits);
xticks(ax_d,1:3);
xticklabels(ax_d,selected_names);
ylabel(ax_d,'Cohen''s d_z');
box(ax_d,'off');
ax_d.TickDir = 'out';
ax_d.FontSize = 7;

for panel = 1:3
    if figure7_p(panel)<0.001
        stat_label = sprintf('d_z=%.3f, P<0.001',figure7_dz(panel));
    else
        stat_label = sprintf('d_z=%.3f, P=%.3f', ...
            figure7_dz(panel),figure7_p(panel));
    end
    if figure7_dz(panel)>=0
        vertical_alignment = 'bottom';
        label_y = figure7_dz(panel)+0.06*diff(cohen_y_limits);
        label_y = min(label_y,cohen_y_limits(2)-0.03*diff(cohen_y_limits));
    else
        vertical_alignment = 'top';
        label_y = figure7_dz(panel)-0.06*diff(cohen_y_limits);
        label_y = max(label_y,cohen_y_limits(1)+0.03*diff(cohen_y_limits));
    end
    text(ax_d,panel,label_y,stat_label, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment',vertical_alignment, ...
        'FontSize',6);
end

%% Three paired custom boxplots
rng(42);
bw = 0.12;
ipsi_color = blue(192,:);
contra_color = gray(192,:);
solid_color = [0 0 0];

for panel = 1:3
    ax = subplot(21,15,box_positions{panel});
    hold(ax,'on');

    ipsi = figure7_ipsi{panel};
    contra = figure7_contra{panel};
    n = numel(ipsi);

    % Light-gray subject-level paired lines behind all other objects.
    for k = 1:n
        plot(ax,[1 2],[ipsi(k) contra(k)],'-', ...
            'Color',[0.84 0.84 0.84], ...
            'LineWidth',0.5);
    end

    for group = 1:2
        if group==1
            values = ipsi;
            color_group = ipsi_color;
        else
            values = contra;
            color_group = contra_color;
        end

        q1 = prctile(values,25);
        q3 = prctile(values,75);
        group_mean = mean(values);
        iqr_value = q3-q1;
        whisker_low = min(values(values>=q1-1.5*iqr_value));
        whisker_high = max(values(values<=q3+1.5*iqr_value));

        % Step4A scatter: size 4, alpha 0.25, no edge, bw*1.8 jitter.
        jitter = (rand(n,1)-0.5)*bw*1.8;
        scatter(ax,group+jitter,values,4,color_group,'filled', ...
            'MarkerFaceAlpha',0.25, ...
            'MarkerEdgeColor','none');

        % Q1-Q3 box.
        fill(ax,[group-bw group+bw group+bw group-bw], ...
            [q1 q1 q3 q3],color_group, ...
            'EdgeColor',solid_color, ...
            'LineWidth',0.5, ...
            'FaceAlpha',0.75);

        % Mean line.
        plot(ax,[group-bw group+bw],[group_mean group_mean],'-', ...
            'Color',solid_color,'LineWidth',2);

        % 1.5-IQR whiskers and caps.
        plot(ax,[group group],[whisker_low q1],'-', ...
            'Color',solid_color,'LineWidth',0.5);
        plot(ax,[group group],[q3 whisker_high],'-', ...
            'Color',solid_color,'LineWidth',0.5);
        plot(ax,[group-bw*0.5 group+bw*0.5], ...
            [whisker_low whisker_low],'-', ...
            'Color',solid_color,'LineWidth',0.5);
        plot(ax,[group-bw*0.5 group+bw*0.5], ...
            [whisker_high whisker_high],'-', ...
            'Color',solid_color,'LineWidth',0.5);
    end

    xlim(ax,[0.6 2.4]);
    panel_y_limits = figure7_y_limits(panel,:);
    ylim(ax,panel_y_limits);
    xticks(ax,[1 2]);
    xticklabels(ax,{'Ipsi','Contra'});
    ylabel(ax,'Mean diameter x volume (mm^4)');
    box(ax,'off');
    ax.TickDir = 'out';
    ax.FontSize = 7;

    % In-axes label does not alter subplot geometry.
    text(ax,1.5,panel_y_limits(2)-0.04*diff(panel_y_limits), ...
        sprintf('%s  n=%d  removed=%d',selected_names(panel), ...
        figure7_n(panel),figure7_removed(panel)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontWeight','bold', ...
        'FontSize',7);
end