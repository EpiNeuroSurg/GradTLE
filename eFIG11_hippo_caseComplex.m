% run after Step 7A2

%% Hippocampal Volume Analysis
load_paper;
load_color;

left_red = red_mean;
right_grey  = gray(16,:);

%--- Subplot 1: Group Comparison SF vs NSF ---%
subplot(10, 7, [1, 2, 3, 8, 9, 10, 15, 16, 17]);

hippo_asym = T.normHippoVol_asym;
sf_vals  = hippo_asym(Prog_bi_num == 1);
nsf_vals = hippo_asym(Prog_bi_num == 0);
sf_vals  = sf_vals(~isnan(sf_vals));
nsf_vals = nsf_vals(~isnan(nsf_vals));

[p_rank, ~] = ranksum(sf_vals, nsf_vals);

groups = [ones(size(sf_vals)); 2*ones(size(nsf_vals))];
vals   = [sf_vals; nsf_vals];

boxplot(vals, groups, ...
    'Labels', {'SF','NSF'}, ...
    'Colors', 'k', ...
    'Symbol', '', ...
    'Widths', 0.4);

hold on;

rng(42);
jitter_amt = 0.12;
x_sf  = 1 + (rand(size(sf_vals))  - 0.5) * 2 * jitter_amt;
x_nsf = 2 + (rand(size(nsf_vals)) - 0.5) * 2 * jitter_amt;
scatter(x_sf,  sf_vals,  30, left_red,  'filled', 'MarkerFaceAlpha', 0.65, ...
'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
scatter(x_nsf, nsf_vals, 30, right_grey,  'filled', 'MarkerFaceAlpha', 0.65, ...
'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
for grp = 1:2
    g_vals = vals(groups == grp);
    mu  = nanmean(g_vals);
    sem = nanstd(g_vals) / sqrt(sum(~isnan(g_vals)));
    errorbar(grp, mu, sem, 'k-', 'LineWidth', 2.5, 'CapSize', 10);
    plot(grp, mu, 'k.', 'MarkerSize', 18);
end
y_top = max(vals(~isnan(vals))) * 1.05;
line([1 2], [y_top y_top], 'Color', 'k', 'LineWidth', 1);
if p_rank < 0.001
    pstr = 'p < 0.001';
else
    pstr = sprintf('p = %.3f', p_rank);
end
text(1.5, y_top * 1.02, pstr, 'HorizontalAlignment', 'center', 'FontSize', 8);
ylabel('Hippocampal Volume Asymmetry Index', 'FontSize', 8);
title('Group comparison of normHippoVolAsymIdx', 'FontSize', 8);
xlim([0.5 2.5]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'SF','NSF'}, 'FontSize', 8, 'Box', 'off');
hold off;

% --- Cohen's d ---
n_sf  = sum(~isnan(sf_vals));
n_nsf = sum(~isnan(nsf_vals));
pooled_sd = sqrt(((n_sf - 1) * std(sf_vals)^2 + (n_nsf - 1) * std(nsf_vals)^2) / ...
                 (n_sf + n_nsf - 2));
cohens_d  = (mean(sf_vals) - mean(nsf_vals)) / pooled_sd;

fprintf('\n--- Group Comparison: SF vs NSF ---\n');
fprintf('SF  (n=%d): mean=%.4f, median=%.4f, SD=%.4f\n', numel(sf_vals), mean(sf_vals), median(sf_vals), std(sf_vals));
fprintf('NSF (n=%d): mean=%.4f, median=%.4f, SD=%.4f\n', numel(nsf_vals), mean(nsf_vals), median(nsf_vals), std(nsf_vals));
fprintf('Wilcoxon rank-sum: p = %.4f\n', p_rank);
fprintf('Cohen''s d: %.4f\n', cohens_d);

%--- Subplot 2: Hippocampal Asymmetry vs Time Since Surgery ---%
subplot(10, 7, [5, 6, 7, 12, 13, 14, 19, 20, 21]);

surgTime_all = surgTime;
valid = ~isnan(hippo_asym) & ~isnan(surgTime_all);
x_time  = surgTime_all(valid);
y_hippo = hippo_asym(valid);
prog_valid = Prog_bi_num(valid);

hold on;

yline(0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);

scatter(x_time(prog_valid==1), y_hippo(prog_valid==1), 30, left_red,  'filled', ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'SF');
scatter(x_time(prog_valid==0), y_hippo(prog_valid==0), 30, right_grey,  'filled', ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'NSF');

p_fit   = polyfit(x_time, y_hippo, 1);
x_range = linspace(min(x_time), max(x_time), 200);
y_fit   = polyval(p_fit, x_range);
plot(x_range, y_fit, 'k-', 'LineWidth', 2);

[r_pearson,  p_pearson]  = corr(x_time, y_hippo, 'Type', 'Pearson');
[r_spearman, p_spearman] = corr(x_time, y_hippo, 'Type', 'Spearman');

txt = sprintf('r = %.3f, p = %.3f\n\\rho = %.3f, p = %.3f', ...
    r_pearson, p_pearson, r_spearman, p_spearman);
text(0.05, 0.95, txt, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 8, ...
    'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.7 0.7 0.7]);

xlabel('Time Since Surgery (years)', 'FontSize', 8);
ylabel('Hippocampal Volume Asymmetry Index', 'FontSize', 8);
title('normHippoVolAsymIdx vs surgTime', 'FontSize', 8);
legend('SF', 'NSF', 'Location', 'best', 'FontSize', 8);
set(gca, 'FontSize', 8, 'Box', 'off');
hold off;

fprintf('\n--- Correlation: Hippocampal Asymmetry vs surgTime ---\n');
fprintf('Pearson:  r = %.4f, p = %.4f\n', r_pearson, p_pearson);
fprintf('Spearman: rho = %.4f, p = %.4f\n', r_spearman, p_spearman);