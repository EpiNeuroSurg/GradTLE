%% run Step 5D at first

%% 5. Ridge plot
x1 = dataX_clean(:,8);
x2 = dataX_clean(:,7);
x3 = dataX_clean(:,5);
x4 = dataX_clean(:,6);
y  = dataY_clean;

% Clean covariate to match dataX_clean
petP975_ExtAll_clean = petP975(validMask);
petP975_ExtAll_clean(187) = [];

data  = [x1, x2, x3, x4];
varNames = {'normHippoVol_asym','wmapVolSurgProp','surgVolProp','border'};

group_colors = [cmap_G2(128,:);
                cmap_mean(128,:);
                cmap_G1(128,:);
                blue(128,:)];

binW  = 0.01;
dotR  = 0.020;
dotSz = 4;
spread = 0.08;

load_paper
ax = subplot(10, 7, [1,2,8,9,15,16]);
hold on;
ax.FontSize  = 10;
ax.LineWidth = 1;
ax.XColor    = 'k';
ax.YColor    = 'k';
ax.YTick     = 1:4;
ax.YTickLabel = fliplr(varNames);
ylim([0.3 4.7]);
set(gca,'XScale','log');
xline(1,'--','LineWidth',1,'Color',[0.5 0.5 0.5]);

pred_prob = nan(size(y,1),4);
for v = 1:4
    ypos = 5 - v;
    vals_raw = data(:, v);
    vals     = zscore(vals_raw);
    col      = group_colors(v, :);

    % --- logistic regression with non-temporal hypometabolism covariate ---
    covar = petP975_ExtAll_clean;
    tbl = table(vals, covar, y, 'VariableNames', {'x','covar','y'});
    mdl = fitglm(tbl, 'y ~ x + covar', 'Distribution', 'binomial');
    beta  = mdl.Coefficients.Estimate(2);   % index 2 = x (predictor of interest)
    se    = mdl.Coefficients.SE(2);
    OR    = exp(beta);
    CI1   = exp(beta - 1.96*se);
    CI2   = exp(beta + 1.96*se);
    pred_prob(:,v) = mdl.Fitted.Probability;

    % Compute p-value from z-score
    zStat = beta / se;
    pVal  = 2 * (1 - normcdf(abs(zStat)));

    % --- center cloud on OR (in log scale) ---
    vals_plot = zscore(vals_raw) * spread + OR;

    % --- Half violin ---
    [fk, xi] = ksdensity(vals_plot, 'NumPoints', 256);
    fk = fk(:) / max(fk) * 0.35;
    xi = xi(:);
    py = [repmat(ypos, size(xi)); flipud(ypos + fk)];
    px = [xi; flipud(xi)];
    fill(px, py, col, 'FaceAlpha', 0.65, 'EdgeColor', col*0.7);

    % --- Stacked dots ---
    edges = (min(vals_plot)-binW) : binW : (max(vals_plot)+binW);
    for b = 1:numel(edges)-1
        inBin = vals_plot(vals_plot >= edges(b) & vals_plot < edges(b+1));
        nDots = numel(inBin);
        if nDots == 0, continue; end
        xCenter = mean([edges(b), edges(b+1)]);
        yDots   = ypos - dotR*(2*(1:nDots) - 1);
        scatter(repmat(xCenter,1,nDots), yDots, dotSz, col, 'filled', 'MarkerFaceAlpha', 0.80, 'MarkerEdgeAlpha', 0);
    end

    % --- OR + CI ---
    by = ypos + 0.0;
    plot([CI1 CI2], [by by], 'k-', 'LineWidth', 1);
    plot(OR, by, 'o', 'MarkerSize', 6, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k');

    fprintf('%-25s  OR = %.3f  95%% CI [%.3f, %.3f]  p = %.4f\n', varNames{v}, OR, CI1, CI2, pVal);
end

xlim([0.6 2.3]);
xticks([0.6 0.8 1 1.5 2]);
xticklabels({'0.6','0.8','1','1.5','2'});
xlabel('Odds Ratio');