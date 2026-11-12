clear; clc;

%% BTTH1 MRI-positive/negative and PET-positive/negative comparisons
% Fully inline script: no local or user-defined functions.

load_surface
load_color

% Diverging colormap for Positive-minus-Negative contrast maps.
% The variables blue and red are supplied by load_color.
contrast_colormap = [flipud(blue); red];

CSV_FILE = '/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/dis_BTTH1_BFTH.csv';
PET_FILE = '/host/verges/tank/data/mo/demo/02Voxel_featsPET_hmz_flip_asym_wPDB.mat';

Demo_BTTH1 = readtable(CSV_FILE,'VariableNamingRule','preserve','TextType','string');
S = load(PET_FILE);
nBTTH1 = height(Demo_BTTH1);

required_vars = {'MRIDx','PETDx','AgeSurg','Sex','Duration','SZFreq'};
missing_vars = setdiff(required_vars,Demo_BTTH1.Properties.VariableNames);
if ~isempty(missing_vars)
    error('Missing CSV columns: %s',strjoin(missing_vars,', '));
end

if isfield(S,'Voxel_disPTs_neofeats_hmz_flip_asym_wPDB')
    MidPET_BTTH1 = S.Voxel_disPTs_neofeats_hmz_flip_asym_wPDB.MidPET;
elseif isfield(S,'Voxel_allPTs_neofeats_hmz_flip_asym_wPDB')
    MidPET_all = S.Voxel_allPTs_neofeats_hmz_flip_asym_wPDB.MidPET;
    MidPET_BTTH1 = MidPET_all(1:nBTTH1,:);
else
    error('No recognized PET structure found in %s.',PET_FILE);
end

if size(MidPET_BTTH1,1) ~= nBTTH1
    error('PET has %d rows but CSV has %d patients.',size(MidPET_BTTH1,1),nBTTH1);
end

fprintf('\n===== BTTH1 cohort =====\n');
fprintf('CSV patients: %d\n',nBTTH1);
fprintf('PET matrix: %d subjects x %d vertices\n',size(MidPET_BTTH1));

%% Convert MRIDx (1/positive = positive; 0/negative = negative)
MRI_raw = Demo_BTTH1.MRIDx;
MRI_group = strings(size(MRI_raw));
MRI_group(:) = missing;
if isnumeric(MRI_raw) || islogical(MRI_raw)
    MRI_group(double(MRI_raw)==1) = "Positive";
    MRI_group(double(MRI_raw)==0) = "Negative";
else
    MRI_text = lower(strtrim(string(MRI_raw)));
    MRI_text = regexprep(MRI_text,'[ _]','');
    MRI_pos = ismember(MRI_text,["1","positive","pos","+","mri+", ...
        "abnormal","lesional","visible","yes","y"]);
    MRI_neg = ismember(MRI_text,["0","negative","neg","-","mri-", ...
        "normal","nonlesional","non-lesional","invisible","no","n"]);
    MRI_group(MRI_pos) = "Positive";
    MRI_group(MRI_neg) = "Negative";
end

%% Convert PETDx
PET_raw = Demo_BTTH1.PETDx;
PET_group = strings(size(PET_raw));
PET_group(:) = missing;
if isnumeric(PET_raw) || islogical(PET_raw)
    PET_group(double(PET_raw)==1) = "Positive";
    PET_group(double(PET_raw)==0) = "Negative";
else
    PET_text = lower(strtrim(string(PET_raw)));
    PET_text = regexprep(PET_text,'[ _]','');
    PET_pos = ismember(PET_text,["1","positive","pos","+","pet+", ...
        "abnormal","hypometabolism","yes","y"]);
    PET_neg = ismember(PET_text,["0","negative","neg","-","pet-", ...
        "normal","nohypometabolism","no","n"]);
    PET_group(PET_pos) = "Positive";
    PET_group(PET_neg) = "Negative";
end

fprintf('\nMRI-positive: %d; MRI-negative: %d; missing/unrecognized: %d\n', ...
    sum(MRI_group=="Positive"),sum(MRI_group=="Negative"),sum(ismissing(MRI_group)));
fprintf('PET-positive: %d; PET-negative: %d; missing/unrecognized: %d\n', ...
    sum(PET_group=="Positive"),sum(PET_group=="Negative"),sum(ismissing(PET_group)));

%% Covariates
Age = double(Demo_BTTH1.AgeSurg);
Duration = double(Demo_BTTH1.Duration);
Sex = upper(strtrim(string(Demo_BTTH1.Sex)));
Sex(ismissing(Sex) | Sex=="") = missing;

if isnumeric(Demo_BTTH1.SZFreq) || islogical(Demo_BTTH1.SZFreq)
    SZFreq = double(Demo_BTTH1.SZFreq);
else
    SZFreq_text = lower(strtrim(string(Demo_BTTH1.SZFreq)));
    SZFreq = nan(size(SZFreq_text));
    SZFreq(SZFreq_text=="daily") = 1;
    SZFreq(SZFreq_text=="weekly") = 2;
    SZFreq(SZFreq_text=="monthly") = 3;
    SZFreq(SZFreq_text=="yearly") = 4;
    SZFreq_num = str2double(SZFreq_text);
    idx_num = isnan(SZFreq) & isfinite(SZFreq_num);
    SZFreq(idx_num) = SZFreq_num(idx_num);
end

%% Ipsilateral temporal mask
nV = size(MidPET_BTTH1,2);
nHalf = nV/2;
tempMask = false(1,nV);
tempMask(1:nHalf) = ismember(economo_parc(1:nHalf),allTempLbls);
FDR_ALPHA = 0.05;
MEAN_MAP_CLIM = [-4 0];  % identical range for Figures 1, 2, 5, and 6
fprintf('Temporal vertices used for FDR: %d\n',sum(tempMask));

%% ================= MRI analysis =================
keep_MRI = ~ismissing(MRI_group) & isfinite(Age) & isfinite(Duration) & ...
    isfinite(SZFreq) & ~ismissing(Sex);
MRI_group_use = MRI_group(keep_MRI);
MRI_Y = MidPET_BTTH1(keep_MRI,:);

fprintf('\n===== MRI-positive vs MRI-negative =====\n');
fprintf('Included: %d; positive: %d; negative: %d\n',sum(keep_MRI), ...
    sum(MRI_group_use=="Positive"),sum(MRI_group_use=="Negative"));

MRI_Group_Cov = term(cellstr(MRI_group_use));
MRI_Age_Cov = term(Age(keep_MRI));
MRI_Sex_Cov = term(cellstr(Sex(keep_MRI)));
MRI_Dur_Cov = term(Duration(keep_MRI));
MRI_Freq_Cov = term(SZFreq(keep_MRI));
MRI_Model = 1 + MRI_Group_Cov + MRI_Age_Cov + MRI_Sex_Cov + MRI_Dur_Cov + MRI_Freq_Cov;

slm_MRI = SurfStatLinMod(MRI_Y,MRI_Model);
slm_MRI = SurfStatT(slm_MRI,MRI_Group_Cov.Negative-MRI_Group_Cov.Positive);
MRI_p = 2*tcdf(-abs(slm_MRI.t),slm_MRI.df);
MRI_effect = slm_MRI.ef;
[MRI_h_temp,MRI_crit_p,MRI_q_temp] = fdr_bh(MRI_p(tempMask),FDR_ALPHA,'pdep','yes');
MRI_h = false(size(MRI_p)); MRI_h(tempMask) = logical(MRI_h_temp);
MRI_q = nan(size(MRI_p)); MRI_q(tempMask) = MRI_q_temp;

fprintf('Minimum uncorrected P: %.6g\n',min(MRI_p(tempMask),[],'omitnan'));
fprintf('Minimum FDR-corrected P: %.6g\n',min(MRI_q_temp,[],'omitnan'));
fprintf('FDR critical P: %.6g\n',MRI_crit_p);
fprintf('FDR-significant vertices: %d\n',sum(MRI_h));

MRI_pos_idx = MRI_group_use=="Positive";
MRI_neg_idx = MRI_group_use=="Negative";
MRI_mean_pos = mean(MRI_Y(MRI_pos_idx,:),1,'omitnan'); MRI_mean_pos(~tempMask)=0;
MRI_mean_neg = mean(MRI_Y(MRI_neg_idx,:),1,'omitnan'); MRI_mean_neg(~tempMask)=0;
MRI_temporal_mean_pos = mean(MRI_Y(MRI_pos_idx,tempMask),'all','omitnan');
MRI_temporal_mean_neg = mean(MRI_Y(MRI_neg_idx,tempMask),'all','omitnan');
fprintf('MRI POS mean temporal PET: %.3f\n',MRI_temporal_mean_pos);
fprintf('MRI NEG mean temporal PET: %.3f\n',MRI_temporal_mean_neg);
fprintf('MRI NEG - POS: %.3f\n',MRI_temporal_mean_neg-MRI_temporal_mean_pos);
MRI_effect_temp = MRI_effect; MRI_effect_temp(~tempMask)=0;
MRI_sig = zeros(size(MRI_effect)); MRI_sig(MRI_h)=-1;

plot_hemispheres(MRI_mean_pos(1:nHalf)',{surf_lh});
fig_MRI_pos = gcf; drawnow;
colormap(fig_MRI_pos,flipud(blue));
set(findall(fig_MRI_pos,'type','axes'),'CLim',MEAN_MAP_CLIM);
set(fig_MRI_pos,'Name','MRI-positive mean PET','NumberTitle','off');
annotation(fig_MRI_pos,'textbox',[0 0.94 1 0.05],'String','MRI-positive mean PET', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold','FontSize',12);
plot_hemispheres(MRI_mean_neg(1:nHalf)',{surf_lh});
fig_MRI_neg = gcf; drawnow;
colormap(fig_MRI_neg,flipud(blue));
set(findall(fig_MRI_neg,'type','axes'),'CLim',MEAN_MAP_CLIM);
set(fig_MRI_neg,'Name','MRI-negative mean PET','NumberTitle','off');
annotation(fig_MRI_neg,'textbox',[0 0.94 1 0.05],'String','MRI-negative mean PET', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold','FontSize',12);
plot_hemispheres(MRI_effect_temp(1:nHalf)',{surf_lh}); colormap(contrast_colormap);
MRI_lim=max(abs(MRI_effect_temp(tempMask)),[],'omitnan'); if MRI_lim==0, MRI_lim=1; end
set(findall(gcf,'type','axes'),'CLim',[-MRI_lim MRI_lim]);
set(gcf,'Name','MRI negative minus positive adjusted contrast','NumberTitle','off');
annotation(gcf,'textbox',[0 0.94 1 0.05], ...
    'String','MRI negative - positive: adjusted contrast', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold','FontSize',12);
plot_hemispheres(MRI_sig(1:nHalf)',{surf_lh}); colormap(flipud(gray));
set(findall(gcf,'type','axes'),'CLim',[-1 0]);
set(gcf,'Name','MRI comparison FDR-significant vertices','NumberTitle','off');
annotation(gcf,'textbox',[0 0.94 1 0.05], ...
    'String','MRI comparison: temporal P_{FDR} < 0.05', ...
    'Interpreter','tex','EdgeColor','none','HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',12);

Results.MRI.keep=keep_MRI; Results.MRI.t=slm_MRI.t; Results.MRI.df=slm_MRI.df;
Results.MRI.effect_negative_minus_positive=MRI_effect;
Results.MRI.p=MRI_p; Results.MRI.q=MRI_q;
Results.MRI.significant=MRI_h;

%% ================= PET analysis =================
keep_PET = ~ismissing(PET_group) & isfinite(Age) & isfinite(Duration) & ...
    isfinite(SZFreq) & ~ismissing(Sex);
PET_group_use = PET_group(keep_PET);
PET_Y = MidPET_BTTH1(keep_PET,:);

fprintf('\n===== PET-positive vs PET-negative =====\n');
fprintf('Included: %d; positive: %d; negative: %d\n',sum(keep_PET), ...
    sum(PET_group_use=="Positive"),sum(PET_group_use=="Negative"));

PET_Group_Cov = term(cellstr(PET_group_use));
PET_Age_Cov = term(Age(keep_PET));
PET_Sex_Cov = term(cellstr(Sex(keep_PET)));
PET_Dur_Cov = term(Duration(keep_PET));
PET_Freq_Cov = term(SZFreq(keep_PET));
PET_Model = 1 + PET_Group_Cov + PET_Age_Cov + PET_Sex_Cov + PET_Dur_Cov + PET_Freq_Cov;

slm_PET = SurfStatLinMod(PET_Y,PET_Model);
slm_PET = SurfStatT(slm_PET,PET_Group_Cov.Negative-PET_Group_Cov.Positive);
PET_p = 2*tcdf(-abs(slm_PET.t),slm_PET.df);
PET_effect = slm_PET.ef;
[PET_h_temp,PET_crit_p,PET_q_temp] = fdr_bh(PET_p(tempMask),FDR_ALPHA,'pdep','yes');
PET_h = false(size(PET_p)); PET_h(tempMask) = logical(PET_h_temp);
PET_q = nan(size(PET_p)); PET_q(tempMask) = PET_q_temp;

fprintf('Minimum uncorrected P: %.6g\n',min(PET_p(tempMask),[],'omitnan'));
fprintf('Minimum FDR-corrected P: %.6g\n',min(PET_q_temp,[],'omitnan'));
fprintf('FDR critical P: %.6g\n',PET_crit_p);
fprintf('FDR-significant vertices: %d\n',sum(PET_h));

PET_pos_idx = PET_group_use=="Positive";
PET_neg_idx = PET_group_use=="Negative";
PET_mean_pos = mean(PET_Y(PET_pos_idx,:),1,'omitnan'); PET_mean_pos(~tempMask)=0;
PET_mean_neg = mean(PET_Y(PET_neg_idx,:),1,'omitnan'); PET_mean_neg(~tempMask)=0;
PET_temporal_mean_pos = mean(PET_Y(PET_pos_idx,tempMask),'all','omitnan');
PET_temporal_mean_neg = mean(PET_Y(PET_neg_idx,tempMask),'all','omitnan');
fprintf('PET POS mean temporal PET: %.3f\n',PET_temporal_mean_pos);
fprintf('PET NEG mean temporal PET: %.3f\n',PET_temporal_mean_neg);
fprintf('PET NEG - POS: %.3f\n',PET_temporal_mean_neg-PET_temporal_mean_pos);
PET_effect_temp = PET_effect; PET_effect_temp(~tempMask)=0;
PET_sig = zeros(size(PET_effect)); PET_sig(PET_h)=-1;

plot_hemispheres(PET_mean_pos(1:nHalf)',{surf_lh});
fig_PET_pos = gcf; drawnow;
colormap(fig_PET_pos,flipud(blue));
set(findall(fig_PET_pos,'type','axes'),'CLim',MEAN_MAP_CLIM);
set(fig_PET_pos,'Name','PET-positive mean PET','NumberTitle','off');
annotation(fig_PET_pos,'textbox',[0 0.94 1 0.05],'String','PET-positive mean PET', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold','FontSize',12);
plot_hemispheres(PET_mean_neg(1:nHalf)',{surf_lh});
fig_PET_neg = gcf; drawnow;
colormap(fig_PET_neg,flipud(blue));
set(findall(fig_PET_neg,'type','axes'),'CLim',MEAN_MAP_CLIM);
set(fig_PET_neg,'Name','PET-negative mean PET','NumberTitle','off');
annotation(fig_PET_neg,'textbox',[0 0.94 1 0.05],'String','PET-negative mean PET', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold','FontSize',12);
plot_hemispheres(-PET_effect_temp(1:nHalf)',{surf_lh}); colormap(contrast_colormap);
PET_lim=max(abs(PET_effect_temp(tempMask)),[],'omitnan'); if PET_lim==0, PET_lim=1; end
set(findall(gcf,'type','axes'),'CLim',[-PET_lim PET_lim]);
set(gcf,'Name','PET negative minus positive adjusted contrast','NumberTitle','off');
annotation(gcf,'textbox',[0 0.94 1 0.05], ...
    'String','PET negative - positive: adjusted contrast', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold','FontSize',12);
plot_hemispheres(PET_sig(1:nHalf)',{surf_lh}); colormap(flipud(gray));
set(findall(gcf,'type','axes'),'CLim',[-1 0]);
set(gcf,'Name','PET comparison FDR-significant vertices','NumberTitle','off');
annotation(gcf,'textbox',[0 0.94 1 0.05], ...
    'String','PET comparison: temporal P_{FDR} < 0.05', ...
    'Interpreter','tex','EdgeColor','none','HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',12);

Results.PET.keep=keep_PET; Results.PET.t=slm_PET.t; Results.PET.df=slm_PET.df;
Results.PET.effect_negative_minus_positive=PET_effect;
Results.PET.p=PET_p; Results.PET.q=PET_q;
Results.PET.significant=PET_h;
