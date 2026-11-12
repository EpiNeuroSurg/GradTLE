clear;clc; 

%% 1. load surface
load_demo
load_surface
load_color

%% 2. demo
load('/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat')

Age_PTs    = Demo_disPTs.AgeSurg;   
Sex_PTs    = Demo_disPTs.Sex;     Sex_num = double(strcmp(Sex_PTs,'M'));
Dur_PTs    = Demo_disPTs.Duration;
SZFreq_PTs = Demo_disPTs.SZFreq;   SZFreq_PTs_num = zeros(size(SZFreq_PTs));
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Daily')) = 1;   SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Weekly')) = 2;
SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Monthly')) = 3; SZFreq_PTs_num(strcmp(SZFreq_PTs, 'Yearly')) = 4;
Prog_PTs   = Demo_disPTs.x1YProgI; Prog_bi = repmat({'NSF'}, size(Prog_PTs)); Prog_bi(strcmp(Prog_PTs, 'ILAE 1')) = {'SF'};

Grp_Cov = term(Prog_bi);
Age_Cov = term(Age_PTs); Sex_Cov    = term(Sex_PTs);
Dur_Cov = term(Dur_PTs); SZFreq_Cov = term(SZFreq_PTs_num);

%% 3. Surface-wise stats of PET between NSF and SF
MidPET_PTs     = Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET;
MidPET_PTs_SF  = mean(MidPET_PTs(strcmp(Prog_bi, 'SF'),:),1);
MidPET_PTs_NSF = mean(MidPET_PTs(strcmp(Prog_bi, 'NSF'),:),1);

Model = 1 + Grp_Cov + Age_Cov + Sex_Cov + Dur_Cov + SZFreq_Cov;

slm_prog = SurfStatLinMod(MidPET_PTs, Model);
slm_prog = SurfStatT(slm_prog, Grp_Cov.SF - Grp_Cov.NSF); 

p_prog = tcdf(slm_prog.t,slm_prog.df);
[h_prog, crit_p_prog, adj_p_prog] = fdr_bh(p_prog, 0.05, 'pdep', 'yes');

D_prog = slm_prog.ef;

% 3.1 display SF and NSF respectively
MidPET_PTs_SF(economo_parc == 0) = max(MidPET_PTs_SF);
plot_hemispheres(MidPET_PTs_SF(1,1:end/2)', {surf_lh}); colormap(flipud(blue)); 
set(findall(gcf,'type','axes'),'CLim',[-6 0]);

MidPET_PTs_NSF(economo_parc == 0) = max(MidPET_PTs_NSF);
plot_hemispheres(MidPET_PTs_NSF(1,1:end/2)', {surf_lh}); colormap(flipud(blue)); 
set(findall(gcf,'type','axes'),'CLim',[-6 0]);

% raw cohen's D, SF - NSF
sigMask_Prog_uncorr = (p_prog < 0.05);
D_Prog_uncorr = D_prog .* sigMask_Prog_uncorr;
plot_hemispheres(D_Prog_uncorr(1,1:end/2)', {surf_lh}); colormap(flipud(blue)); 
set(findall(gcf,'type','axes'),'CLim',[-6 0]);

% after FDR, SF - NSF
sigMask_Prog_FDR = (h_prog == 1);
D_Prog_FDR = D_prog .* sigMask_Prog_FDR;
plot_hemispheres(D_Prog_FDR(1,1:end/2)', {surf_lh}); colormap(flipud(gray));
set(findall(gcf, 'type', 'axes'), 'CLim', [-1 0]);

%% 4. Surface-wise stats pf surigcal cavity between NSF and SF
PATH_PTs = '/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/';

%  4.1 load surgMask data
surgMask_PTs   = zeros(length(Demo_disPTs.ID),64984);
for i = 1:length(Demo_disPTs.ID)
    File_PTs_LH_surgMask = strcat(PATH_PTs,'surgMask/',Demo_disPTs.ID{i},'/ses-post/',Demo_disPTs.ID{i},'_ses-post_hemi-L_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    File_PTs_RH_surgMask = strcat(PATH_PTs,'surgMask/',Demo_disPTs.ID{i},'/ses-post/',Demo_disPTs.ID{i},'_ses-post_hemi-R_space-nativepro_surf-fsLR-32k_label-cavity.func.gii');
    surgMask_raw = [gifti(File_PTs_LH_surgMask).cdata', gifti(File_PTs_RH_surgMask).cdata'];

    if strcmp(Demo_disPTs.Side{i}, 'R')
        surgMask_temp = [surgMask_raw(1, end/2+1:end), surgMask_raw(1, 1:end/2)];
    else
        surgMask_temp = surgMask_raw;
    end
    
    in_mask = (ismember(economo_parc, allTempLbls))';  
    surgMask_temp(~in_mask) = 0;               
    surgMask_PTs(i,:) = surgMask_temp;
end

surgMask_PTs_SF  = mean(surgMask_PTs(strcmp(Prog_bi, 'SF'),:),1);
surgMask_PTs_NSF = mean(surgMask_PTs(strcmp(Prog_bi, 'NSF'),:),1);

Model = 1 + Grp_Cov + Age_Cov + Sex_Cov + Dur_Cov + SZFreq_Cov;

slm_surgMask = SurfStatLinMod(surgMask_PTs, Model);
slm_surgMask = SurfStatT(slm_surgMask, Grp_Cov.SF - Grp_Cov.NSF); 

p_surgMask = tcdf(slm_surgMask.t,slm_surgMask.df);
[h_surgMask, crit_p_surgMask, adj_p_surgMask] = fdr_bh(p_surgMask, 0.05, 'pdep', 'yes');

D_prog = slm_surgMask.ef;

% 3.1 display SF and NSF respectively
surgMask_PTs_SF(economo_parc == 0) = min(surgMask_PTs_SF);
plot_hemispheres(surgMask_PTs_SF(1,1:end/2)', {surf_lh}); colormap(cmap_purple); 
set(findall(gcf,'type','axes'),'CLim',[0 1]);

surgMask_PTs_NSF(economo_parc == 0) = min(surgMask_PTs_NSF);
plot_hemispheres(surgMask_PTs_NSF(1,1:end/2)', {surf_lh}); colormap(cmap_purple); 
set(findall(gcf,'type','axes'),'CLim',[0 1]);

% raw cohen's D, SF - NSF
sigMask_Prog_uncorr = (D_prog < 0.05);
D_Prog_uncorr = D_prog .* sigMask_Prog_uncorr;
plot_hemispheres(D_Prog_uncorr(1,1:end/2)', {surf_lh}); colormap(cmap_purple); 
set(findall(gcf,'type','axes'),'CLim',[0 1]);

% after FDR, SF - NSF
sigMask_Prog_FDR = (h_surgMask == 1);
D_Prog_FDR = D_prog .* sigMask_Prog_FDR;
plot_hemispheres(D_Prog_FDR(1,1:end/2)', {surf_lh}); colormap(flipud(cmap_red));

%% 5. Group comparison: SF vs NSF (subject-level)
load_paper;
left_blue = blue(200,:);
right_col = gray(200,:);
left_alpha  = 1.0;
right_alpha = 0.30;
cx=1; v_off=0.2; s_off=0.2; box_w=0.08; v_w=0.15; lw=0.5;

subplots = [7, 14];
titles   = {'MidPET', 'Surgical Cavity'};

% -------- subject-level summaries --------
PET_data  = Voxel_PTs_neofeats_hmz_flip_asym_wPDB.MidPET(:, 1:end/2); % 263 x 32492
Surg_data = surgMask_PTs(:, 1:end/2);

subj_PET  = NaN(size(PET_data,1),  1);
subj_Surg = NaN(size(Surg_data,1), 1);
for s = 1:size(PET_data,1)
    v = PET_data(s,:);
    if any(v ~= 0), subj_PET(s) = mean(v(v ~= 0)); end

    v = Surg_data(s,:);
    if any(v ~= 0), subj_Surg(s) = sum(v == 1) / numel(v); end
end

subj_Surg(subj_Surg < 0.02) = NaN;
data_all = {subj_PET, subj_Surg};

% for k = 1:2
%     subplot(10,7,subplots(k)); hold on;
% 
%     d_SF  = data_all{k}(strcmp(Prog_bi, 'SF'));
%     d_NSF = data_all{k}(strcmp(Prog_bi, 'NSF'));
% 
%     % remove NaN
%     d_SF  = d_SF(~isnan(d_SF));
%     d_NSF = d_NSF(~isnan(d_NSF));
% 
%     % -------- stats --------
%     [~, p] = ttest2(d_SF, d_NSF);
%     n1 = numel(d_SF); n2 = numel(d_NSF);
%     s1 = std(d_SF);   s2 = std(d_NSF);
%     pooledSD = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2));
%     cohen_d  = (mean(d_SF) - mean(d_NSF)) / pooledSD;
% 
%     % -------- left violin+box: SF --------
%     vx = cx - v_off;
%     [f,xi] = ksdensity(d_SF,'Bandwidth',0.3); f = v_w*(f/max(f));
%     fill([vx-f, vx*ones(1,numel(xi))],[xi, fliplr(xi)], left_blue, ...
%         'FaceAlpha',left_alpha,'EdgeColor','k','LineWidth',lw,'EdgeAlpha',1);
%     bx = cx - s_off; q = prctile(d_SF,[25 50 75]); iqr_v = q(3)-q(1);
%     wlo = min(d_SF(d_SF >= q(1)-1.5*iqr_v)); whi = max(d_SF(d_SF <= q(3)+1.5*iqr_v));
%     fill([bx-box_w bx+box_w bx+box_w bx-box_w],[q(1) q(1) q(3) q(3)], left_blue, ...
%         'FaceAlpha',left_alpha,'EdgeColor','k','LineWidth',lw,'EdgeAlpha',1);
%     plot([bx-box_w bx+box_w],[q(2) q(2)],'-','Color','w','LineWidth',lw);
%     plot([bx bx],[wlo q(1)],'k-','LineWidth',lw);
%     plot([bx bx],[q(3) whi],'k-','LineWidth',lw);
% 
%     % -------- right violin+box: NSF --------
%     vx = cx + v_off;
%     [f,xi] = ksdensity(d_NSF,'Bandwidth',0.3); f = v_w*(f/max(f));
%     fill([vx+f, vx*ones(1,numel(xi))],[xi, fliplr(xi)], right_col, ...
%         'FaceAlpha',right_alpha,'EdgeColor','k','LineWidth',lw,'EdgeAlpha',1);
%     bx = cx + s_off; q = prctile(d_NSF,[25 50 75]); iqr_v = q(3)-q(1);
%     wlo = min(d_NSF(d_NSF >= q(1)-1.5*iqr_v)); whi = max(d_NSF(d_NSF <= q(3)+1.5*iqr_v));
%     fill([bx-box_w bx+box_w bx+box_w bx-box_w],[q(1) q(1) q(3) q(3)], right_col, ...
%         'FaceAlpha',right_alpha,'EdgeColor','k','LineWidth',lw,'EdgeAlpha',1);
%     plot([bx-box_w bx+box_w],[q(2) q(2)],'-','Color','w','LineWidth',lw);
%     plot([bx bx],[wlo q(1)],'k-','LineWidth',lw);
%     plot([bx bx],[q(3) whi],'k-','LineWidth',lw);
% 
%     % -------- axes & annotation --------
%     set(gca,'XTick',[cx-v_off cx+v_off],'XTickLabel',{'SF','NSF'}, ...
%         'FontSize',5,'Box','off');
%     xlim([cx-v_off-v_w-0.08, cx+v_off+v_w+0.08]);
%     title(titles{k},'FontSize',6);
%     text(0.02, 0.95, sprintf('p=%.3e, d=%.3f', p, cohen_d), ...
%         'Units','normalized','FontSize',6,'VerticalAlignment','top');
% end
