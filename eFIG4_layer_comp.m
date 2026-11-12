clear;clc;

%% 1. load data, surface, color
load_color
load_surface

%% 2. correlation
load('/host/verges/tank/data/mo/atlas/scholtens2018neuroimage_economo.mat')
load('/host/verges/tank/data/mo/demo/02Parcel_feats_hmz_flip_asym_wPDB.mat')

% HCs
Parcel_PET_HCs = mean(Parcel_HCs_feats_hmz_flip_asym_wPDB.MidPET)';

X= {};
X{1,1} = ECONOMO_43.layerItotal_thickness_overall;
X{2,1} = ECONOMO_43.layerIItotal_thickness_overall;
X{3,1} = ECONOMO_43.layerIIItotal_thickness_overall;
X{4,1} = ECONOMO_43.layerIVtotal_thickness_overall;
X{5,1} = ECONOMO_43.layerVtotal_thickness_overall;
X{6,1} = ECONOMO_43.layerVItotal_thickness_overall;
X{7,1} = ECONOMO_43.total_thickness_overall;

all_HCs_uni_r = zeros(1,7); all_HCs_uni_p = zeros(1,7);
figure; hold on;
for i = 1:7
    subplot(1,7,i)
    X_temp = X{i,1}(allTempLbls); Y_temp = Parcel_PET_HCs(allTempLbls);
    valid = ~isnan(X_temp) & ~isnan(Y_temp);
    [r, p] = corr(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    scatter(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    lsline;
    title(sprintf('r=%.3f, p=%.4f', r, p));
    all_HCs_uni_p(1,i) = p;
    all_HCs_uni_r(1,i) = r;
end
FDR_p = mafdr(all_HCs_uni_p, 'BHFDR', true); display(FDR_p)

X{1,2} = ECONOMO_43.layerItotal_cell_content_cellsize;
X{2,2} = ECONOMO_43.layerIItotal_cell_content_cellsize;
X{3,2} = ECONOMO_43.layerIIItotal_cell_content_cellsize;
X{4,2} = ECONOMO_43.layerIVtotal_cell_content_cellsize;
X{5,2} = ECONOMO_43.layerVtotal_cell_content_cellsize;
X{6,2} = ECONOMO_43.layerVItotal_cell_content_cellsize;
X{7,2} = ECONOMO_43.total_cell_content_cellsize;

figure; hold on;
for i = 1:7
    subplot(1,7,i)
    X_temp = X{i,2}(allTempLbls); Y_temp = Parcel_PET_HCs(allTempLbls);
    valid = ~isnan(X_temp) & ~isnan(Y_temp);
    [r, p] = corr(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    scatter(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    lsline;
    title(sprintf('r=%.3f, p=%.4f', r, p));
    all_HCs_uni_p(1,i) = p;
    all_HCs_uni_r(2,i) = r;
end
FDR_p = mafdr(all_HCs_uni_p, 'BHFDR', true); display(FDR_p)

X{1,3} = ECONOMO_43.layerItotal_cell_content_mm3;
X{2,3} = ECONOMO_43.layerIItotal_cell_content_mm3;
X{3,3} = ECONOMO_43.layerIIItotal_cell_content_mm3;
X{4,3} = ECONOMO_43.layerIVtotal_cell_content_mm3;
X{5,3} = ECONOMO_43.layerVtotal_cell_content_mm3;
X{6,3} = ECONOMO_43.layerVItotal_cell_content_mm3;
X{7,3} = ECONOMO_43.total_cell_content_mm3;

figure; hold on;
for i = 1:7
    subplot(1,7,i)
    X_temp = X{i,3}(allTempLbls); Y_temp = Parcel_PET_HCs(allTempLbls);
    valid = ~isnan(X_temp) & ~isnan(Y_temp);
    [r, p] = corr(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    scatter(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    lsline;
    title(sprintf('r=%.3f, p=%.4f', r, p));
    all_HCs_uni_p(1,i) = p;
    all_HCs_uni_r(3,i) = r;
end
FDR_p = mafdr(all_HCs_uni_p, 'BHFDR', true); display(FDR_p)

for i = 1:7
    X_temp = [ ...
        X{i,1}(allTempLbls), ... % thickness
        X{i,2}(allTempLbls), ... % cell size
        X{i,3}(allTempLbls)  ... % cell density1
    ];
    Y_temp = Parcel_PET_HCs(allTempLbls);
    valid = all(~isnan(X_temp), 2) & ~isnan(Y_temp);

    [XL,~,XS,YS,beta,PCTVAR] = plsregress(zscore(X_temp(valid,:)), zscore(Y_temp(valid)),1);
    plsScore_HCs = XS(:,1);    
    [r,p] = corr(plsScore_HCs, zscore(Y_temp(valid)));
    all_HCs_mul_p(1,i) = p;
    all_HCs_mul_r(1,i) = r;
end
FDR_p = mafdr(all_HCs_mul_p, 'BHFDR', true); display(FDR_p)

% PTs
Parcel_PET_PTs = mean(Parcel_PTs_feats_hmz_flip_asym_wPDB.MidPET)';

X{1,1} = ECONOMO_43.layerItotal_thickness_overall;
X{2,1} = ECONOMO_43.layerIItotal_thickness_overall;
X{3,1} = ECONOMO_43.layerIIItotal_thickness_overall;
X{4,1} = ECONOMO_43.layerIVtotal_thickness_overall;
X{5,1} = ECONOMO_43.layerVtotal_thickness_overall;
X{6,1} = ECONOMO_43.layerVItotal_thickness_overall;
X{7,1} = ECONOMO_43.total_thickness_overall;

all_PTs_uni_r = zeros(1,7); all_PTs_uni_p = zeros(1,7);
figure; hold on;
for i = 1:7
    subplot(1,7,i)
    X_temp = X{i,1}(allTempLbls); Y_temp = Parcel_PET_PTs(allTempLbls);
    valid = ~isnan(X_temp) & ~isnan(Y_temp);
    [r, p] = corr(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    scatter(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    lsline;
    title(sprintf('r=%.3f, p=%.4f', r, p));
    all_PTs_uni_p(1,i) = p;
    all_PTs_uni_r(1,i) = r;
end
FDR_p = mafdr(all_PTs_uni_p, 'BHFDR', true); display(FDR_p)

X{1,2} = ECONOMO_43.layerItotal_cell_content_cellsize;
X{2,2} = ECONOMO_43.layerIItotal_cell_content_cellsize;
X{3,2} = ECONOMO_43.layerIIItotal_cell_content_cellsize;
X{4,2} = ECONOMO_43.layerIVtotal_cell_content_cellsize;
X{5,2} = ECONOMO_43.layerVtotal_cell_content_cellsize;
X{6,2} = ECONOMO_43.layerVItotal_cell_content_cellsize;
X{7,2} = ECONOMO_43.total_cell_content_cellsize;

figure; hold on;
for i = 1:7
    subplot(1,7,i)
    X_temp = X{i,2}(allTempLbls); Y_temp = Parcel_PET_PTs(allTempLbls);
    valid = ~isnan(X_temp) & ~isnan(Y_temp);
    [r, p] = corr(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    scatter(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    lsline;
    title(sprintf('r=%.3f, p=%.4f', r, p));
    all_PTs_uni_p(1,i) = p;
    all_PTs_uni_r(2,i) = r;
end
FDR_p = mafdr(all_PTs_uni_p, 'BHFDR', true); display(FDR_p)

X{1,3} = ECONOMO_43.layerItotal_cell_content_mm3;
X{2,3} = ECONOMO_43.layerIItotal_cell_content_mm3;
X{3,3} = ECONOMO_43.layerIIItotal_cell_content_mm3;
X{4,3} = ECONOMO_43.layerIVtotal_cell_content_mm3;
X{5,3} = ECONOMO_43.layerVtotal_cell_content_mm3;
X{6,3} = ECONOMO_43.layerVItotal_cell_content_mm3;
X{7,3} = ECONOMO_43.total_cell_content_mm3;

figure; hold on;
for i = 1:7
    subplot(1,7,i)
    X_temp = X{i,3}(allTempLbls); Y_temp = Parcel_PET_PTs(allTempLbls);
    valid = ~isnan(X_temp) & ~isnan(Y_temp);
    [r, p] = corr(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    scatter(zscore(X_temp(valid)), zscore(Y_temp(valid)));
    lsline;
    title(sprintf('r=%.3f, p=%.4f', r, p));
    all_PTs_uni_p(1,i) = p;
    all_PTs_uni_r(3,i) = r;
end
FDR_p = mafdr(all_PTs_uni_p, 'BHFDR', true); display(FDR_p)

for i = 1:7
    X_temp = [ ...
        X{i,1}(allTempLbls), ... % thickness
        X{i,2}(allTempLbls), ... % cell size
        X{i,3}(allTempLbls)  ... % cell density
    ];
    Y_temp = Parcel_PET_PTs(allTempLbls);
    valid = all(~isnan(X_temp), 2) & ~isnan(Y_temp);   

    [XL,~,XS,YS,beta,PCTVAR] = plsregress(zscore(X_temp(valid,:)), zscore(Y_temp(valid)),1);
    plsScore_PTs = XS(:,1);    
    [r,p] = corr(plsScore_PTs, zscore(Y_temp(valid)));
    all_PTs_mul_p(1,i) = p;
    all_PTs_mul_r(1,i) = r;
end
FDR_p = mafdr(all_PTs_mul_p, 'BHFDR', true); display(FDR_p)

%% plot correlation matrix
%  HCs
load_paper
subplot(10, 7, [4 5 11 12])
imagesc(all_HCs_uni_r(:,1:end-1)');
colormap(gca, gray);
axis equal tight;
xticks(1:4); yticks(1:6);

subplot(10, 7, [6 13])
imagesc(all_HCs_mul_r(:,1:end-1)')
clim([-1 1]);
colormap(gca, gray);
axis equal tight;
cb = colorbar;
xticks(1:4); yticks(1:6);

%  PTs
load_paper
subplot(10, 7, [4 5 11 12])
imagesc(all_PTs_uni_r(:,1:end-1)');
colormap(gca, blue);
axis equal tight;
xticks(1:4); yticks(1:6);

subplot(10, 7, [6 13])
imagesc(all_PTs_mul_r(:,1:end-1)');
clim([-1 1]);
colormap(gca, blue);
axis equal tight;
cb = colorbar;
xticks(1:4); yticks(1:6);

%% plot correlation of layers features
% for i = 1:100
%     x1(i) = rand*(10-20)+20;
%     y1(i) = rand*(15-25)+25;
%     z1(i) = rand*(30-40)+40;
% end
layer2Feats = [ECONOMO_43.layerIItotal_thickness_overall, ECONOMO_43.layerIItotal_cell_content_cellsize, ECONOMO_43.layerIItotal_cell_content_mm3];
z1 = zscore(layer2Feats(allTempLbls,1)') + 3;
y1 = zscore(layer2Feats(allTempLbls,2)');
x1 = zscore(layer2Feats(allTempLbls,3)');

[XX1,YY1,Z1]=griddata(x1,y1,z1,linspace(min(x1),max(x1))',linspace(min(y1),max(y1)),'cubic');

load_paper
subplot(10,7, [4 5 11 12])
mesh(XX1,YY1,Z1,'edgecolor','none','facecolor','interp');
hold on;
% mesh(XX1,YY1,2+0*Z1,Z1)
hold on;
contour(XX1,YY1,Z1,10)
colormap(slanCM('BuPu'))
cb = colorbar;
cb.Position = [0.2 0.5 0.02 0.1];
xlabel('Cell Intensity'); % X-axis label
ylabel('Cell Size'); % Y-axis label
zlabel('Layer Thickness'); % Z-axis label

fig = gcf;
fig.Renderer = 'painters';
exportgraphics(fig, 'figure.pdf', 'ContentType', 'vector');

%% plot correlation of lyrs and PET
%  remeber set i = 2, and rerun the pls_score
load_paper
ax = subplot(10,7,[6 7 13 14]);
corr_2grps_dist(zscore(plsScore_HCs),zscore(Parcel_PET_HCs(allTempLbls,1)),zscore(plsScore_PTs),zscore(Parcel_PET_PTs(allTempLbls,1)), ax)
[rho_HCs_ips, pval_HCs_ips] = corr(plsScore_HCs, Parcel_PET_HCs(allTempLbls,1));
[rho_PTs_ips, pval_PTs_ips] = corr(plsScore_PTs, Parcel_PET_PTs(allTempLbls,1));
