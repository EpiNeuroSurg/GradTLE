%  load surfaces and color
surface = fetch_template_surface('conte69', 'join', true);
[surf_lh, surf_rh] = load_conte69();

inflated_file_LH = gifti(fullfile('/host/verges/tank/data/mo/atlas/S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'));
inflated_file_RH = gifti(fullfile('/host/verges/tank/data/mo/atlas/S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii'));
inflated_surf_lh.tri = inflated_file_LH.faces; inflated_surf_lh.coord = inflated_file_LH.vertices'; 
inflated_surf_rh.tri = inflated_file_RH.faces; inflated_surf_rh.coord = inflated_file_RH.vertices'; 

very_inflated_file_LH = gifti(fullfile('/host/verges/tank/data/mo/atlas/S1200.L.very_inflated_MSMAll.32k_fs_LR.surf.gii'));
very_inflated_file_RH = gifti(fullfile('/host/verges/tank/data/mo/atlas/S1200.R.very_inflated_MSMAll.32k_fs_LR.surf.gii'));
very_inflated_surf_lh.tri = very_inflated_file_LH.faces; very_inflated_surf_lh.coord = very_inflated_file_LH.vertices'; 
very_inflated_surf_rh.tri = very_inflated_file_RH.faces; very_inflated_surf_rh.coord = very_inflated_file_RH.vertices'; 

lh = fullfile('/host/verges/tank/data/mo/atlas/lh.economo.32k_fsLR.label.gii');
rh = fullfile('/host/verges/tank/data/mo/atlas/rh.economo.32k_fsLR.label.gii');
test_lh = gifti(lh).cdata; test_rh = gifti(rh).cdata;
uni_lh = unique(test_lh);  uni_rh = unique(test_rh);

test_lh = test_lh - 1;  test_lh(test_lh == -1) = 0;
test_rh = test_rh + 42; test_rh(test_rh == 42) = 0;
economo_parc = single([test_lh; test_rh]);

economo_parc_34layers = readmatrix('/host/verges/tank/data/mo/software/ENIGMA-2.0.0/enigmatoolbox/datasets/parcellations/economo_koskinas_conte69_43regions.csv');
economo_parc_5layers  = readmatrix('/host/verges/tank/data/mo/software/ENIGMA-2.0.0/enigmatoolbox/datasets/parcellations/economo_koskinas_conte69.csv');

periTempLbls = [12 17 18];
mesTempLbls = [14 15 16];
latTempLbls = [37 38 39 40 41 42 43];
pureTempLbls = [mesTempLbls, latTempLbls];
allTempLbls = [periTempLbls, mesTempLbls, latTempLbls];

julich_lh = gifti(fullfile('/host/verges/tank/data/mo/atlas/lh.JulichBrainAtlas_3.1.fsLR32k.label.gii'));
julich_rh = gifti(fullfile('/host/verges/tank/data/mo/atlas/rh.JulichBrainAtlas_3.1.fsLR32k.label.gii'));
julich_parc1000 = [julich_lh.cdata;julich_rh.cdata*1000];
julich_parc = [julich_lh.cdata;julich_rh.cdata];