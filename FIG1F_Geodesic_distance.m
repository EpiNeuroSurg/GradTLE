% clear;clc;
% 
% %% 1. load surface, color
% load_surface
% load_color
% load('/host/verges/tank/data/mo/demo/11idxHBI.mat')
% 
% %% 2. set origin and display contour
% % Find the 'Bridgehead'
% unfold_surf_path = "/host/verges/tank/data/mo/BIDS_20251201_TIANTAN_TLE/derivatives/hippunfold_v1.4.1/hippunfold/sub-MNI152NLin2009cAsym/surf/sub-MNI152NLin2009cAsym_hemi-L_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii";
% hipp_mni_surf_path = "/host/verges/tank/data/mo/software/hippunfold-templateflow-main/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii";
% fslr_surf_path = "/data_/mica1/01_programs/micapipe-v0.2.0/surfaces/fsLR-32k.L.surf.gii";
% bridgehead_axis = 2;           % MATLAB column index (1=x,2=y,3=z)
% bridgehead_thresh = -199.7;    % match your python code
% medial_dist_mm  =  2.5;
% lateral_dist_mm = 15.5;
% contour_tol_mm = 0.75; % Contour band tolerance (mm): picks vertices with |d - target| <= tol, % Increase if the band is too sparse/disconnected. 
% ap_axis = 2; % Posterior direction definition: In MNI, anterior-posterior is typically Y (2nd coord), anterior = larger Y.
% fprintf("Loading unfolded HippUnfold surface...\n");
% [U_V, ~] = load_gifti_surf(unfold_surf_path);
% fprintf("Loading MNI-space HippUnfold surface (for mapping to fsLR)...\n");
% [H_V, ~] = load_gifti_surf(hipp_mni_surf_path);
% fprintf("Loading fsLR-32k MNI-space surface...\n");
% [F_V, F_F] = load_gifti_surf(fslr_surf_path);
% fprintf("Unfolded mins: [%g %g %g]\n", min(U_V,[],1));
% fprintf("Unfolded maxs: [%g %g %g]\n", max(U_V,[],1));
% 
% is_TL = ismember(economo_parc, allTempLbls);
% %% 3. Build ROI-restricted geodesic 
% % Build full graph
% G_full = mesh_graph(F_V, F_F);
% % Induce subgraph on TL ROI ONLY
% roi_vidx = find(is_TL);
% G_roi = induced_subgraph(G_full, roi_vidx);  % has Nodes.OriginalIndex mapping
% % Map seed vertices into ROI:
% 
% 
% seed_in_roi1 = idxVtxSeed(is_TL(idxVtxSeed));
% seed_in_roi2 = find(economo_parc == 14);
% seed_in_roi = intersect(seed_in_roi1,seed_in_roi2);
% 
% % seed_in_roi = idxVtxSeed;
% if isempty(seed_in_roi)
%     error("No bridgehead seeds fall inside TL ROI. Check mapping or TL_labels.");
% end
% % Choose origin = most anterior bridgehead seed within ROI (Y max)
% [~, i0] = max(F_V(seed_in_roi, ap_axis));
% origin_vidx = seed_in_roi(i0);
% fprintf("Origin (bridgehead in TL ROI): vidx=%d, xyz=[%g %g %g]\n", ...
%     origin_vidx, F_V(origin_vidx,:));
% % Compute ROI-only distance: distance to nearest seed (multi-source) BUT constrained to TL
% dist_roi = nan(size(F_V,1),1);
% % Run multi-source Dijkstra on the ROI graph
% % Need ROI node ids for the seeds
% roi_orig = G_roi.Nodes.OriginalIndex;
% seed_roi_nodes = arrayfun(@(v)find(roi_orig==v,1), seed_in_roi);
% d_roi_nodes = multi_source_dijkstra(G_roi, seed_roi_nodes); % distances in ROI node space
% dist_roi(roi_orig) = d_roi_nodes;
% % Outside ROI set to NaN
% dist_roi(~is_TL) = NaN;
% fprintf("ROI distance stats: min=%g, max=%g (mm)\n", nanmin(dist_roi), nanmax(dist_roi));

clear;clc;

%% 1. load surface, color
load_surface
load_color

%% 2. set origin and display contour
% Find the 'Bridgehead'
unfold_surf_path = "/host/verges/tank/data/mo/BIDS_20251201_TIANTAN_TLE/derivatives/hippunfold_v1.4.1/hippunfold/sub-MNI152NLin2009cAsym/surf/sub-MNI152NLin2009cAsym_hemi-L_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii";
hipp_mni_surf_path = "/host/verges/tank/data/mo/software/hippunfold-templateflow-main/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii";
fslr_surf_path = "/data_/mica1/01_programs/micapipe-v0.2.0/surfaces/fsLR-32k.L.surf.gii";
bridgehead_axis = 2;           % MATLAB column index (1=x,2=y,3=z)
bridgehead_thresh = -199.7;    % match your python code
ap_axis = 2; % Posterior direction definition: In MNI, anterior-posterior is typically Y (2nd coord), anterior = larger Y.
fprintf("Loading unfolded HippUnfold surface...\n");
[U_V, ~] = load_gifti_surf(unfold_surf_path);
fprintf("Loading MNI-space HippUnfold surface (for mapping to fsLR)...\n");
[H_V, ~] = load_gifti_surf(hipp_mni_surf_path);
fprintf("Loading fsLR-32k MNI-space surface...\n");
[F_V, F_F] = load_gifti_surf(fslr_surf_path);
fprintf("Unfolded mins: [%g %g %g]\n", min(U_V,[],1));
fprintf("Unfolded maxs: [%g %g %g]\n", max(U_V,[],1));

bridgehead = find(U_V(:, bridgehead_axis) < bridgehead_thresh);
fprintf("Found %d bridgehead vertices in unfolded space.\n", numel(bridgehead));
if isempty(bridgehead)
    error("No bridgehead vertices found. Check threshold/axis or input surface.");
end

% ---- 2) Map bridgehead vertices to fsLR vertices (MNI space) -------
%   Vertex ordering is consistent between unfolded HippUnfold and MNI-space HippUnfold,
%   so "bridgehead" indices can be used to index H_V.
H_bridge_xyz = H_V(bridgehead, :);

% Map HippUnfold bridgehead xyz -> nearest fsLR vertex
fprintf("Mapping bridgehead points to fsLR vertices via nearest-neighbor...\n");
% knnsearch requires Statistics and Machine Learning Toolbox
fslr_bridge_vidx = knnsearch(F_V, H_bridge_xyz);  % returns fsLR vertex indices (1-based)

% Unique seed set (multi-source distances)
idxVtxSeed = unique(fslr_bridge_vidx);
fprintf("Mapped to %d unique fsLR seed vertices.\n", numel(idxVtxSeed));

bridge_labels = economo_parc(idxVtxSeed);
in_mask = ismember(bridge_labels, pureTempLbls);
idxVtxSeed_masked = idxVtxSeed(in_mask);
fprintf("Bridgehead vertices inside TL mask: %d / %d\n", sum(in_mask), numel(idxVtxSeed));

is_TL = ismember(economo_parc, pureTempLbls);

% 3) Build ROI-restricted geodesic 
% Build full graph
G_full = mesh_graph(F_V, F_F);
% Induce subgraph on TL ROI ONLY
roi_vidx = find(is_TL);
G_roi = induced_subgraph(G_full, roi_vidx);  % has Nodes.OriginalIndex mapping
% Map seed vertices into ROI:

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seed_in_roi = idxVtxSeed;
seed_in_roi1 = idxVtxSeed(is_TL(idxVtxSeed));
seed_in_roi2 = find(economo_parc == 14);
seed_in_roi = intersect(seed_in_roi1,seed_in_roi2)
% aa = false(64984, 1);          % <-- logical, not zeros()
% aa(seed_in_roi1) = true;
% [~, k] = max(F_V(aa, 2));         % column 2 = Y
% temp_vidx = find(aa);
% seed_in_roi = temp_vidx(k)

% temporal_mask = ismember(economo_parc, pureTempLbls);
% [~, k] = max(F_V(temporal_mask, 2));         % column 2 = Y
% temporal_vidx = find(temporal_mask);
% seed_in_roi = temporal_vidx(k);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(seed_in_roi)
    error("No bridgehead seeds fall inside TL ROI. Check mapping or TL_labels.");
end
% Choose origin = most anterior bridgehead seed within ROI (Y max)
[~, i0] = max(F_V(seed_in_roi, ap_axis));
origin_vidx = seed_in_roi(i0);
fprintf("Origin (bridgehead in TL ROI): vidx=%d, xyz=[%g %g %g]\n", ...
    origin_vidx, F_V(origin_vidx,:));

% Compute ROI-only distance: distance to nearest seed (multi-source) BUT constrained to TL
dist_roi = nan(size(F_V,1),1);
% Run multi-source Dijkstra on the ROI graph
% Need ROI node ids for the seeds
roi_orig = G_roi.Nodes.OriginalIndex;
seed_roi_nodes = arrayfun(@(v)find(roi_orig==v,1), seed_in_roi);
d_roi_nodes = multi_source_dijkstra(G_roi, seed_roi_nodes); % distances in ROI node space
dist_roi(roi_orig) = d_roi_nodes;
% Outside ROI set to NaN
dist_roi(~is_TL) = NaN;
fprintf("ROI distance stats: min=%g, max=%g (mm)\n", nanmin(dist_roi), nanmax(dist_roi));

% Contour levels (mm): choose step that makes sense for your resolution
contour_step_mm = 1.0;
max_dist_mm     = 150;     % or set [] to use max within ROI
contour_tol_mm  = 0.5;    % band half-width for averaging PET, e.g. 1mm-wide bands

%% 4. After you have dist_roi (NaN outside TL) and surf_left/right ======
%  4.1 Blank brain for distance
blank_brain = zeros(32492,1);
plot_hemispheres(blank_brain, inflated_surf_lh, 'views', 'lm');
colormap([1 1 1])

dist_plot = dist_roi;            % Nx1, concatenated L+R
dist_plot(~is_TL) = 0;

% --- Get axes handles from the figure created by plot_hemispheres ---
fig = gcf;
axs = findobj(fig, 'Type','axes');

% Some versions return multiple axes (L lat, L med, R lat, R med). We'll draw on ALL axes to be safe.
hold(axs, 'on');

% --- Define isoline levels you want ---
levels = [2.5 5 7.5 10 12.5 15.5];
max_level = 20;  % set to your expected max distance (adjust as needed)
min_level = 0;

% --- Build TL face mask (faces fully inside TL) ---
F = F_F; V = F_V;              % fsLR vertices/faces (concatenated)
face_tl = is_TL(F(:,1)) & is_TL(F(:,2)) & is_TL(F(:,3));
F_tl = F(face_tl,:);

% --- Draw each isoline as stitched polylines ---
for i = 1:numel(levels)
    L = levels(i);
    
    % Map L -> colormap index [1, 256]
    idx = round( (L - min_level) / (max_level - min_level) * 255 ) + 1;
    idx = max(1, min(256, idx));  % clamp
    line_color = red_flip(idx, :);
    
    % 1) triangle interpolation -> segments
    segs = mesh_isoline_segments(V, F_tl, dist_roi, L);
    % 2) stitch segments -> polylines
    polylines = stitch_isoline_segments(segs, 1e-3);
    % 3) keep only long curves
    polylines = keep_long_polylines(polylines, 30);
    % 4) smooth
    polylines = smooth_polylines(polylines, 7);
    % 5) plot on every axis
    for k = 1:numel(polylines)
        P = polylines{k};
        for a = 1:numel(axs)
            plot3(axs(a), P(:,1), P(:,2), P(:,3), '--', ...
                  'LineWidth', 1, 'Color', line_color);
        end
    end
end

% mark origin
oxyz = V(origin_vidx,:);
for a = 1:numel(axs)
    plot3(axs(a), oxyz(1), oxyz(2), oxyz(3), 'o', ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', red(1,:), ...
    'MarkerEdgeColor', [0 0 0], ...
    'LineWidth', 0.5);
end

%% distance plot
obj = plot_hemispheres(dist_plot(1:end/2,1), surf_lh, 'views', 'lm');
red_flip(1,:) = [1 1 1];
colormap(red_flip)

fig = gcf;
axs = findobj(fig, 'Type','axes');
hold(axs, 'on');
% mark origin
oxyz = V(origin_vidx,:);
for a = 1:numel(axs)
    plot3(axs(a), oxyz(1), oxyz(2), oxyz(3), 'o', ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', red(1,:), ...
    'MarkerEdgeColor', [0 0 0], ...
    'LineWidth', 0.5);
end