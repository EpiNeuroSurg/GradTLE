clear;
clc;
close all;

%% Step N5: bilateral arterial-territory masks
% Selected territory: Step2A blue
% Remaining cortical surface: Step2A gray
% Medial wall: white
% One figure is created for each territory.

ARTERIAL_DIR = ...
    '/host/verges/tank/data/mo/atlas/arterial/fsLR32k';

HBT_FILE = [ ...
    '/host/verges/tank/data/mo/atlas/hippocampus/HBT/' ...
    'MNI152NLin2009cAsym/' ...
    'HippoHBT_space-HippUnfold_den-2mm_hemi-L.label.gii'];

HIPPO_FOLDED_FILE = [ ...
    '/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/' ...
    'derivatives/hippunfold_v1.4.1/hippunfold/' ...
    'sub-MNI152NLin2009cAsym/surf/' ...
    'sub-MNI152NLin2009cAsym_hemi-L_space-T1w_' ...
    'den-2mm_label-hipp_midthickness.surf.gii'];

arterial_left_files = {
    fullfile(ARTERIAL_DIR,'ACTPL_hemi-L_fsLR32k.func.gii')
    fullfile(ARTERIAL_DIR,'MCATL_hemi-L_fsLR32k.func.gii')
    fullfile(ARTERIAL_DIR,'PCATL_hemi-L_fsLR32k.func.gii')
};

arterial_right_files = {
    fullfile(ARTERIAL_DIR,'ACTPR_hemi-R_fsLR32k.func.gii')
    fullfile(ARTERIAL_DIR,'MCATR_hemi-R_fsLR32k.func.gii')
    fullfile(ARTERIAL_DIR,'PCATR_hemi-R_fsLR32k.func.gii')
};

vessel_names = {'ACTP','MCAT','PCAT'};
HBT_names = {'Head','Body','Tail'};

if ~isfile(HBT_FILE)
    error('File not found:\n%s',HBT_FILE);
end
if ~isfile(HIPPO_FOLDED_FILE)
    error('File not found:\n%s',HIPPO_FOLDED_FILE);
end

for i = 1:numel(vessel_names)
    if ~isfile(arterial_left_files{i})
        error('File not found:\n%s',arterial_left_files{i});
    end
    if ~isfile(arterial_right_files{i})
        error('File not found:\n%s',arterial_right_files{i});
    end
end

load_surface;
load_color;
[mask_lh,mask_rh] = load_mask('midline');

mask_lh = logical(mask_lh(:));
mask_rh = logical(mask_rh(:));

if numel(mask_lh) ~= 32492 || numel(mask_rh) ~= 32492
    error('Each midline mask must contain 32,492 vertices.');
end

white_wall = [1 1 1];
gray_cortex = gray(16,:);
blue_mask = colors_cyto_nor(2,:);

for region = 1:numel(vessel_names)
    left_gii = gifti(arterial_left_files{region});
    right_gii = gifti(arterial_right_files{region});

    left_region = isfinite(left_gii.cdata(:)) & left_gii.cdata(:) > 0;
    right_region = isfinite(right_gii.cdata(:)) & right_gii.cdata(:) > 0;

    if numel(left_region) ~= 32492 || numel(right_region) ~= 32492
        error('%s masks must contain 32,492 vertices.',vessel_names{region});
    end

    % 0 = medial wall (white), 1 = cortex (gray), 2 = territory (blue)
    bilateral_display = ones(64984,1);
    bilateral_display(~[mask_lh;mask_rh]) = 0;
    bilateral_display([left_region;right_region] & [mask_lh;mask_rh]) = 2;

    % Do not call figure() here. plot_hemispheres creates the figure itself.
    plot_hemispheres( ...
        bilateral_display, ...
        {surf_lh,surf_rh}, ...
        'views','lm', ...
        'labeltext',vessel_names{region});

    region_figure = gcf;
    set(region_figure,'Color','w');

    brain_axes = findall(region_figure,'Type','axes');
    for a = 1:numel(brain_axes)
        colormap(brain_axes(a),[
            white_wall
            gray_cortex
            blue_mask
        ]);
        set(brain_axes(a),'CLim',[0 2]);
    end

    colormap(region_figure,[
        white_wall
        gray_cortex
        blue_mask
    ]);

    % Draw territory boundary directly, as in Step2A.
    outline_axes = findobj(region_figure,'Type','axes');
    outline_axes = outline_axes(~strcmp(get(outline_axes,'Tag'),'Colorbar'));

    for a = 1:numel(outline_axes)
        surface_patches = findobj(outline_axes(a),'Type','patch');

        for p = 1:numel(surface_patches)
            patch_faces = get(surface_patches(p),'Faces');
            patch_vertices = get(surface_patches(p),'Vertices');
            patch_values = get(surface_patches(p),'FaceVertexCData');

            if isempty(patch_faces) || isempty(patch_vertices) || ...
                    numel(patch_values) ~= size(patch_vertices,1)
                continue
            end

            selected_vertices = patch_values(:) == 2;
            mesh_edges = [patch_faces(:,[1 2]); ...
                          patch_faces(:,[2 3]); ...
                          patch_faces(:,[3 1])];

            is_boundary = selected_vertices(mesh_edges(:,1)) ~= ...
                          selected_vertices(mesh_edges(:,2));
            boundary_edges = unique( ...
                sort(mesh_edges(is_boundary,:),2),'rows');

            if isempty(boundary_edges)
                continue
            end

            v1 = patch_vertices(boundary_edges(:,1),:);
            v2 = patch_vertices(boundary_edges(:,2),:);

            X = [v1(:,1) v2(:,1) nan(size(v1,1),1)]';
            Y = [v1(:,2) v2(:,2) nan(size(v1,1),1)]';
            Z = [v1(:,3) v2(:,3) nan(size(v1,1),1)]';

            hold(outline_axes(a),'on');
            plot3(outline_axes(a),X(:),Y(:),Z(:), ...
                'k-','LineWidth',0.8);
        end
    end

    sgtitle(region_figure, ...
        sprintf('%s: bilateral lateral and medial views', ...
        vessel_names{region}), ...
        'FontSize',14, ...
        'FontWeight','bold');
end

%% Folded hippocampal HBT masks
% One figure with four panels: Head, Body, Tail, and black hippocampus.
% Selected HBT region = blue; remaining hippocampus = gray.

HBT_gii = gifti(HBT_FILE);
HBT_labels = double(HBT_gii.cdata(:));
hippo_folded = gifti(HIPPO_FOLDED_FILE);

if numel(HBT_labels) ~= 419 || size(hippo_folded.vertices,1) ~= 419
    error('HBT labels and folded hippocampal surface must have 419 vertices.');
end

figure('Color','w','Position',[100 100 1750 430]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');

for region = 1:3
    ax = nexttile;

    % Classify complete triangular faces to avoid the spiky vertex border.
    selected_vertices = HBT_labels == region;
    selected_faces = sum( ...
        selected_vertices(hippo_folded.faces),2) >= 2;

    % 1 = gray remaining hippocampus; 2 = selected blue HBT region
    HBT_face_display = 1 + double(selected_faces);

    patch(ax, ...
        'Faces',hippo_folded.faces, ...
        'Vertices',hippo_folded.vertices, ...
        'FaceVertexCData',HBT_face_display, ...
        'FaceColor','flat', ...
        'EdgeColor','none', ...
        'AmbientStrength',0.35, ...
        'DiffuseStrength',0.75, ...
        'SpecularStrength',0.15);

    axis(ax,'equal');
    axis(ax,'off');
    view(ax,[-120 25]);
    colormap(ax,[gray_cortex; blue_mask]);
    clim(ax,[0.5 2.5]);
    camlight(ax,'headlight');
    lighting(ax,'gouraud');

    % Draw HBT boundary directly, as in Step2A.
    mesh_edges = sort([hippo_folded.faces(:,[1 2]); ...
                       hippo_folded.faces(:,[2 3]); ...
                       hippo_folded.faces(:,[3 1])],2);

    edge_selected_face = repmat(selected_faces,3,1);
    [unique_edges,~,edge_group] = unique(mesh_edges,'rows');

    total_face_count = accumarray(edge_group,1);
    selected_face_count = accumarray( ...
        edge_group,double(edge_selected_face));

    % Include interfaces with gray faces and the external edge of the ROI.
    is_boundary = selected_face_count > 0 & ...
        (selected_face_count < total_face_count | total_face_count == 1);
    boundary_edges = unique_edges(is_boundary,:);

    v1 = hippo_folded.vertices(boundary_edges(:,1),:);
    v2 = hippo_folded.vertices(boundary_edges(:,2),:);

    X = [v1(:,1) v2(:,1) nan(size(v1,1),1)]';
    Y = [v1(:,2) v2(:,2) nan(size(v1,1),1)]';
    Z = [v1(:,3) v2(:,3) nan(size(v1,1),1)]';

    hold(ax,'on');
    plot3(ax,X(:),Y(:),Z(:),'k-','LineWidth',0.8);

    title(ax,HBT_names{region}, ...
        'FontSize',14, ...
        'FontWeight','bold');
end

%% Fourth panel: entirely black hippocampus
ax = nexttile;
black_hippocampus = ones(419,1);

patch(ax, ...
    'Faces',hippo_folded.faces, ...
    'Vertices',hippo_folded.vertices, ...
    'FaceVertexCData',black_hippocampus, ...
    'FaceColor','flat', ...
    'EdgeColor','none', ...
    'AmbientStrength',0.35, ...
    'DiffuseStrength',0.75, ...
    'SpecularStrength',0.15);

axis(ax,'equal');
axis(ax,'off');
view(ax,[-120 25]);
colormap(ax,[0 0 0]);
clim(ax,[0.5 1.5]);
camlight(ax,'headlight');
lighting(ax,'gouraud');

title(ax,'Whole hippocampus', ...
    'FontSize',14, ...
    'FontWeight','bold');