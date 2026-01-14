%% IT Neuron Archetype Analysis - Main Script
% This script consolidates ParTI fitting and 3D visualization for 
% cross-species IT neuron archetype analysis (Figure 2G-I, K, O-P, S3A-B)
%
% Prerequisites:
%   - ParTI toolbox: https://github.com/AlonLabWIS/ParTI (put this on the path)
%   - CSV exports from R notebook: 2GHIJKOP_OpossumMouse_ITArchetypes.Rmd
%
% Directory Structure:
%   opossum-V1-omics/figure 2/data/
%     - 2G_Mouse_SubITPCEmbeddings.csv
%     - 2G_Mouse_SubProjITPCEmbeddings.csv
%     - 2H_Opossum_SubITPCEmbeddings.csv
%     - 2H_Opossum_SubProjITPCEmbeddings.csv
%     - 2I_OpossumMouse_MouseSubProjITPCEmbeddings.csv
%
% Output:
%   - *_Data.mat files containing arcOrig (tetrahedron vertices)
%   - 3D scatter plots with tetrahedron overlays and silhouette heatmaps

%% Configuration
addpath(genpath(fullfile("E:/Ryan/GitHub/ParTI/")));

data_dir = 'E:/Ryan/GitHub/opossum-V1-omics/figure 2/data/';
output_dir = 'E:/Ryan/GitHub/opossum-V1-omics/figure 2/plots/';

% Consistent color definitions
colors.mouse_subclass = {'#FFB3B3'; '#FF7F50'; '#FFA07A'; '#FF6347'};  % L2/3, L4, L5IT, L6IT
colors.opossum_subclass = {'#FFB3B3'; '#FF7F50'; '#FF6347'; '#FFA07A'};  % IT_A, IT_B, IT_C, IT_D
colors.species = {'#AAAAAA'; '#C692B8'};  % Mouse, Opossum

% View angles (manually tuned)
views.mouse_3d = [168, 16];
views.mouse_flat = [8, 9];
views.opossum_3d = [444, 19];
views.opossum_flat = [-57, 12];
views.shared_subclass = [153, 16];

%% ========================================================================
%  PART 1: ParTI FITTING
%  ========================================================================

%% 2G: ParTI -> Mouse Species-Specific PC Space
input = readtable(fullfile(data_dir, '2G_Mouse_SubITPCEmbeddings.csv'));
data = [];
for pc = 1:10
    data = [data, input.(sprintf('PC_%i', pc))];
end
labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 10, labels, ...
    ident, 0, [], [], [], 0.2, 'data/Mouse_IT_10PCs');

save(fullfile(data_dir, 'Mouse_IT_10PCs_Data.mat'), 'arc', 'arcOrig', 'pc', 'errs', 'pval', 'coefs');

%% 2H: ParTI -> Opossum Species-Specific PC Space
input = readtable(fullfile(data_dir, '2H_Opossum_SubITPCEmbeddings.csv'));
data = [];
for pc = 1:10
    data = [data, input.(sprintf('PC_%i', pc))];
end
labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 10, labels, ...
    ident, 0, [], [], [], 0.2, 'data/Opossum_IT_10PCs');

save(fullfile(data_dir, 'Opossum_IT_10PCs_Data.mat'), 'arc', 'arcOrig', 'pc', 'errs', 'pval', 'coefs');

%% ========================================================================
%  PART 2: 3D VISUALIZATION
%  ========================================================================

%% Figure 2G: Mouse IT in Species-Specific Space
data = readtable(fullfile(data_dir, '2G_Mouse_SubProjITPCEmbeddings.csv'));

pc1 = data.pca_1 * -1;
pc2 = data.pca_2;
pc3 = data.pca_3 * -1;

[subclass_groups, ~, subclass_idx] = unique(data.subclass, 'stable');
cmap = hex2rgb(colors.mouse_subclass);

f = figure;
set(gcf, 'Renderer', 'Painters');
hold on;
for i = 1:length(subclass_groups)
    idx = (subclass_idx == i);
    scatter3(pc1(idx), pc2(idx), pc3(idx), 15, cmap(i, :), 'filled', 'DisplayName', subclass_groups{i});
end

% Load and plot tetrahedron
load(fullfile(data_dir, 'Mouse_IT_10PCs_Data.mat'), 'arcOrig');
plot_tetrahedron(arcOrig(:, 1:3) .* [-1, 1, -1]);

axis equal; axis square;
view(views.mouse_3d);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
set(gcf, 'Color', 'w');
grid on;
SavePNGandSVG(f, output_dir, '2G_Mouse-ITSubclass-Top3PCs');
hold off;

%% Figure 2O (Top): Mouse IT in Species-Specific Space (Flat View)
data = readtable(fullfile(data_dir, '2G_Mouse_SubProjITPCEmbeddings.csv'));

pc1 = data.pca_1 * -1;
pc2 = data.pca_2;
pc3 = data.pca_3 * -1;

[subclass_groups, ~, subclass_idx] = unique(data.subclass, 'stable');
cmap = hex2rgb(colors.mouse_subclass);

f = figure;
set(gcf, 'Renderer', 'Painters');
hold on;
for i = 1:length(subclass_groups)
    idx = (subclass_idx == i);
    scatter3(pc1(idx), pc2(idx), pc3(idx), 15, cmap(i, :), 'filled', 'DisplayName', subclass_groups{i});
end

view(views.mouse_flat);
set(gcf, 'Color', 'w');
axis square; axis off;
SavePNGandSVG(f, output_dir, '2O_Mouse-ITSubclass-Top3PCsFlat');
hold off;

%% Figure 2H: Opossum IT in Species-Specific Space
data = readtable(fullfile(data_dir, '2H_Opossum_SubProjITPCEmbeddings.csv'));

% COORDINATE TRANSFORMATION: Opossum requires axis swap for consistent orientation
% [PC1 * -1, PC3, PC2] maps to [x, y, z] for plotting
pc1 = data.pca_1 * -1;
pc2 = data.pca_3;
pc3 = data.pca_2 * -1;

[subclass_groups, ~, subclass_idx] = unique(data.subclass, 'stable');
cmap = hex2rgb(colors.opossum_subclass);

f = figure;
set(gcf, 'Renderer', 'Painters');
hold on;
for i = 1:length(subclass_groups)
    idx = (subclass_idx == i);
    scatter3(pc1(idx), pc2(idx), pc3(idx), 15, cmap(i, :), 'filled', 'DisplayName', subclass_groups{i});
end

% Load and plot tetrahedron with matching transformation
load(fullfile(data_dir, 'Opossum_IT_10PCs_Data.mat'), 'arcOrig');
plot_tetrahedron(arcOrig(:, 1:3) .* [-1, 1, -1]);

axis equal; axis square;
view(views.opossum_3d);
xlabel('PC1'); ylabel('PC3'); zlabel('PC2');
set(gcf, 'Color', 'w');
grid on;
SavePNGandSVG(f, output_dir, '2H_Opossum-ITSubclass-Top3PCs');
hold off;

%% Figure 2O (Bottom): Opossum IT in Species-Specific Space (Flat View)
data = readtable(fullfile(data_dir, '2H_Opossum_SubProjITPCEmbeddings.csv'));

pc1 = data.pca_1 * -1;
pc2 = data.pca_3;
pc3 = data.pca_2 * -1;

[subclass_groups, ~, subclass_idx] = unique(data.subclass, 'stable');
cmap = hex2rgb(colors.opossum_subclass);

f = figure;
set(gcf, 'Renderer', 'Painters');
hold on;
for i = 1:length(subclass_groups)
    idx = (subclass_idx == i);
    scatter3(pc1(idx), pc2(idx), pc3(idx), 15, cmap(i, :), 'filled', 'DisplayName', subclass_groups{i});
end

view(views.opossum_flat);
set(gcf, 'Color', 'w');
axis square; axis off;
SavePNGandSVG(f, output_dir, '2O_Opossum-ITSubclass-Top3PCsFlat');
hold off;

%% Figure 2H: Opossum in Shared Space (Mouse Gray, Opossum by Subclass)
data = readtable(fullfile(data_dir, '2I_OpossumMouse_MouseSubProjITPCEmbeddings.csv'));

% Replace Mouse subclass labels with 'Mouse' for uniform coloring
data(ismember(data.species, 'Mouse'), "subclass") = repmat({'Mouse'}, [nnz(ismember(data.species, 'Mouse')), 1]);

pc1 = data.pca_1 * -1;
pc2 = data.pca_2;
pc3 = data.pca_3 * -1;  % Flip for consistent orientation

[subclass_groups, ~, subclass_idx] = unique(data.subclass, 'stable');
% Colors: Mouse gray first, then opossum subclass colors
hex_colors = {'#AAAAAA'; '#FFB3B3'; '#FFA07A'; '#FF6347'; '#FF7F50'};
cmap = hex2rgb(hex_colors);

f = figure;
set(gcf, 'Renderer', 'Painters');
hold on;
% Plot in reverse order so opossum subclasses appear on top
for i = length(subclass_groups):-1:1
    idx = (subclass_idx == i);
    scatter3(pc1(idx), pc2(idx), pc3(idx), 5, cmap(i, :), 'filled', 'DisplayName', subclass_groups{i});
end

% Load and plot tetrahedron
load(fullfile(data_dir, 'Mouse_IT_10PCs_Data.mat'), 'arcOrig');
plot_tetrahedron(arcOrig(:, 1:3) .* [-1, 1, -1]);

axis equal; axis square;
view(views.shared_subclass);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
set(gcf, 'Color', 'w');
grid on;
SavePNGandSVG(f, output_dir, '2I_OpossumMouse-ITSubclass-Top3PCs');
hold off;

%% ========================================================================
%  PART 3: SILHOUETTE ANALYSIS
%  ========================================================================

%% 2K. Cross-Species Silhouette Scores
data = readtable(fullfile(data_dir, '2I_OpossumMouse_MouseSubProjITPCEmbeddings.csv'));
mouse_data = data(ismember(data.species, 'Mouse'), :);
opossum_data = data(ismember(data.species, 'Opossum'), :);

mouse_pc_matrix = [mouse_data.pca_1, mouse_data.pca_2, mouse_data.pca_3];
opossum_pc_matrix = [opossum_data.pca_1, opossum_data.pca_2, opossum_data.pca_3];

combined_pc_matrix = [mouse_pc_matrix; opossum_pc_matrix];
combined_labels = [mouse_data.subclass; opossum_data.subclass];

mouse_subclasses = unique(mouse_data.subclass);
opossum_subclasses = unique(opossum_data.subclass);

cross_species_silhouette = compute_cross_species_silhouette(...
    combined_pc_matrix, combined_labels, mouse_subclasses, opossum_subclasses);

f = figure;
set(f, 'Color', 'w', 'Position', [100, 100, 250, 200]);
% Define white-red-blue colormap
custom_colormap = [linspace(1,1,256)', linspace(1,0,256)', linspace(1,0,256)';  % White to red
                   linspace(1,0,256)', zeros(256,1), linspace(0,1,256)'];       % Red to blue
plot_cross_species_heatmap(cross_species_silhouette, mouse_subclasses, opossum_subclasses, custom_colormap);
SavePNGandSVG(f, output_dir, '2K_OpossumMouse-ITSubclassSilhouetteCross-Heatmap');

%% Figure 2P. Within-Species Silhouette Scores
mouse_data = readtable(fullfile(data_dir, '2G_Mouse_SubProjITPCEmbeddings.csv'));
opossum_data = readtable(fullfile(data_dir, '2H_Opossum_SubProjITPCEmbeddings.csv'));

mouse_pc_matrix = [mouse_data.pca_1, mouse_data.pca_2, mouse_data.pca_3];
opossum_pc_matrix = [opossum_data.pca_1, opossum_data.pca_2, opossum_data.pca_3];

mouse_subclasses = unique(mouse_data.subclass);
opossum_subclasses = unique(opossum_data.subclass);

% Compute pairwise silhouette for mouse
mouse_silhouette_matrix = compute_pairwise_silhouette(mouse_pc_matrix, mouse_data.subclass, mouse_subclasses);

% Compute pairwise silhouette for opossum
opossum_silhouette_matrix = compute_pairwise_silhouette(opossum_pc_matrix, opossum_data.subclass, opossum_subclasses);

% Plot heatmaps
f = figure;
set(f, 'Color', 'w', 'Position', [100, 100, 600, 200]);
custom_colormap = [linspace(1,0,256)', zeros(256,1), linspace(0,1,256)'];
subplot(1,2,1);
plot_silhouette_heatmap(mouse_silhouette_matrix, mouse_subclasses, 'Mouse IT Silhouette', custom_colormap);
subplot(1,2,2);
plot_silhouette_heatmap(opossum_silhouette_matrix, opossum_subclasses, 'Opossum IT Silhouette', custom_colormap);
% Save
SavePNGandSVG(f, output_dir, '2P_OpossumMouse-ITSubclassSilhouetteWithin-Heatmap');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function cmap = hex2rgb(hex_colors)
    % Convert cell array of hex colors to RGB matrix
    cmap = cell2mat(cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x')'/255, hex_colors, 'UniformOutput', false));
end

function plot_tetrahedron(vertices)
    % Plot tetrahedron edges given 4x3 vertex matrix
    tetra_edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    for k = 1:size(tetra_edges, 1)
        plot3(vertices(tetra_edges(k, :), 1), ...
              vertices(tetra_edges(k, :), 2), ...
              vertices(tetra_edges(k, :), 3), '-k', 'LineWidth', 1.5);
    end
end

function sil_matrix = compute_pairwise_silhouette(pc_matrix, labels, subclasses)
    % Compute silhouette scores for all pairs of subclasses
    n = length(subclasses);
    sil_matrix = zeros(n);
    
    for i = 1:n
        for j = 1:n
            if i ~= j
                idx = strcmp(labels, subclasses{i}) | strcmp(labels, subclasses{j});
                labels_subset = labels(idx);
                pc_subset = pc_matrix(idx, :);
                numeric_labels = grp2idx(categorical(labels_subset));
                scores = silhouette(pc_subset, numeric_labels);
                sil_matrix(i, j) = median(scores);
            end
        end
    end
end

function cross_sil = compute_cross_species_silhouette(pc_matrix, labels, subclasses1, subclasses2)
    % Compute cross-species silhouette matrix
    cross_sil = zeros(length(subclasses1), length(subclasses2));
    
    for i = 1:length(subclasses1)
        for j = 1:length(subclasses2)
            idx = strcmp(labels, subclasses1{i}) | strcmp(labels, subclasses2{j});
            labels_subset = labels(idx);
            pc_subset = pc_matrix(idx, :);
            numeric_labels = grp2idx(categorical(labels_subset));
            scores = silhouette(pc_subset, numeric_labels);
            cross_sil(i, j) = median(scores);
        end
    end
end

function plot_silhouette_heatmap(sil_matrix, subclasses, title_str, custom_colormap)
    % Plot silhouette heatmap
    h = heatmap(subclasses, subclasses, round(sil_matrix, 2), ...
        'Colormap', custom_colormap, 'ColorLimits', [0.75, 1]);
    title(title_str);
    h.CellLabelFormat = '%.2f';
    h.GridVisible = 'off';
end

function plot_cross_species_heatmap(cross_sil, subclasses1, subclasses2, custom_colormap)
    % Plot cross-species silhouette heatmap
    h = heatmap(subclasses2, subclasses1, round(cross_sil, 2)', ...
        'Colormap', custom_colormap, 'ColorLimits', [0, 1]);
    title('Cross-Species IT Silhouette');
    xlabel('Opossum Subclass');
    ylabel('Mouse Subclass');
    h.CellLabelFormat = '%.2f';
    h.GridVisible = 'off';
end

function SavePNGandSVG(fig, out_dir, base_name)
% SavePNGandSVG Save a figure as PNG + SVG (with fallbacks).
%
% Usage:
%   SavePNGandSVG(fig, out_dir, base_name)
%   SavePNGandSVG(gcf, output_dir, 'my_figure')

    png_path = fullfile(out_dir, [base_name '.png']);
    svg_path = fullfile(out_dir, [base_name '.svg']);

    % PNG
    if exist('exportgraphics','file') == 2
        exportgraphics(fig, png_path, 'Resolution', 300);
    else
        print(fig, png_path, '-dpng', '-r300');
    end

    % SVG
    set(fig, 'Renderer', 'painters');
    try
        print(fig, svg_path, '-dsvg');
    catch
        eps_path = fullfile(out_dir, [base_name '.eps']);
        print(fig, eps_path, '-depsc');
        warning('SVG export unsupported in this MATLAB. Saved EPS instead: %s', eps_path);
    end
end
