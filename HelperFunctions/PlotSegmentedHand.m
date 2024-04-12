function f = PlotSegmentedHand(ChannelMap, grid_tags, pbt_regions, palm_seg_idx, dbt_regions, dors_seg_idx, text_size) 
    if nargin < 7
        text_size = 0;
    end
    % Load and scale ref images
    [a,~,~] = fileparts(mfilename('fullpath'));
    orig_size = [1200, 1050];
    % Palmar
    [palm_ref_img, ~, palm_ref_alpha] = imread(fullfile(a, '..', 'ReferenceImages', 'TopLayer-handpcontour.png'));
    palm_ref_img = imresize(palm_ref_img, orig_size);
    palm_ref_alpha = imresize(palm_ref_alpha, orig_size);
    % Dorsum
    [dors_ref_img, ~, dors_ref_alpha] = imread(fullfile(a, '..', 'ReferenceImages', 'TopLayer-contour.png'));
    dors_ref_img = imresize(dors_ref_img, orig_size);
    dors_ref_alpha = imresize(dors_ref_alpha, orig_size);

    kw = 4;
    kernel = ones(kw) / kw ^ 2;

    sensory_names = ChannelMap.ArrayNames(ChannelMap.IsSensory);
    array_rotations = ChannelMap.ArrayRotations(ChannelMap.IsSensory);
    grid_locations = ChannelMap.ChannelNumbers(ChannelMap.IsSensory);
    if contains(sensory_names{1}, 'Medial')
        a1_h = .55;
        a2_h = .075;
    elseif contains(sensory_names{1}, 'Lateral')
        a2_h = .55;
        a1_h = .075;
    end

    inset_lims = sqrt(sum(size(grid_locations{1}).^2)) / 2;
    
    [palm_colored, dorsum_colored] = deal(ones(size(palm_ref_img)));
    f = figure; 
    ax(1) = axes('Position', [.0 0.05 .35 .9]); hold on
        % Set background image colors
        c_idx = find(palm_seg_idx);
        for c = 1:length(c_idx)
            palm_colored(pbt_regions(c_idx(c)).PixelIdxList) = pbt_regions(c_idx(c)).Color(1);
            palm_colored(pbt_regions(c_idx(c)).PixelIdxList + prod(orig_size)) = pbt_regions(c_idx(c)).Color(2);
            palm_colored(pbt_regions(c_idx(c)).PixelIdxList + 2 * prod(orig_size)) = pbt_regions(c_idx(c)).Color(3);
        end
        % Blur the image to de-emphasize segment boundaries
        palm_colored = imfilter(palm_colored, kernel);
        imshow(palm_colored)
        % Overlay ref img
        image(palm_ref_img, 'AlphaData', palm_ref_alpha);
        set(ax(1), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', [80, 1000], 'YLim', [5 1150])
    
    
    inset(1) = axes('Position', [.275 a1_h .175 .35]); hold on
        PlotGridSegments(grid_locations{1}, grid_tags{1,1}, pbt_regions, text_size, inset(1), array_rotations(1))    
        set(inset(1), 'DataAspectRatio', [1 1 1], ...
                      'XTick', [], 'XColor', 'none', 'XLim', [-inset_lims, inset_lims], ...
                      'YTick', [], 'YColor', 'none', 'YLim', [-inset_lims, inset_lims], ...
                      'XDir', 'normal', 'YDir', 'reverse', 'Color', 'none')
    
    inset(2) = axes('Position', [.275 a2_h .175 .35]); hold on
        PlotGridSegments(grid_locations{2}, grid_tags{1,2}, pbt_regions, text_size, inset(2), array_rotations(2))    
        set(inset(2), 'DataAspectRatio', [1 1 1], ...
                      'XTick', [], 'XColor', 'none', 'XLim', [-inset_lims, inset_lims], ...
                      'YTick', [], 'YColor', 'none', 'YLim', [-inset_lims, inset_lims], ...
                      'XDir', 'normal', 'YDir', 'reverse', 'Color', 'none')
    
    
    ax(2) = axes('Position', [.67 0.05 .35 .9]); hold on
        % Set background image colors
        c_idx = find(dors_seg_idx);
        for c = 1:length(c_idx)
            dorsum_colored(dbt_regions(c_idx(c)).PixelIdxList) = dbt_regions(c_idx(c)).Color(1);
            dorsum_colored(dbt_regions(c_idx(c)).PixelIdxList + prod(orig_size)) = dbt_regions(c_idx(c)).Color(2);
            dorsum_colored(dbt_regions(c_idx(c)).PixelIdxList + 2 * prod(orig_size)) = dbt_regions(c_idx(c)).Color(3);
        end
        % Blur the image to de-emphasize segment boundaries
        dorsum_colored = imfilter(dorsum_colored, kernel);
        imshow(dorsum_colored)
        % Overlay ref img
        image(dors_ref_img, 'AlphaData', dors_ref_alpha);
        set(ax(2), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', [80, 1000], 'YLim', [5 1150])
    
    
    inset(3) = axes('Position', [.525 a1_h .175 .35]); hold on
        PlotGridSegments(grid_locations{1}, grid_tags{2,1}, dbt_regions, text_size, inset(3), array_rotations(1))
        set(inset(3), 'DataAspectRatio', [1 1 1], ...
                      'XTick', [], 'XColor', 'none', 'XLim', [-inset_lims, inset_lims], ...
                      'YTick', [], 'YColor', 'none', 'YLim', [-inset_lims, inset_lims], ...
                      'XDir', 'normal', 'YDir', 'reverse', 'Color', 'none')
    
    
    inset(4) = axes('Position', [.525 a2_h .175 .35]); hold on
        PlotGridSegments(grid_locations{2}, grid_tags{2,2}, dbt_regions, text_size, inset(4), array_rotations(2))
        set(inset(4), 'DataAspectRatio', [1 1 1], ...
                      'XTick', [], 'XColor', 'none', 'XLim', [-inset_lims, inset_lims], ...
                      'YTick', [], 'YColor', 'none', 'YLim', [-inset_lims, inset_lims], ...
                      'XDir', 'normal', 'YDir', 'reverse', 'Color', 'none')
    
    set(gcf, 'Position', [50, 50, 1500, 700])
end