function [pbt_regions, dbt_regions] = GetHandSegments(plot_segments)
    if nargin == 0
        plot_segments = false;
    end
    orig_size = [1200, 1050];
    % Palm segments
    [a,b,~] = fileparts(mfilename('fullpath'));
    palm_boundary_template = imread(fullfile(a, '..', 'ReferenceImages','BoundaryMapPalm.png'));
    pbt_resize = imresize(palm_boundary_template, orig_size); % Resize to survey hand map image
    pbt_flat = ~logical(mean(pbt_resize,3));
    pbt_thin = ~bwmorph(pbt_flat, 'skel', Inf);
    pbt_thin = ~bwmorph(~pbt_thin, 'spur', Inf);
    pbt_connections = bwconncomp(pbt_thin, 4);
    pbt_regions = regionprops(pbt_connections, 'Centroid', 'PixelIdxList');
    pbt_regions = pbt_regions(2:end); % Skip outline

    % Manually add regions names and colors - yes this was tedious
     region_tag_col = { ...
        'D5d-du', [30, 69, 31]; ...
        'D5d-dr', [30, 69, 31]; ...
        'D5d-pu', [40 , 73, 41]; ...
        'D5m-u', [51, 78 52]; ...
        'D5d-pr', [40 , 73, 41]; ...
        'D5m-r', [51, 78 52]; ...
        'D5p-u', [65, 84, 66]; ...
        'D5p-r', [65, 84, 66]; ...
        'P5-mcp', [78, 90, 79]; ...
        'P0-du', [70, 87, 86]; ...
        'P0-pu', [88, 97, 98]; ...
        'D4d-du', [1, 66, 96]; ...
        'D4d-pu', [16, 71, 97]; ...
        'D4m-u', [31, 77, 97]; ...
        'W',    [94, 92, 91]; ...
        'D4p-u', [51, 83, 98]; ...
        'D4d-dr', [1, 66, 96]; ...
        'D4d-pr', [16, 71, 97]; ...
        'P4-mcp', [70, 90, 99]; ...
        'D4m-r', [31, 77, 97]; ...
        'D4p-r', [51, 83, 98]; ...
        'P0-d', [88, 75, 91]; ...
        'P3-mcp', [82, 77, 91]; ...
        'D3p-u', [70, 62, 86]; ...
        'D3d-du', [40, 23, 72]; ...
        'D3d-pu', [50, 34, 76]; ...
        'D3m-u', [58, 46, 80]; ...
        'D3d-dr', [40, 23, 72]; ...
        'D3d-pr', [50, 34, 76]; ...
        'D3m-r', [58, 46, 80]; ...
        'P0-dr', [97, 73, 82]; ...
        'P0-pr', [100, 99, 90]; ...
        'D3p-r', [70, 62, 86]; ...
        'P2-mcp', [100, 80, 82]; ...
        'D2p-u', [94, 60, 60]; ...
        'D2m-u', [90, 45, 45]; ...
        'D2p-r', [94, 60, 60]; ...
        'D2d-pu', [94, 33, 31]; ...
        'D2m-r', [90, 45, 45]; ...
        'D2d-du', [96, 26, 21]; ...
        'D2d-pr', [94, 33, 31]; ...
        'D2d-dr', [96, 26, 21]; ...
        'D1p-u', [100, 80, 50]; ...
        'D1p-r', [100, 80, 50]; ...
        'D1d-pu', [100, 72, 30]; ...
        'D1d-pr', [100, 72, 30]; ...
        'D1d-du', [100, 66, 15]; ...
        'D1d-dr', [100, 66, 15]};   

    for i = 1:size(pbt_regions,1)
        pbt_regions(i).Tag = region_tag_col{i,1};
        pbt_regions(i).Color = region_tag_col{i,2} ./ 100;
    end
    
    % Dorsum segments
    dorsum_boundary_template = imread(fullfile(a, '..', 'ReferenceImages','BoundaryMapDorsum.png'));
    dbt_resize = imresize(dorsum_boundary_template, orig_size); % Resize to survey hand map image
    dbt_flat = ~logical(mean(dbt_resize,3));
    dbt_thin = ~bwmorph(dbt_flat, 'skel', Inf);
    dbt_thin = ~bwmorph(~dbt_thin, 'spur', Inf);
    dbt_connections = bwconncomp(dbt_thin, 4);
    dbt_regions = regionprops(dbt_connections,'Centroid', 'PixelIdxList');
    dbt_regions = dbt_regions(2:end); % Skip outline

    % Manually add regions names
    region_tag_col = {
        'D1d-r', [100, 66, 15]; ....
        'D1d-u', [100, 66, 15]; ....
        'D1p', [100, 72, 30]; ....
        'P-pr', [100, 99, 91]; ....
        'D2d-r', [94, 33, 31]; ....
        'D2m', [90, 45, 45]; ....
        'D2p', [94, 60, 60]; ....
        'D2d-u', [94, 33, 31]; ....
        'P-dr', [97, 73, 82]; ....
        'W', [94, 92, 91]; ....
        'D3d-r', [50, 34, 76]; ....
        'D3m', [58, 46, 80]; ....
        'D3p', [70, 62, 86]; ....
        'D3d-u', [50, 34, 76]; ....
        'P-pu', [88, 97, 98]; ....
        'P-du', [70, 87, 86]; ....
        'D4p', [51, 83, 98]; ....
        'D4m', [31, 77, 97]; ....
        'D4d-r', [16, 71, 97]; ....
        'D4d-u', [16, 71, 97]; ....
        'D5p', [65, 84, 66]; ....
        'D5m', [51, 78, 52]; ....
        'D5d-r', [40, 73, 42]; ....
        'D5d-u', [40, 73, 42]};

    for i = 1:size(dbt_regions,1)
        dbt_regions(i).Tag = region_tag_col{i,1};
        dbt_regions(i).Color = region_tag_col{i,2} ./ 100;
    end
    
    if plot_segments
        figure;
        ax(1) = axes('Position', [.05 0.05 .35 .9]); hold on
            % Convert logical to image
            pbt_segment_image = repmat(double(pbt_thin), [1,1,3]);
            % Assign segment colors
            for i = 1:size(pbt_regions,1)
                  pbt_segment_image(pbt_regions(i).PixelIdxList) = pbt_regions(i).Color(1);
                  pbt_segment_image(pbt_regions(i).PixelIdxList + prod(orig_size)) = pbt_regions(i).Color(2);
                  pbt_segment_image(pbt_regions(i).PixelIdxList + 2 * prod(orig_size)) = pbt_regions(i).Color(3);
            end
            imshow(pbt_segment_image);

            % Text overlay
            for i = 1:size(pbt_regions,1)
                text(pbt_regions(i).Centroid(1), pbt_regions(i).Centroid(2), pbt_regions(i).Tag,...
                    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 9)
            end
                
            set(ax(1), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', [80, 1000], 'YLim', [0 1200], 'Box', 'off')
        
        ax(2) = axes('Position', [.5 0.05 .35 .9]); hold on
            % Convert logical to image
            dbt_segment_image = repmat(double(dbt_thin), [1,1,3]);
            % Assign segment colors
            for i = 1:size(dbt_regions,1)
                  dbt_segment_image(dbt_regions(i).PixelIdxList) = dbt_regions(i).Color(1);
                  dbt_segment_image(dbt_regions(i).PixelIdxList + prod(orig_size)) = dbt_regions(i).Color(2);
                  dbt_segment_image(dbt_regions(i).PixelIdxList + 2 * prod(orig_size)) = dbt_regions(i).Color(3);
            end
            imshow(dbt_segment_image);

            % Text overlay
            for i = 1:size(dbt_regions,1)
                text(dbt_regions(i).Centroid(1), dbt_regions(i).Centroid(2), dbt_regions(i).Tag,...
                    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 9)
            end
            set(ax(2), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', [80, 1000], 'YLim', [0 1150], 'Box', 'off')

    end
end