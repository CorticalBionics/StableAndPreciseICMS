function pf_table = CombineSensoryFields(pf_data, segmenter, orig_size, threshold)
    % Setup table
    variable_names_types = [["SegmentName", "string"]; ...
                            ["PixelIdxList", "cell"]; ...
                            ["PixelIds", "cell"]; ...
                            ["Count", "cell"]; ...
                            ["NormalizedCount", "cell"];...
                            ["Area", "double"];...
                            ["WeightedCentroid", "cell"];...
                            ["Present", "logical"];...
                            ["TotalWeight", "double"];...
                            ["Dominant", "logical"];...
                            ["GlobalWeightRatio", "double"]];
    
    pf_table = table('Size', [length(segmenter), size(variable_names_types, 1)],...
                     'VariableNames', variable_names_types(:,1),...
                     'VariableTypes', variable_names_types(:,2));
                
    if isempty(pf_data) % Skip completely empty
        return
    end
    
    % Remove missing data
    pf_data = pf_data(~cellfun(@(c) isempty(c), {pf_data(:).Area}));
    
    centroid_loc = {pf_data(:).CentroidLocation};
    pixel_idx_list = {pf_data(:).PixelIdxList};
    [pixel_counts, ~] = groupcounts(vertcat(pixel_idx_list{:})); % Count unique values
    global_pixels = sum(pixel_counts);
    total_pixels = sum(pixel_counts(pixel_counts >= threshold));

    % Threshold cut off calculations
    for sl = 1:length(segmenter)
        pf_table.SegmentName(sl) = segmenter{sl};
        seg_idx = find(contains(centroid_loc, segmenter{sl})); % Find matching names with segment
        if ~isempty(seg_idx)
            tmp_pixel_list = vertcat(pixel_idx_list{seg_idx});  % Concatenate all PixelIdxs
            [pixel_counts, pixel_ids] = groupcounts(tmp_pixel_list); % Count unique values
            above_thresh_idx = pixel_counts >= threshold; % Determine which pixels cross threshold
            pixel_counts = pixel_counts(above_thresh_idx);
            pixel_ids = pixel_ids(above_thresh_idx);
            normalized_counts = pixel_counts ./ total_pixels;
    
            if any(above_thresh_idx) % If there are still pixels present after thresholding
                pf_table.PixelIdxList{sl} = tmp_pixel_list(ismember(tmp_pixel_list, pixel_ids));
                pf_table.Count{sl} = pixel_counts;
                pf_table.PixelIds{sl} = pixel_ids;
                pf_table.NormalizedCount{sl} = normalized_counts;
                pf_table.Area(sl) = length(pixel_ids) / 25; % to convert to mm^2
                if pf_table.Area(sl) > 0
                    pf_table.Present(sl) = true;
                    pf_table.TotalWeight(sl)= sum(pf_table.NormalizedCount{sl});
                    pf_table.GlobalWeightRatio(sl) = sum(pf_table.Count{sl}) / global_pixels;
                end
    
                % Centroid Calculation
                [pixel_y, pixel_x] = ind2sub(orig_size, pf_table.PixelIds{sl});
                centroid_x = sum(pf_table.Count{sl} .* pixel_x) ./ sum(pf_table.Count{sl});
                centroid_y = sum(pf_table.Count{sl} .* pixel_y) ./ sum(pf_table.Count{sl});
                pf_table.WeightedCentroid{sl} = [centroid_x, centroid_y];
            end
        end
    end
    if any(pf_table.TotalWeight)
        [~,max_idx] = max(pf_table.TotalWeight);
        pf_table.Dominant(max_idx) = true;
    end
end
