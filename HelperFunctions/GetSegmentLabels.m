function [grid_tags, palm_seg_present, dors_seg_present] = GetSegmentLabels(SurveyStruct, ChannelMap, pbt_regions, dbt_regions, method, orig_size)
    if nargin == 5
        orig_size = [1200, 1050];
    end
    % Tags for palmar and dorsum
    grid_locations = ChannelMap.ChannelNumbers(ChannelMap.IsSensory);
    [grid_tags{1, 1}, grid_tags{1, 2}, grid_tags{2, 1}, grid_tags{2, 2}] = deal(cell(size(grid_locations{1})));
    % To keep track of which hand segments to highlight
    palm_seg_present = false(size(pbt_regions));
    dors_seg_present = false(size(dbt_regions));

    if contains(method, 'Weighted', 'IgnoreCase', true) && ~isfield(SurveyStruct,'ProcessedSFs')
        error('When using Weighted* methods, the SurveyStruct must contain the ProcessedSFs field')
    end
    
    for i = 1:length(SurveyStruct)
        n_sessions = length(SurveyStruct(i).SensoryField);
        % Workout where to put tag
        if SurveyStruct(i).Channel <= 32
            [row_idx,col_idx] = find(grid_locations{1} == SurveyStruct(i).Channel);
            array_idx = 1;
        else
            [row_idx,col_idx] = find(grid_locations{2} == SurveyStruct(i).Channel);
            array_idx = 2;
        end

        if strcmpi(method, 'DailyCentroid') % Plot by median centroid across days
                % Loop through sessions and get each centroid
                temp_palm_centroid = [];
                temp_dorsum_centroid = [];
                for s = 1:n_sessions
                    % Palm loop
                    n_p_rfs = size(SurveyStruct(i).SensoryField(s).PalmSF);
                    for r = 1:n_p_rfs % Number rfs
                        temp_palm_centroid = cat(1, temp_palm_centroid, SurveyStruct(i).SensoryField(s).PalmSF(r).Centroid);
                    end
            
                    % Dorsum loop
                    n_d_rfs = size(SurveyStruct(i).SensoryField(s).DorsSF);
                    for r = 1:n_d_rfs % Number rfs
                        temp_dorsum_centroid = cat(1, temp_dorsum_centroid, SurveyStruct(i).SensoryField(s).DorsSF(r).Centroid);
                    end
                end
    
                % Across-session centroids - Palm
                if ~isempty(temp_palm_centroid)
                    palm_centroid = mean(temp_palm_centroid,1);
                    segment_idx = GetSegment_Centroid(palm_centroid, orig_size, pbt_regions);
                    palm_seg_present(segment_idx) = true;
                    grid_tags{1,array_idx}{row_idx,col_idx} = pbt_regions(segment_idx).Tag;
                end
                
                % Across-session centroids - dorsum
                if ~isempty(temp_dorsum_centroid)
                    dorsum_centroid = mean(temp_dorsum_centroid,1);
                    segment_idx = GetSegment_Centroid(dorsum_centroid, orig_size, dbt_regions);
                    dors_seg_present(segment_idx) = true;
                    grid_tags{2,array_idx}{row_idx,col_idx} = dbt_regions(segment_idx).Tag;
                end
                
    
        elseif strcmpi(method, 'DailyOverlap') % Plot by the overlap across days    
                % Loop through sessions and check individual overlap
                temp_palm_idx = [];
                temp_dors_idx = [];
                for s = 1:n_sessions
                    % Palm loop
                    n_p_rfs = size(SurveyStruct(i).SensoryField(s).PalmSF);
                    for r = 1:n_p_rfs % Number rfs
                        temp_palm_idx = cat(1, temp_palm_idx, SurveyStruct(i).SensoryField(s).PalmSF(r).PixelIdxList);
                    end
            
                    % Dorsum loop
                    n_d_rfs = size(SurveyStruct(i).SensoryField(s).DorsSF);
                    for r = 1:n_d_rfs % Number rfs
                        temp_dors_idx = cat(1, temp_dors_idx, SurveyStruct(i).SensoryField(s).DorsSF(r).PixelIdxList);
                    end
                end
                
                % Assign palms
                palm_seg_list = GetSegment_PixelIdx(temp_palm_idx, pbt_regions);
                if ~isempty(palm_seg_list)
                    palm_seg_present(palm_seg_list) = true;
                    grid_tags{1,array_idx}{row_idx,col_idx} = {pbt_regions(palm_seg_list).Tag};
                end
    
                % Assign dorsum
                dors_seg_list = GetSegment_PixelIdx(temp_dors_idx, dbt_regions);
                if ~isempty(dors_seg_list)
                    dors_seg_present(dors_seg_list) = true;
                    grid_tags{2,array_idx}{row_idx,col_idx} = {dbt_regions(dors_seg_list).Tag};
                end
    
        elseif strcmpi(method, 'WeightedCentroid') % Plot by the weighted centroid
            palm_centroids = cat(1,SurveyStruct(i).ProcessedSFs.PalmSF.WeightedCentroid{:});
            for p = 1:size(palm_centroids,1) % Find relevant regions
                segment_idx = GetSegment_Centroid(palm_centroids(p,:), orig_size, pbt_regions);
                palm_seg_present(segment_idx) = true;
                grid_tags{1,array_idx}{row_idx,col_idx}{p} = pbt_regions(segment_idx).Tag;
            end

            dorsum_centroids = cat(1,SurveyStruct(i).ProcessedSFs.DorsSF.WeightedCentroid{:});
            for d = 1:size(dorsum_centroids,1) % Find relevant regions
                segment_idx = GetSegment_Centroid(dorsum_centroids(d,:), orig_size, dbt_regions);
                dors_seg_present(segment_idx) = true;
                grid_tags{2,array_idx}{row_idx,col_idx}{d} = dbt_regions(segment_idx).Tag;
            end
    
        elseif strcmpi(method, 'WeightedOverlap') % Plot by the weighted overlap
            temp_palm_idx = cat(1,SurveyStruct(i).ProcessedSFs.PalmSF.PixelIds{:});
            % Assign palms
            palm_seg_list = GetSegment_PixelIdx(temp_palm_idx, pbt_regions);
            if ~isempty(palm_seg_list)
                palm_seg_present(palm_seg_list) = true;
                grid_tags{1,array_idx}{row_idx,col_idx} = {pbt_regions(palm_seg_list).Tag};
            end
            
            temp_dors_idx = cat(1,SurveyStruct(i).ProcessedSFs.DorsSF.PixelIds{:});
            % Assign dorsum
            dors_seg_list = GetSegment_PixelIdx(temp_dors_idx, dbt_regions);
            if ~isempty(dors_seg_list)
                dors_seg_present(dors_seg_list) = true;
                grid_tags{2,array_idx}{row_idx,col_idx} = {dbt_regions(dors_seg_list).Tag};
            end
    
        else
            error('Invalid method. Must be: DailyCentroid, DailyOverlap, WeightedCentroid, or WeightedOverlap')
        end

    end
end