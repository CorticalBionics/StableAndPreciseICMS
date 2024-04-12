function segment_idx = GetSegment_Centroid(centroid, orig_size, segmenter)
    centroid = round(centroid); % Ensure is idx
    centroid_idx = sub2ind(orig_size, centroid(2), centroid(1)); % xy is reversed
    % Segment - within segment
    for segment_idx = 1:size(segmenter,1)
        if ismember(centroid_idx, segmenter(segment_idx).PixelIdxList)
            return
        end
    end
    % Segment - nearest segment - won't get here if already returned
    [~, segment_idx] = min(sum(abs(centroid - cat(1,segmenter.Centroid)),2));
end

