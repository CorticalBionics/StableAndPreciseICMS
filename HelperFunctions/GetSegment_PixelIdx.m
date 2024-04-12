function segment_list = GetSegment_PixelIdx(PixelIdxList, segmenter)
    segment_list = [];
    for segment_idx = 1:size(segmenter,1)
        if any(ismember(PixelIdxList, segmenter(segment_idx).PixelIdxList))
            segment_list = cat(1, segment_list, segment_idx);
        end
    end
end