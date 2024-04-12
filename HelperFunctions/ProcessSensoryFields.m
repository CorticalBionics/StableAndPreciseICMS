function InputStruct = ProcessSensoryFields(InputStruct, SegmentationStrings, varargin)
    % Defaults
    Threshold = 0;
    ExportPlots = false;
    ExportFolder = "";
    ShowCentroids = false;
    ParseVarargin()

    % Create output folder
    if ExportPlots
        if ~exist(ExportFolder, 'dir')
            mkdir(ExportFolder)
        end
    end
    
    % Things that should not be touched
    orig_size = [1200, 1050];
    
    if ExportPlots
        [palmar_mask, palmar_template, dorsum_mask, dorsum_template] = GetHandMasks();
        fig = figure('Name', 'Processed Sensory Fields', 'Position', [50, 50, 1500, 700]);
        cmap = [linspace(1,.83,255); linspace(.92,.18,255); linspace(.93,.18,255)]';

        ax(1) = subplot(1,2,1); hold on
            % Add template
            imshow(palmar_template, 'Parent',ax(1))
            % Prepare empty handles
            palmar_image = zeros(orig_size);
            i1 = imagesc(palmar_image, 'AlphaData', palmar_image ~= 0, 'Parent', ax(1));
            colormap(ax(1), cmap)
            s1 = scatter([],[], 50, 'kx', 'Parent', ax(1), 'LineWidth', 1.5);
            set(ax(1), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', [80, 1000], 'YLim', [0 1150])

        ax(2) = subplot(1,2,2); hold on
            % Add template
            imshow(dorsum_template, 'Parent',ax(2))
            % Prepare empty handles
            dorsum_image = zeros(orig_size);
            i2 = imagesc(dorsum_image, 'AlphaData', dorsum_image ~= 0, 'Parent', ax(2));
            colormap(ax(2), cmap)
            s2 = scatter([],[], 50, 'kx', 'Parent', ax(2), 'LineWidth', 1.5);
            set(ax(2), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', [80, 1000], 'YLim', [0 1150])
    end
    
    for ch = 1:size(InputStruct, 2)
        ch_threshold = Threshold * size(InputStruct(ch).SensoryField,2);
        InputStruct(ch).ProcessedSFs.RawThreshold = Threshold;
        InputStruct(ch).ProcessedSFs.AdjustedThreshold = ch_threshold;
        InputStruct(ch).ProcessedSFs.PalmSF = CombineSensoryFields(vertcat(InputStruct(ch).SensoryField(:).PalmSF), ...
            SegmentationStrings, orig_size, ch_threshold);
        InputStruct(ch).ProcessedSFs.DorsSF = CombineSensoryFields(vertcat(InputStruct(ch).SensoryField(:).DorsSF), ...
            SegmentationStrings, orig_size, ch_threshold);
    
        if ExportPlots
            % Figure plotting
            palmar_image = zeros(orig_size);
            dorsum_image = zeros(orig_size);
            palmar_centroid = NaN(length(SegmentationStrings), 2);
            dorsum_centroid = NaN(length(SegmentationStrings), 2);
    
            for sl = 1:length(SegmentationStrings)
                if ~isempty(InputStruct(ch).ProcessedSFs.PalmSF.PixelIds{sl})
                    palmar_image(InputStruct(ch).ProcessedSFs.PalmSF.PixelIds{sl}) = InputStruct(ch).ProcessedSFs.PalmSF.Count{sl};
                    palmar_centroid(sl, :) = InputStruct(ch).ProcessedSFs.PalmSF.WeightedCentroid{sl};
                end
                if ~isempty(InputStruct(ch).ProcessedSFs.DorsSF.PixelIds{sl})
                    dorsum_image(InputStruct(ch).ProcessedSFs.DorsSF.PixelIds{sl}) = InputStruct(ch).ProcessedSFs.DorsSF.Count{sl};
                    dorsum_centroid(sl, :) = InputStruct(ch).ProcessedSFs.DorsSF.WeightedCentroid{sl};
                end
            end
            % Rescale 
            palmar_image = palmar_image ./ size(InputStruct(ch).SensoryField,2) .* 255;
            dorsum_image = dorsum_image ./ size(InputStruct(ch).SensoryField,2) .* 255;

            % Format
            if length(InputStruct(ch).Channel) > 1
                fname = sprintf(['Subject%s_Channel%d', repmat('.%d', [1, length(InputStruct(ch).Channel)-1])],...
                    InputStruct(ch).Subject, InputStruct(ch).Channel);
            else
                fname = sprintf('Subject%s_Channel%d', InputStruct(ch).Subject, InputStruct(ch).Channel);
            end
            title(ax(1),strrep(fname, '_', ' '), 'FontWeight', 'bold')
            set(i1, 'CData', palmar_image, 'AlphaData', (palmar_image ~= 0) .* 0.8, 'CDataMapping', 'direct')
            set(i2, 'CData', dorsum_image, 'AlphaData', (dorsum_image ~= 0) .* 0.8, 'CDataMapping', 'direct')

            if ShowCentroids
                set(s1, 'XData', palmar_centroid(:,1), 'YData', palmar_centroid(:,2))
                set(s2, 'XData', dorsum_centroid(:,1), 'YData', dorsum_centroid(:,2))
            end
            
            saveas(fig, fullfile(ExportFolder, [fname, '.png']))
        end
    end

    if ExportPlots
        close(fig)
    end

    function ParseVarargin()
        if ~isempty(varargin)
            nargin = ceil(length(varargin)/2);
            varargin = reshape(varargin, [2, nargin]);
            for n = 1:nargin
                if strcmpi(varargin{1,n},'ExportPlots')
                    ExportPlots = varargin{2,n};
                elseif strcmpi(varargin{1,n},'ExportFolder')
                    ExportFolder = varargin{2,n};
                elseif strcmpi(varargin{1,n},'ShowCentroids')
                    ShowCentroids = varargin{2,n};
                elseif strcmpi(varargin{1,n},'Threshold')
                    Threshold = varargin{2,n};
                else
                    error('%s is an unrecognized input.', varargin{1,n})
                end
            end

            if ExportPlots && isempty(ExportFolder)
                    error('If ExportPlots is true you must give the ExportFolder')
            end
        end
    end
end