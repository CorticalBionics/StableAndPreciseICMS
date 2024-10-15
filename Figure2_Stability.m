%%% Figure 2 - Stability
addpath("HelperFunctions")
load(fullfile(DataPath, 'SurveyMainData_Processed.mat'))
load(fullfile(DataPath, 'Stability.mat'))
[palmar_mask, palmar_template, dorsal_mask, dorsal_template] = GetHandMasks();
[pbt_regions, dbt_regions] = GetHandSegments();
orig_size = [1200, 1050];
SetFont('Arial', 9)
SubjectIDs = {'BCI02', 'CRS02', 'CRS07'};
SubjColors = SubjectColors(SubjectIDs);
subj_str = {'C1', 'P2', 'P3'};

xl = [80, 1000];
yl = [0, 800];

[palmar_mask, palmar_template, dorsal_mask, dorsal_template] = GetHandMasks();
palm_thick = mean(palmar_template,3);
palm_thick = bwmorph(~palm_thick, 'thicken', 3);
palm_thick = uint8(repmat(~palm_thick,[1,1,3])) .* 255;
SegmentationStrings = {'D1', 'D2', 'D3', 'D4', 'D5', 'W', 'P'};
LT = ProcessSensoryFields(SurveyMainData, SegmentationStrings, 'Threshold', 0);
HT = ProcessSensoryFields(SurveyMainData, SegmentationStrings, 'Threshold', .33);
si = [1,3,2];

%% Total surveys + survey per-electrode
num_surveys = zeros(size(SurveyMainData));
for i = 1:length(SurveyMainData)
    num_surveys(i) = length(SurveyMainData(i).SensoryField);
end
surveys_per_subject = zeros(3,1);
for s = 1:3
    s_idx = find(strcmp({SurveyMainData.Subject}, SubjectIDs{s}));
    surveys_per_subject(s) = sum(num_surveys(s_idx));
end

%% Get all areas
areas_all = zeros(size(SurveyMainData));
for i = 1:length(SurveyMainData)
    areas_all(i) = sum(SurveyMainData(i).ProcessedSFs.PalmSF.Area) + sum(SurveyMainData(i).ProcessedSFs.DorsSF.Area);
end
areas_all(areas_all < 1) = NaN; 
areas_all = areas_all ./ 100; % Convert to cm

[rf_density_low, rf_density_high] = deal(cell(3,1));
for s = 1:3
    s_idx = find(strcmp({SurveyMainData.Subject}, SubjectIDs{s}));
    % Low
    pix_cell = cell(length(s_idx),1);
    for i = 1:length(s_idx)
        pix_cell{i} = cat(1,LT(s_idx(i)).ProcessedSFs.PalmSF.PixelIds{:});
    end
    pix_all = cat(1,pix_cell{:});
    [GC,GR,GP] = groupcounts(pix_all);
    rf_mat = zeros(1200, 1050);
    rf_mat(GR) = GP;
    rf_density_low{s}= rf_mat;
    % High
    pix_cell = cell(length(s_idx),1);
    for i = 1:length(s_idx)
        pix_cell{i} = cat(1,SurveyMainData(s_idx(i)).ProcessedSFs.PalmSF.PixelIds{:});
    end
    pix_all = cat(1,pix_cell{:});
    [GC,GR,GP] = groupcounts(pix_all);
    rf_mat = zeros(1200, 1050);
    rf_mat(GR) = GP;
    rf_density_high{s}= rf_mat;
end

%% Main 1
sub_ch = [42; 34; 39];
if exist('fig', 'var')
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax processed_image_axis processed_ax
end
fig = figure('Units', 'inches', 'Position', [27.5, 1, 6.48, 5.5]);

ax(1) = axes('Position', [.0 .5 .325 .6]);
    imshow(palmar_template); hold on % Background
    for c = 1:size(SubjectIDs,2)
        sub_ch_idx = find(strcmp({SurveyMainData.Subject}, SubjectIDs{c}) & [SurveyMainData.Channel] == sub_ch(c));
        last_obs = find(cellfun(@(s) size(s,1), {SurveyMainData(sub_ch_idx).SensoryField.PalmSF}),1,'last');
        xy = SurveyMainData(sub_ch_idx).SensoryField(last_obs).PalmSF(1).Boundary;
        patch(xy(:,1), xy(:,2), SubjColors(c,:), 'FaceAlpha', .5, 'EdgeColor', SubjColors(c,:))
        scatter(SurveyMainData(sub_ch_idx).SensoryField(last_obs).PalmSF(1).Centroid(1),...
                SurveyMainData(sub_ch_idx).SensoryField(last_obs).PalmSF(1).Centroid(2),...
                50, SubjColors(c,:), 'x', 'LineWidth', 2)
    end

    set(ax(1), 'DataAspectRatio', [1 1 1], 'XLim', xl, 'YLim', yl, 'XColor', 'none', 'YColor', 'none')
    % Legend
    annotation("textbox", [.25 .75 .1 .1], 'String', ColorText(subj_str, SubjColors), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')


ax(2) = axes('Position', [.335 .5 .325 .6]);
    imshow(palmar_template); hold on % Background
    for c = 1:size(SubjectIDs,2)
        sub_ch_idx = find(strcmp({SurveyMainData.Subject}, SubjectIDs{c}) & [SurveyMainData.Channel] == sub_ch(c));
        obs = find(cellfun(@(s) size(s,1), {SurveyMainData(sub_ch_idx).SensoryField.PalmSF}) > 0);
        for o = 1:length(obs)
            xy = SurveyMainData(sub_ch_idx).SensoryField(obs(o)).PalmSF(1).Boundary;
            patch(xy(:,1), xy(:,2), SubjColors(c,:), 'FaceAlpha', .1, 'EdgeColor', SubjColors(c,:))
        end
    end

    set(ax(2), 'DataAspectRatio', [1 1 1], 'XLim', xl, 'YLim', yl, 'XColor', 'none', 'YColor', 'none')


ax(3) = axes('Position', [.675 .5 .325 .6], 'XColor', 'none', 'YColor', 'none');
    imshow(palmar_template,'Parent', ax(3)); hold on
    thresh_values = [0:.1:1];
    for c = 1:3
        temp_cmap = ColorGradient(SubjColors(c,:), SubjColors(c,:));
        sub_ch_idx = find(strcmp({SurveyMainData.Subject}, SubjectIDs{c}) & [SurveyMainData.Channel] == sub_ch(c));
        processed_ax(c) = axes('Position', ax(3).Position);
        processed_image = zeros(size(palmar_mask));
        for t = thresh_values
            temp = ProcessSensoryFields(SurveyMainData(sub_ch_idx), SegmentationStrings, 'Threshold', t);
            pix_ids = cat(1, temp.ProcessedSFs.PalmSF.PixelIds{:});
            processed_image(pix_ids) = processed_image(pix_ids) + 1;
        end
        processed_image = processed_image / (length(thresh_values)-3); processed_image(processed_image > 0 & processed_image < 0.1) = 0.1;
        processed_image_axis(c) = imagesc(ones(size(palmar_mask)), 'AlphaData',  processed_image); hold on %#ok<*SAGROW> 
        set(processed_ax(c), 'DataAspectRatio', [1 1 1], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', xl, 'YLim', yl)
        colormap(processed_ax(c), temp_cmap)
    end
    set(ax(3), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', xl, 'YLim', yl)

ax(4) = axes('Position', [0 .3 .225 .325], 'XColor', 'none', 'YColor', 'none');
    imshow(palm_thick,'Parent', ax(4)); hold on
    processed_ax = axes('Position', ax(4).Position);
    processed_image = zeros(size(palmar_mask));
    % Low threshold
    for i = 1:length(SurveyMainData)
        if ~strcmp(SurveyMainData(i).Subject, SubjectIDs{1})
            continue
        end
        temp_px_idx = cat(1, LT(i).ProcessedSFs.PalmSF.PixelIds{:});
        processed_image(temp_px_idx) = 0.25;
    end

    % High threshold
    for i = 1:length(SurveyMainData)
        if ~strcmp(SurveyMainData(i).Subject, SubjectIDs{1})
            continue
        end
        temp_px_idx = cat(1, SurveyMainData(i).ProcessedSFs.PalmSF.PixelIds{:});
        processed_image(temp_px_idx) = .75;
    end
    % Plot
    imagesc(processed_image, 'AlphaData',  processed_image)
    colormap(processed_ax, ColorGradient(SubjColors(1,:), SubjColors(1,:)))
    set(processed_ax, 'DataAspectRatio', [1 1 1], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', xl, 'YLim', yl)
    set(ax(4), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', xl, 'YLim', yl)
        

ax(5) = axes('Position', [.22 .3 .225 .325], 'XColor', 'none', 'YColor', 'none');
    imshow(palm_thick,'Parent', ax(5)); hold on
    processed_ax = axes('Position', ax(5).Position);
    processed_image = zeros(size(palmar_mask));
    % Low threshold
    for i = 1:length(SurveyMainData)
        if ~strcmp(SurveyMainData(i).Subject, SubjectIDs{3})
            continue
        end
        temp_px_idx = cat(1, LT(i).ProcessedSFs.PalmSF.PixelIds{:});
        processed_image(temp_px_idx) = 0.25;
    end

    % High threshold
    for i = 1:length(SurveyMainData)
        if ~strcmp(SurveyMainData(i).Subject, SubjectIDs{3})
            continue
        end
        temp_px_idx = cat(1, SurveyMainData(i).ProcessedSFs.PalmSF.PixelIds{:});
        processed_image(temp_px_idx) = .75;
    end
    % Plot
    imagesc(processed_image, 'AlphaData',  processed_image)
    colormap(processed_ax, ColorGradient(SubjColors(3,:), SubjColors(3,:)))
    set(processed_ax, 'DataAspectRatio', [1 1 1], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', xl, 'YLim', yl)
    set(ax(5), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', xl, 'YLim', yl)

ax(6) = axes('Position', [.44 .3 .225 .325], 'XColor', 'none', 'YColor', 'none');
    imshow(palm_thick,'Parent', ax(6)); hold on
    processed_ax = axes('Position', ax(6).Position);
    processed_image = zeros(size(palmar_mask));
    % Low threshold
    for i = 1:length(SurveyMainData)
        if ~strcmp(SurveyMainData(i).Subject, SubjectIDs{2})
            continue
        end
        temp_px_idx = cat(1, LT(i).ProcessedSFs.PalmSF.PixelIds{:});
        processed_image(temp_px_idx) = 0.25;
    end

    % High threshold
    for i = 1:length(SurveyMainData)
        if ~strcmp(SurveyMainData(i).Subject, SubjectIDs{2})
            continue
        end
        temp_px_idx = cat(1, SurveyMainData(i).ProcessedSFs.PalmSF.PixelIds{:});
        processed_image(temp_px_idx) = .75;
    end
    % Plot
    imagesc(processed_image, 'AlphaData',  processed_image)
    colormap(processed_ax, ColorGradient(SubjColors(2,:), SubjColors(2,:)))
    set(processed_ax, 'DataAspectRatio', [1 1 1], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', xl, 'YLim', yl)
    set(ax(6), 'DataAspectRatio', [1 1 1], 'XColor', 'w', 'YColor', 'w', 'YDir', 'reverse', 'XLim', xl, 'YLim', yl)

ax(7) = axes('Position', [.735 .34 .075 .225]); hold on
    for s = 1:3
        Swarm(s, sum(rf_density_high{si(s)} > 0, 'all') / 25 / 100, 'Color', SubjColors(si(s),:), 'DS', 'Bar', 'SPL', 0, 'EW', 0)
    end
    ylabel('Area (cm^2)')
    set(gca, 'XLim', [.25 3.75], 'YLim', [0 40], 'YTick', [0:10:800], 'XTick', [], 'XTickLabels', {})

ax(8) = axes('Position', [.85 .34 .1 .225]);
    for s = 1:3
        s_idx = strcmp({SurveyMainData.Subject}, SubjectIDs{si(s)});
        Swarm(s, areas_all(s_idx), SubjColors(si(s),:), 'DS', 'Box', 'SPL', 0, 'ShowStats', true)
    end
    set(gca, 'XTick', [], 'XLim', [.25 3.75], 'YLim', [0 15], 'YTick', [0:5:15])

ax(9) = axes('Position', [.075 .075 .2 .175]); hold on
    for s = 1:3
        seq_xy = cat(1, StabilityData.WeightedDistances{s}{:});
        x = seq_xy.Datenum - min(seq_xy.Datenum);
        y = seq_xy.Distance;
        bins = linspace(0,max(x),6);
        bin_x = bins(1:end-1) + diff(bins([1,2]))/2; bin_x = bin_x - min(bin_x) + 1;
        [r,p] = corr(x,y, 'Rows', 'complete');
        fprintf('Weighted Correlation, Subject %s: r = %0.2f, %s\n', subj_str{s}, r, pStr(p*3))
        % Discretize
        [~,~,bin_idx] = histcounts(x, bins);
        [y_mean, y_std] = deal(zeros(1,length(bins)-1));
        for i = 1:length(bin_x)
            y_mean(i) = mean(y(bin_idx == i), 'omitnan');
            y_std(i) = std(y(bin_idx == i));
        end
        nan_idx = ~isnan(y_mean);
        fill([bin_x(nan_idx), fliplr(bin_x(nan_idx))], [y_mean(nan_idx) + y_std(nan_idx), ...
            fliplr(y_mean(nan_idx) - y_std(nan_idx))], SubjColors(s,:), 'EdgeColor', SubjColors(s,:),...
            'FaceAlpha', .1, 'EdgeAlpha', .1)
        plot(bin_x(nan_idx), y_mean(nan_idx), 'Color', SubjColors(s,:), 'LineWidth', 1.5);
    end

    set(gca, 'YLim', [0 20], 'XLim', [1 2000], 'XTick', [1, 10, 100, 1000], 'YTick', [0:10:20], 'XScale', 'log', 'XTickLabelRotation', 0)
    xlabel('Time from 1st Survey (Days)'); ylabel('Distance (mm)')

ax(10) = axes('Position', [.4 .075 .2 .2]); hold on
    y = cat(1, StabilityData.WeightedDistanceMeans{:});
    f = fitlm(log(areas_all), log(y));
    [r,p] = corr(log(areas_all)', log(y), 'rows', 'complete');
    fprintf('Size x variablility correlation: r = %0.2f, %s\n', r, pStr(p))
    plot([.1 75], feval(f, [.1 75]), 'Color', [.6 .6 .6], 'LineStyle', '--')
    for s = 1:3
        s_idx = strcmp({SurveyMainData.Subject}, SubjectIDs{si(s)});
        scatter(areas_all(s_idx), StabilityData.WeightedDistanceMeans{s}, 30, SubjColors(s,:), 'filled')
    end
    set(gca, 'YLim', [1 20], 'YTick', [2, 20], 'XLim', [.1 25], 'XTick', [.1, 1, 10], 'XScale', 'log', 'YScale', 'log')
    x = xlabel('Area (cm^2)'); ylabel('Distance (mm)')
    x.Position(2) = x.Position(2) + 0.1;
    x = cat(1, StabilityData.PerceptArea{:});
    [r,p] = corr(x, y, 'rows', 'complete');

ax(11) = axes('Position', [.7 .075 .2 .2]); hold on
    plot([0 20], [0 20], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for s = 1:3
        s_idx = strcmp(StabilityData.IntraSubj, SubjectIDs{s});
        scatter(StabilityData.IntraDayDist(s_idx), StabilityData.IntraDayDistComp(s_idx), 50, SubjColors(s,:), 'filled')
    end
    set(gca, 'XLim', [0 20], 'YLim', [0 20])
    xlabel('Within (mm)'); ylabel('Across (mm)')
    [r,p] = corr(StabilityData.IntraDayDist, StabilityData.IntraDayDistComp, 'Rows', 'complete');
    fprintf('Intra vs Inter day correlation: r = %0.2f, %s\n', r, pStr(p))
    [~, p, ~, st] = ttest(StabilityData.IntraDayDist, StabilityData.IntraDayDistComp);
    fprintf('Intra vs Inter variability t-test: t = %0.2f, %s\n',  st.tstat, pStr(p))

% Labels
y1 = .93; y2 = .55; y3 = .265;
char_offset = 64;
annotation("textbox", [0 y1 .05 .05], 'String', char(char_offset+1), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [.33 y1 .05 .05], 'String', char(char_offset+2), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [.65 y1 .05 .05], 'String', char(char_offset+3), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0 y2 .05 .05], 'String', char(char_offset+4), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [.665 y2 .05 .05], 'String', char(char_offset+5), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [.79 y2 .05 .05], 'String', char(char_offset+6), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0 y3 .05 .05], 'String', char(char_offset+7), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [.325 y3 .05 .05], 'String', char(char_offset+8), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [.625 y3 .05 .05], 'String', char(char_offset+9), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
shg

% Export
export_path = fullfile(DataPath, 'Figure2_Stability');
% Export
print(gcf, export_path, '-dpng', '-r300')
print(gcf, export_path, '-depsc', '-r300')