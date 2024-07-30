% Somatotopy figure
addpath("HelperFunctions")
load(fullfile(DataPath, 'RFData.mat'))
load(fullfile(DataPath, 'Discriminability.mat'))
load(fullfile(DataPath, 'SurveyMainData_Processed.mat'))
load(fullfile(DataPath, 'RFMappingData_Processed.mat'))

export_path = fullfile(DataPath, 'Figure3_Somatotopy');
SetFont('Arial', 9)

SubjectIDs = {'BCI02', 'CRS02', 'CRS07'};
subj_str = {'C1', 'P2', 'P3'};
SubjColors = SubjectColors(SubjectIDs);
[palmar_mask, palmar_template, dorsal_mask, dorsal_template] = GetHandMasks();
[pbt_regions, dbt_regions] = GetHandSegments();

%% Load CWRU data
cwru_color = rgb(2, 136, 209);
CWRU_raw = readcell(fullfile(DataPath, 'CWRU_RFPF_Digits.xlsx'));
CWRU_RFPF = cell2mat(CWRU_raw(2:end,[2,3]));
% Filter by both channels having RF/PF
both_idx = all(CWRU_RFPF > 0, 2);
CWRU_RFPF = CWRU_RFPF(both_idx,:);

%% Main
if exist('fig', 'var') 
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax
end

s = 1; % Use BCI02 as example
ChannelMap = LoadSubjectChannelMap(SubjectIDs{s}); 
grid_locations = ChannelMap.ChannelNumbers(ChannelMap.IsSensory);
sidx = strcmp({SurveyMainData.Subject}, SubjectIDs{s});
[grid_tags, palm_seg_present, dors_seg_present] = GetSegmentLabels(SurveyMainData(sidx), ...
    ChannelMap, pbt_regions, dbt_regions, 'DailyCentroid');
grid_cols = [pbt_regions(48).Color; ...
             pbt_regions(41).Color; ...
             pbt_regions(28).Color; ...
             pbt_regions(17).Color; ...
             pbt_regions(1).Color];


fig = figure('Units', 'inches', 'Position', [27, 1, 6.25, 5.5]);

ax(1) = axes('Position', [.0 .675 .25 .3]); hold on
    PlotGridSegments(grid_locations{1}, grid_tags{1,1}, pbt_regions, 0, ax(1), 0)
    % Get grid_loc/grid_tag
    ids = NaN(size(grid_tags{1,1}));
    for i = 1:numel(grid_tags{1,1})
        if ~isempty(grid_tags{1,1}{i})
            ids(i) = str2double(grid_tags{1,1}{i}(2));
        end
    end

    % Plot to the right of the array
    for i = 1:size(ids, 1)
        jj = 4;
        for j = 1:size(ids, 2)
            if ~isnan(ids(i,j))
                scatter(jj, i - 5.5, 30, grid_cols(ids(i,j),:), 'filled')
                jj = jj + 0.7;
            end
        end
    end

    % Plot below the array
    for j = 1:size(ids, 2)
        ii = -6;
        for i = 1:size(ids, 1)
            if ~isnan(ids(i,j))
                scatter(j - 3.5, ii, 30, grid_cols(ids(i,j),:), 'filled')
                ii = ii - 0.7;
            end
        end
    end

    set(ax(1), 'DataAspectRatio', [1,1,1], 'YColor', 'none', 'XColor', 'none', 'YLim', [-10 7], 'XLim', [-5 6.5])


ax(2) = axes('Position', [.375 .725 .25 .225]); hold on
    for s = 1:length(SubjectIDs)
        for a = 1:2
            ads = smoothdata(DiscriminabilityResults.AngDiscrim{s,a});
            [~, max_idx] = max(ads);
            ads = circshift(ads, -round(max_idx - length(DiscriminabilityResults.Angles)/4));
            plot(DiscriminabilityResults.Angles(1:floor(length(DiscriminabilityResults.Angles)/2)), ...
                 ads(1:floor(length(DiscriminabilityResults.Angles)/2)), 'Color', SubjectColors(SubjectIDs{s}), 'LineWidth', 1.5)
        end
    end
    set(gca, 'XLim', [0 180], 'YLim', [.45 1], 'XTick', [0:60:180], 'YTick', [0:.2:1])
    text(75, 0.45, ColorText(subj_str, SubjColors), 'VerticalAlignment', 'bottom', 'FontWeight', 'bold')
    xlabel(sprintf('Rotation (%s)',GetUnicodeChar('Degree')))
    ylabel('Classification Performance')

ax(3) = axes('Position', [.725 .725 .225 .225]); hold on
    plot([0 1], [0 1], 'Color', [.6 .6 .6], 'LineStyle', '--')
    y = zeros(6,1);
    c = zeros(6,3);
    for s = 1:3
        for a = 1:2
            if isempty(DiscriminabilityResults.AngDiscrim{s,a})
                continue
            end
            scatter(DiscriminabilityResults.DigitClassification(s,a), max(DiscriminabilityResults.AngDiscrim{s,a}), 30, SubjColors(s,:),...
                'LineWidth', 2)

        end
    end
    set(gca, 'XLim', [.5 1], 'YLim', [.5 1], 'XTick', [.5 1], 'YTick', [1])
    ylabel('Classification (1D)'); xlabel('Classification (2D)')
    text(.55, 1, ColorText({'C1', 'P2', 'P3'}, SubjColors), 'VerticalAlignment', 'top')

sp_pos = [.0 .325 .25 .275; ...
          .2125 .325 .25 .275; ...
          .4225 .325 .25 .275];
xl = [550, 900; 320, 600; 460, 740];
yl = [700, 1200; 850, 1250; 750, 1150];
idx1 = [2,4];
idx2 = [4,2];
for s = 1:3
    [a,b,c] = ShowSubjectImplant(SubjectIDs{s}, 'SubjectString', '', 'ShowScale', false, 'FontSize', 0, 'SplineStyle', 'none');
    copyobj(b, fig)
    close(a)
    set(fig.Children(1), 'Position', sp_pos(s,:))
    % Add alignment angle thing
    for a = 1:2
        if s == 1
            total_angle = DiscriminabilityResults.OpimAngle(s,a) + c.ArrayRotations(idx1(a));
            x_mean = mean(c.ArrayCoordinates{idx1(a),1});
            y_mean = mean(c.ArrayCoordinates{idx1(a),2});
        else
            total_angle = DiscriminabilityResults.OpimAngle(s,a) + c.ArrayRotations(idx2(a));
            x_mean = mean(c.ArrayCoordinates{idx2(a),1});
            y_mean = mean(c.ArrayCoordinates{idx2(a),2});
        end
        x3 = [0 0];
        y = [-75 75];
        [xr, yr] = rotate_coords(x3,y,total_angle);
        xo = xr + x_mean;
        yo = yr + y_mean;
        
        plot(xo, yo, 'Color', 'k', 'LineStyle', ':', 'Parent', fig.Children(1), 'LineWidth', 2)
        set(fig.Children(1), 'XLim', xl(s,:), 'YLim', yl(s,:), 'DataAspectRatio', [1 1 1])
    end
end

ax(4) = axes('Position', [.725 .375 .225 .2]); hold on
    bins = linspace(0,4,8);
    bin_x = bins(1:end-1) + range(bins([1,2]))/2;
    [x2, y2] = deal(cell(3,2));
    for s = 1:3
        for a = 1:2
            if s == 3 && a == 2
                continue
            end
            x2{s,a} = DiscriminabilityResults.ArrayDists{s,a}  .* 0.4;
            y2{s,a} = mean(cat(3, DiscriminabilityResults.PerceptDistsRaw{s,a,:}), 3, 'omitnan') ./ 5;
        end
    end
    x_cat = cellfun(@(c) c(:), x2, 'UniformOutput', false); x_cat = cat(1, x_cat{:});
    y_cat = cellfun(@(c) c(:), y2, 'UniformOutput', false); y_cat = cat(1, y_cat{:});
    [ps,rs] = deal(length(SubjectIDs),2);
    for s = 1:length(SubjectIDs)
        for a = 1:2
            if s == 3 && a == 2
                continue
            end
            [~,~,bin_idx] = histcounts(DiscriminabilityResults.ArrayDists{s,a}(:) .* 0.4, bins);
            % Discretize
            [b_mean, b_std] = deal(zeros(1,length(bin_x)));
            for b = 1:length(bin_x)
                b_idx = bin_idx == b;
                b_all = y2{s,a}(b_idx);
                b_mean(b) = mean(b_all, 'omitnan');
                b_std(b) = std(b_all, 'omitnan') / sqrt(sum(~isnan(b_all)));
            end
            nan_idx = isnan(b_mean);
            plot(bin_x(~nan_idx), b_mean(~nan_idx), 'Color', SubjectColors(SubjectIDs{s}), 'LineWidth', 1, 'LineStyle', '--')
            [rs(s,a),ps(s,a)] = corr(DiscriminabilityResults.ArrayDists{s,a}(:),...
                         DiscriminabilityResults.PerceptDistsRaw{s,a}(:), 'rows','complete');
        end
    end
    % Discretize
    [~,~,bin_idx] = histcounts(x_cat, bins);
    [b_mean, b_std] = deal(zeros(1,length(bin_x)));
    for b = 1:length(bin_x)
        b_idx = bin_idx == b;
        b_all = y_cat(b_idx);
        b_mean(b) = mean(b_all, 'omitnan');
        b_std(b) = std(b_all, 'omitnan');
    end
    nan_idx = isnan(b_mean);
    fill([bin_x(~nan_idx), fliplr(bin_x(~nan_idx))], [b_mean(~nan_idx) + b_std(~nan_idx),...
          fliplr(b_mean(~nan_idx) - b_std(~nan_idx))], [.6 .6 .6], 'EdgeColor', [.6 .6 .6],...
                 'FaceAlpha', .1, 'EdgeAlpha', .1)
    plot(bin_x(~nan_idx), b_mean(~nan_idx), 'Color', [.6 .6 .6], 'LineWidth', 3)    
    set(gca, 'XLim', [0 4.1], 'YLim', [0 40], 'XTick', [0:1:4], 'YTick', [0:10:50])
    xlabel('Cortical Distance (mm)'); ylabel('Sensation Distance (mm)')
    ps = HolmBonferroni(ps);

subj = 'CRS07';
ch = 18;
axis_positions = [.01 .025 .275 .275];
ax(5) = axes('Position', axis_positions(1,:), 'Color', 'none', 'XColor','none','YColor','none');
        pf_idx = find([SurveyMainData.Channel] == ch & strcmp({SurveyMainData.Subject}, subj));
        PFTable = SurveyMainData(pf_idx).ProcessedSFs.PalmSF;
        rf_idx = find([RFMappingData.Channel] == ch & strcmp({RFMappingData.Subject}, subj));
        RFTable = RFMappingData(rf_idx).ProcessedSFs.PalmSF;
        PlotRFPF(palmar_template, axis_positions, PFTable, RFTable);

annotation(fig, "textbox", [.225 .125 .075 .1], 'String', ColorText({'PF', 'RF'}, [.83, .18, .18; .18, .18, .83]), 'FontSize', 9,...
    'EdgeColor', 'none')

subj = {'BCI02', 'CRS07'};
ax(6) = axes('Position', [.375 .085 .15 .175]); hold on
    plot([0.1 200], [.1 200], 'Color', [.6 .6 .6], 'LineStyle','--')
    [area_r, area_rp, area_t, area_tp] = deal(zeros(2,1));
    for i = 1:2
        subj_idx = strcmp({RFMappingData.Subject}, subj{i});
        scatter(RFAnalysis.PFArea(subj_idx,1) / 100, RFAnalysis.RFArea(subj_idx,1) / 100, 30, SubjectColors(subj{i}), 'filled')
        [area_r(i), area_rp(i)] = corr(RFAnalysis.PFArea(subj_idx,1), RFAnalysis.RFArea(subj_idx,1), 'rows','complete');
        [~ ,area_tp(i), tstats] = ttest(RFAnalysis.PFArea(subj_idx,1), RFAnalysis.RFArea(subj_idx,1));

    end
    set(gca, 'XLim', [0.1 200], 'YLim', [0.1 200], 'XTick', [.1, 1, 10, 100], 'YTick', [.1, 1, 10, 100],...
        'XScale', 'log', 'YScale', 'log', 'XTickLabelRotation', 0)
    xlabel('PF Area (cm^2)'); ylabel('RF Area (cm^2)')
 
ax(7) = axes('Position', [.585 .085 .15 .175]); hold on
    p25 = zeros(2,1);
    for i = 1:2
        subj_idx = strcmp({RFMappingData.Subject}, subj{i});
        Swarm(i, RFAnalysis.PFInclusion(subj_idx,1), 'Color', SubjectColors(subj{i}), 'DS', 'Box')
        p25(i) = prctile(RFAnalysis.PFInclusion(subj_idx,1), 25);
    end
    set(gca, 'XLim', [0.5 2.5], 'XTick', [1,2], 'XTickLabel', ColorText(subj_str([1,3]), SubjColors([1,3],:)), ...
        'YTick', [0 1])
    ylabel('PF Inclusion')

ax(7) = axes('Position', [.8 .085 .15 .175]); hold on
    plot([0 4], [0 0], 'Color', [.6 .6 .6], 'LineStyle', '--')
    ii = 1;
    for i = 1:2
        subj_idx = strcmp({RFMappingData.Subject}, subj{i});
        Swarm(ii, RFAnalysis.RFDigits(subj_idx,1) - RFAnalysis.PFDigits(subj_idx,1), 'Color', SubjectColors(subj{i}), 'DS', 'Box')
        [p,h] = ranksum(RFAnalysis.RFDigits(subj_idx,1), RFAnalysis.PFDigits(subj_idx,1));
        if h
            text(ii, 4, '*', 'Color', SubjectColors(subj{i}), ...
                'FontSize', 20, 'VerticalAlignment', 'middle', 'HorizontalAlignment','center')
        end
        ii = ii + 1;
    end

    % Add CWRU
    Swarm(ii, CWRU_RFPF(:,1) - CWRU_RFPF(:,2), 'Color', cwru_color, 'DS', 'Box')
    % Ranksum test
    [p,h] = ranksum(CWRU_RFPF(:,2), CWRU_RFPF(:,1));
    if h
        text(ii, 4, '*', 'Color', cwru_color, ...
            'FontSize', 20, 'VerticalAlignment', 'middle', 'HorizontalAlignment','center')
    end
    set(gca, 'XLim', [0.5, 3.5], 'XTick', [1:3], 'XTickLabel', ColorText({'C1', 'P3', 'R1'}, ...
        [SubjectColors('BCI02'); SubjectColors('CRS07'); cwru_color]), 'YLim', [-4 4])
    ylabel('#RF - #PF')


y1 = .94; y2 = .575; y3 = 0.25;
char_offset = 64;
annotation("textbox", [0 y1 .05 .05], 'String', char(char_offset+1), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.325 y1 .05 .05], 'String', char(char_offset+2), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.65 y1 .05 .05], 'String', char(char_offset+4), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0 y2 .05 .05], 'String', char(char_offset+3), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.68 y2 .05 .05], 'String', char(char_offset+5), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0 y3 .05 .05], 'String', char(char_offset+6), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.3 y3 .05 .05], 'String', char(char_offset+7), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.525 y3 .05 .05], 'String', char(char_offset+8), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.75 y3 .05 .05], 'String', char(char_offset+9), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')

print(fig, export_path, '-dsvg', '-r300')

%% Helper functions
function [xr,yr] = rotate_coords(x,y,r)
    r = deg2rad(r);
    xr = x*cos(r) - y*sin(r);
    yr = y*cos(r) + x*sin(r);
end

function PlotRFPF(palmar_template, axis_position, PFTable, RFTable)
    pf_cmap = [linspace(1,.83,255); linspace(.92,.18,255); linspace(.93,.18,255)]';
    rf_cmap = [linspace(.93,.18,255); linspace(.92,.18,255); linspace(1,.83,255)]';
    xl = [80, 1000];
    yl = [0, 800];
    % Background
    base_axis = axes('Position', axis_position);
    imshow(palmar_template, 'Parent', base_axis)
    
    % RF first (because it's bigger)
    rf_axis = axes('Position', axis_position);
    rf_density = zeros(size(palmar_template,1), size(palmar_template,2));
    pix_ids = cat(1, RFTable.PixelIds{:});
    pix_vals = cat(1, RFTable.NormalizedCount{:});
    rf_density(pix_ids) = pix_vals;
    imagesc(rf_density, 'AlphaData',  (rf_density ~= 0) .* .75);
    colormap(rf_axis, rf_cmap)

    % PF second
    pf_axis = axes('Position', axis_position);
    pf_density = zeros(size(palmar_template,1), size(palmar_template,2));
    pix_ids = cat(1, PFTable.PixelIds{:});
    pix_vals = cat(1, PFTable.NormalizedCount{:});
    pf_density(pix_ids) = pix_vals;
    imagesc(pf_density, 'AlphaData',  (pf_density ~= 0) .* .75);
    colormap(pf_axis, pf_cmap)

    set(base_axis, 'DataAspectRatio', [1 1 1], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', xl, 'YLim', yl)
    set(rf_axis, 'DataAspectRatio', [1 1 1], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', xl, 'YLim', yl)
    set(pf_axis, 'DataAspectRatio', [1 1 1], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', xl, 'YLim', yl)
end