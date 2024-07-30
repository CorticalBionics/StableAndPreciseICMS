%%% Figure 1 - Percept maps
addpath("HelperFunctions")
load(fullfile(DataPath, 'SurveyMainData_Processed.mat'))
export_path = fullfile(DataPath, 'Figure1_Intro');
SetFont('Arial', 9)
SubjectIDs = {'BCI02', 'CRS02', 'CRS07'};
subj_str = {'C1', 'P2', 'P3'};
[pbt_regions, dbt_regions] = GetHandSegments();

%% Get implant axes
if exist('fig', 'var')
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax
end
fig = figure('Units', 'inches', 'Position', [27.5, 1, 6.48, 7.5]);

implant_axes = [];
for s = 1:length(subj_str)
    if s == 1
        ss = true;
    else
        ss = false;
    end
    [f, a] = ShowSubjectImplant(SubjectIDs{s}, 'FontSize', 9, 'SubjectString', '',...
        'ShowScale', ss, 'ArrayColors', rgb(96, 125, 139), 'ShowCorner', false);
    copyobj(a, fig) 
    close(f)
end
set(fig.Children(1), 'Position', [.05 .05 .3 .275])
set(fig.Children(2), 'Position', [.05 .3625 .3 .275])
set(fig.Children(3), 'Position', [.05 .675 .3 .275])
% Labels
annotation('textbox',[.35, .85, .1, .1] , 'String', 'C1', 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none')
annotation('textbox',[.35, .5375, .1, .1] , 'String', 'P2', 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none')
annotation('textbox',[.35, .2250, .1, .1] , 'String', 'P3', 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none')

y = [.6 .725 .675];
for s = 1:length(subj_str)
    ChannelMap = LoadSubjectChannelMap(SubjectIDs{s}); 
    [grid_tags, palm_seg_present, dors_seg_present] = GetSegmentLabels(SurveyMainData(strcmp({SurveyMainData.Subject}, SubjectIDs{s})), ...
        ChannelMap, pbt_regions, dbt_regions, 'DailyCentroid');
    temp_palm_seg_present = palm_seg_present;
    f = PlotSegmentedHand(ChannelMap, grid_tags, pbt_regions, temp_palm_seg_present, dbt_regions, dors_seg_present);
    temp_children = get(f, 'Children');
    temp_children = temp_children(4:6); % Palmar only
    copyobj(temp_children, fig)
    close(f)
    if s == 1
        set(fig.Children(1), 'Position', [0.8    y(1)    0.150    0.3], 'Clipping', 'off')
        set(fig.Children(2), 'Position', [0.8    y(2)    0.150    0.3000], 'Clipping', 'off')
    else
        set(fig.Children(1), 'Position', [0.8    y(2)    0.150    0.3], 'Clipping', 'off')
        set(fig.Children(2), 'Position', [0.8    y(1)    0.150    0.3000], 'Clipping', 'off')
    end
    set(fig.Children(3), 'Position', [0.375 y(3) 0.45 .275], 'YLim', [5, 750], 'XLim', [50, 1000])
    y = y - 0.3125;
end

shg

% Export
print(fig, export_path, '-dsvg', '-r300')