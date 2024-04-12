%%% All subject PFs
load(fullfile(DataPath, 'SurveyMainData_Processed.mat'))
[palmar_mask, palmar_template, dorsal_mask, dorsal_template] = GetHandMasks();
[pbt_regions, dbt_regions] = GetHandSegments();

%% Generate temporary high-resolution frames
u_subjects = unique({SurveyMainData.Subject});
sub_labels = {'C1', 'P2', 'P3'};
[frames, ffnames] = deal(cell(size(u_subjects)));
for s = 1:length(u_subjects)
    ChannelMap = LoadSubjectChannelMap(u_subjects{s}); 
    [grid_tags, palm_seg_present, dors_seg_present] = GetSegmentLabels(SurveyMainData(strcmp({SurveyMainData.Subject}, u_subjects{s})), ...
    ChannelMap, pbt_regions, dbt_regions, 'DailyCentroid');
    PlotSegmentedHand(ChannelMap, grid_tags, pbt_regions, palm_seg_present, dbt_regions, dors_seg_present);
    annotation("textbox", [.01 .9 .2 .1], "String", sub_labels{s}, 'EdgeColor', 'None', 'FontSize', 28, 'FontWeight', 'bold');
    % Make a temporary file to get around screen size limits
    ffnames{s} = fullfile(DataPath, sprintf('Supplement1_%s', sub_labels{s}));
    print(gcf, ffnames{s}, '-dpng', '-r300')
    % Import the file 
    frames{s} = imread([ffnames{s}, '.png']);
    delete(fullfile([ffnames{s}, '.png']))
end
close all

full_frame = cat(1, frames{:}); % Concatenate images vertically
clf; 
a1 = axes('Position', [0 0 1 1]);
image(full_frame, 'Parent', a1)
set(gca, 'XColor', 'none', 'YColor', 'none', 'DataAspectRatio', [1 1 1])
set(gcf, 'Position', [50, 50, 1000, 1300])
