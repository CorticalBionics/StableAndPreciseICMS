%%% Fig 5 - Multi-channel
addpath("HelperFunctions")
load(fullfile(DataPath, 'MultiElecData.mat'))
load(fullfile(DataPath, 'MultiChannelData_Processed.mat'))
load(fullfile(DataPath, 'Localization'))
SetFont('Arial', 9)
[palmar_mask, palmar_template, dorsal_mask, dorsal_template] = GetHandMasks();

%% Plot
if exist('fig', 'var') 
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax
end

single_cmap = ColorGradient([.93, .92, 1], [.18, .18, .83]);
multi_cmap = ColorGradient([1, .92, .93], [.83, .18, .18]);

ex_1 = [53,62];
ex_2 = [10,57];

axis_positions1 = [.025 .7 .275 .35;...
                  .35 .7 .275 .35;...
                  .675 .7 .275 .35];
axis_positions2 = [.025 .45 .275 .35;...
                  .35 .45 .275 .35;...
                  .675 .45 .275 .35];

fig = figure('Units', 'inches', 'Position', [27, 1, 6, 5.5]);
ax(1) = axes('Position', axis_positions1(1,:), 'Color', 'none', 'XColor','none','YColor','none');
    PlotPF(palmar_template, axis_positions1(1,:), MultiChannelData, ex_1(1), single_cmap)
    PlotPF(palmar_template, axis_positions1(2,:), MultiChannelData, ex_1(2), single_cmap)
    PlotPF(palmar_template, axis_positions1(3,:), MultiChannelData, ex_1, multi_cmap)

ax(2) = axes('Position', axis_positions2(1,:), 'Color', 'none', 'XColor','none','YColor','none');
    PlotPF(palmar_template, axis_positions2(1,:), MultiChannelData, ex_2(1), single_cmap)
    PlotPF(palmar_template, axis_positions2(2,:), MultiChannelData, ex_2(2), single_cmap)
    PlotPF(palmar_template, axis_positions2(3,:), MultiChannelData, ex_2, multi_cmap)

annotation(fig, "textbox", [.22 .685 .3 .1], 'String', ColorText({'Single Electrode'}, [.18, .18, .83]), ...
    'FontSize', 12, 'EdgeColor', 'none', 'VerticalAlignment','middle');
annotation(fig, "textbox", [.725 .685 .3 .1], 'String', ColorText({'Multi Electrode'}, [.83, .18, .18]), ...
    'FontSize', 12, 'EdgeColor', 'none', 'VerticalAlignment','middle');


ax(3) = axes('Position', [.1 .065 .1 .35]); hold on
    plot([.5 2.5], [0 0], 'Color', [.6 .6 .6], 'LineStyle','--')
    Swarm(1, MultiElecAnalysis.SumCorr, 'Color', SubjectColors('BCI02'), 'DS', 'Violin')
    Swarm(2, MultiElecAnalysis.NullCorr(:), 'Color', [.6 .6 .6], 'DS', 'Violin', 'SPL', 0)

    set(gca, 'XLim', [.5 2.5], 'XTick', [1,2], 'XTickLabel', {'Obs', 'Null'}, 'YLim', [-1 1], ...
        'YTick', [-1:.5:1], 'XTickLabelRotation', 0)
    ylabel('Correlation (r)')


axes('Position', [0.25 0.05 .4 .4]);
    imshow("B:\UserFolders\CharlesGreenspon\BCI_RFPF\Figures\LocalizationReference.png")
    set(gca, 'XTick', [], 'YTick', [])

ds = 'Box';
spl = 0;
sc_color = rgb(158, 158, 158);
mc_color = rgb(244, 67, 54);
axes('Position', [.75 .065 .2 .35]); hold on
    Swarm(1, LocalizationData.SingleChannel, 'GroupName', 'Multi Channel','Color', sc_color,...
            'DS', ds, 'SPL', spl, 'CenterMethod', 'Median', 'ShowStats', true)
    Swarm(2, LocalizationData.MultiChannel, 'GroupName', 'Single Channel','Color', mc_color,...
            'DS', ds, 'SPL', spl, 'CenterMethod', 'Median', 'ShowStats', true)

    text(1, 0.025, {'Single', 'Electrode'}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment','center')
    text(2, 0.025, {'Multi', 'Electrode'}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment','center')

    set(gca, 'XTick', [1, 2], 'XTickLabel', {''}, 'XLim', [.5 2.5], 'YTick', [0:.2:1], 'YLim', [0 1])
    ylabel('Performance')

y1 = .92; y2 = .68; y3 = 0.4;
char_offset = 64;
annotation("textbox", [0.025 y1 .05 .05], 'String', char(char_offset+1), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.025 y2 .05 .05], 'String', char(char_offset+2), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.025 y3 .05 .05], 'String', char(char_offset+3), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.2 y3 .05 .05], 'String', char(char_offset+4), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.675 y3 .05 .05], 'String', char(char_offset+5), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')


%% Helper function
function PlotPF(palmar_template, axis_position, MCData, COI, cmap)
    xl = [80, 800];
    yl = [0, 500];

    % Find channel indices
    idx = 0;
    for i = 1:length(MCData)
        if length(COI) == 2 && length(MCData(i).Channel) == 2 && all(MCData(i).Channel == COI)
            idx = i;
        end
        if length(COI) == 1 && length(MCData(i).Channel) == 1
            if MCData(i).Channel == COI
                idx = i;
            end
        end
    end

    pix_density = zeros(size(palmar_template,1), size(palmar_template,2));

    % Background
    base_axis = axes('Position', axis_position);
    imshow(palmar_template, 'Parent', base_axis)
    
    % single 1
    ax1 = axes('Position', axis_position);
    pix_ids = cat(1, MCData(idx).ProcessedSFs.PalmSF.PixelIds{:});
    pix_vals = cat(1, MCData(idx).ProcessedSFs.PalmSF.NormalizedCount{:});
    pix_density(pix_ids) = pix_vals;
    imagesc(pix_density, 'AlphaData',  (pix_density ~= 0) .* .75);
    colormap(ax1, cmap)
    
    set(base_axis, 'DataAspectRatio', [1 1 1], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', xl, 'YLim', yl)
    set(ax1, 'DataAspectRatio', [1 1 1], 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'XLim', xl, 'YLim', yl)
end