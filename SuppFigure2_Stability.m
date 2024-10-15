%%% Figure 2 - Stability
addpath("HelperFunctions")
load(fullfile(DataPath, 'Stability.mat'))
SetFont('Arial', 9)
SubjectIDs = {'BCI02', 'CRS02', 'CRS07'};
SubjColors = SubjectColors(SubjectIDs);
subj_str = {'C1', 'P2', 'P3'};

%% Supplement
if exist('fig', 'var')
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax processed_image_axis processed_ax
end
fig = figure('Units', 'inches', 'Position', [27.5, 1, 6.48, 6]);

binx = [0:0.05:1];
h = 3; w = 3;
ax(1) = subplot(h,w,1); hold on
    for i = 1:3
        Swarm(i, StabilityData.SubjectPerceptFrequency{i}, "Color", SubjectColors(SubjectIDs{i}), "DistributionStyle", "Box")
    end
    set(gca, 'XLim', [.5 3.5], ...
             'YLim', [0 1], ...
             'XTick', [1:3], ...
             'YTick', [0:.2:1], ...
             'XTickLabel', ColorText(subj_str, SubjColors))
    ylabel('Proportion Observed') 

ax(2) = subplot(h,w,2); hold on
    plot([.33, .33], [0, 1], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for s = 1:3
        AlphaLine(binx, StabilityData.AggregrateThresholdStability{s}, SubjColors(s,:), 'LineWidth', 2)
    end
    set(gca, 'XTick', [0:.5:1])
    xlabel('Aggregate Threshold'); ylabel('p(Electrodes)')

ax(3) = subplot(h,w,3); hold on
    plot([.33, .33], [0, 1], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for s = 1:3
        AlphaLine(binx, StabilityData.AggregrateSize{s}, SubjColors(s,:), 'LineWidth', 2)
    end
    set(gca, 'XTick', [0:.5:1], 'YLim', [0 1])
    xlabel('Aggregate Threshold'); ylabel('Relative Area')

ax(4) = subplot(h,w,4); hold on
    plot([.33, .33], [0, 30], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for s = 1:3
        AlphaLine(binx, StabilityData.AggregratePWD{s} ./ 5, SubjColors(s,:), 'LineWidth', 2)
    end
    set(gca, 'XTick', [0:.5:1])
    xlabel('Aggregate Threshold'); ylabel('Distance (mm)')

ax(5) = subplot(h,w,5); hold on
    [rs,ps] = deal(NaN(3,1));
    for s = 1:3
        x = StabilityData.MatchedSubjectDTs{s};
        y = StabilityData.SubjectPerceptFrequency{s};
        [rs(s), ps(s)] = corr(x,y');
        f = fitlm(x,y,'linear');
        plot([0 100], feval(f, [0 100]), 'Color', SubjColors(s,:), 'LineStyle', '--')
        scatter(x, y, 30, SubjColors(s,:), 'filled')
    end
    xlabel(sprintf('Detection Threshold (%sA)', GetUnicodeChar('mu')))
    ylabel('Proportion Observed')

ax(6) = subplot(h,w,6); hold on
    for s = 1:3
        seq_xy = cat(1, StabilityData.SequentialDistances{s}{:});
        x = seq_xy.Days;
        y = seq_xy.Distance;
        x0 = x == 0;
        x = x(~x0); y = y(~x0);
        bins = linspace(0,max(x),6);
        bin_x = bins(1:end-1) + diff(bins([1,2]))/2; bin_x = bin_x - min(bin_x) + 1;
        [r,p] = corr(x,y, 'Rows', 'complete');
        fprintf('Sequential Correlation, Subject %s: r = %0.2f, %s\n', subj_str{s}, r, pStr(p*3))
        % Discretize
        [~,~,bin_idx] = histcounts(x, bins);
        [y_mean, y_std] = deal(zeros(1,length(bins)-1));
        for i = 1:length(bin_x)
            y_mean(i) = mean(y(bin_idx == i), 'omitnan');
            y_std(i) = std(y(bin_idx == i));
        end
        nan_idx = ~(isnan(y_mean) | isnan(y_std));
        fill([bin_x(nan_idx), fliplr(bin_x(nan_idx))], [y_mean(nan_idx) + y_std(nan_idx), ...
            fliplr(y_mean(nan_idx) - y_std(nan_idx))], SubjColors(s,:), 'EdgeColor', SubjColors(s,:),...
            'FaceAlpha', .1, 'EdgeAlpha', .1)
        plot(bin_x(nan_idx), y_mean(nan_idx), 'Color', SubjColors(s,:), 'LineWidth', 1.5);
    end

    set(gca, 'YLim', [0 60], 'XLim', [1 2000], 'XTick', [1, 10, 100, 1000], 'XScale', 'log')
    xlabel('Days from 1st Survey'); ylabel('Dist from 1st Cent (mm)')

    % Check P2 for first 2 years only
    s = 2;
    seq_xy = cat(1, StabilityData.SequentialDistances{s}{:});
    x = seq_xy.Days;
    y = seq_xy.Distance;
    x0 = x == 0;
    x = x(~x0); y = y(~x0);
    x2year = x < 365 * 2;
    x = x(x2year); y = y(x2year);
    bins = linspace(0,max(x),6);
    bin_x = bins(1:end-1) + diff(bins([1,2]))/2; bin_x = bin_x - min(bin_x) + 1;
    [r,p] = corr(x,y, 'Rows', 'complete');
    fprintf('Sequential Correlation (2-year), Subject %s: r = %0.2f, %s\n', subj_str{s}, r, pStr(p*3))

ax(7) = subplot(h,w,7); hold on
    plot([.5 3.5], [0 0], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for s = 1:3
        Swarm(s, StabilityData.SequentialSlopes{s}, SubjColors(s,:), 'DS', 'Box', 'SPL', 0)
        [~, p, ~, st] = ttest(StabilityData.SequentialSlopes{s});
        fprintf('Sequential slope t-test, Subject %s: t = %0.2f, %s\n', subj_str{s}, st.tstat, pStr(p*3))
    end
    set(gca, 'XLim', [.5 3.5], 'YLim', [-0.025 0.15], 'XTick', [1:3], 'XTickLabel', ColorText(subj_str, SubjColors))
    ylabel('Slope (mm/day)')

ax(8) = subplot(h,w,8); hold on
    plot([0 1], [0 1], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for s = 1:3
        scatter(StabilityData.VectorStrength{s}, StabilityData.BiVectorStrength{s}, 30, SubjColors(s,:), 'filled')
    end
    set(gca, 'XTick', [0 1], 'YTick', [0 1], 'XLim', [0 1], 'YLim', [0 1])
    xlabel('Vector Strength'); ylabel(sprintf('2%s Vector Strength', GetUnicodeChar('pi')))   

AddFigureLabels(gcf, [0.05, 0.0])

export_path = fullfile(DataPath, 'SuppFig2_Stability');
print(gcf, export_path, '-dpng', '-r300')
print(gcf, export_path, '-depsc', '-r300')
