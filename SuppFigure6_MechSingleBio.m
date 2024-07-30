% Imports
addpath("HelperFunctions")
load(fullfile(DataPath, "WeberData.mat")) % Weber fraction
load(fullfile(DataPath, "AmpDiscrim.mat")) % ICMS amp discrim
load(fullfile(DataPath, "AmpDiscrimBioLin")) % Amp discrim bio vs lin
load(fullfile(DataPath, "MechAmpDiscrim.mat"))
load(fullfile(DataPath, "MagEstBiovsLin.mat"))
load(fullfile(DataPath, "ForceEquivMulti.mat")) % Mech amp est
load(fullfile(DataPath, "DetectionDataProcessed_BCI02.mat")); % Detection data for BCI02

sigfun = @(c,x) 1./(1 + exp(-c(1).*(x-c(2)))); % c(1) = rate of change, c(2) = x-
invsig = @(c,y) (log(1/y - 1)/-c(1)) + c(2);
force_ratio = 0.37;

SetFont('Arial', 9)
scatter_size = 40;
line_width = 1.5;

%% Data extraction
crs02b_dt = [8, 7.6;...
             14, 6;...
             19, 12.8;...
             28, 48.2;...
             40, 12.4]; % electrode, dt
crs07_dt = [2, 10.2;...
            12, 14.4;...
            16, 36.8;...
            46, 29;...
            62, 22.8]; % electrode, dt

% Detection and JND
[icms_constant_jnd, dt_values1] = deal(NaN([length(AS4EData), 1]));
jnd_scatter_colors1 = zeros([length(AS4EData), 3]);

for e = 1:length(AS4EData)
    j1_idx = find(AS4EData(e).ISISigmoids.ISI == 1); % 1 second ISI
    icms_constant_jnd(e) = AS4EData(e).ISISigmoids.JND(j1_idx,1);
    jnd_scatter_colors1(e,:) = SubjectColors(AS4EData(e).ParticipantID);
    
    % Detection threshold
    if strcmp(AS4EData(e).ParticipantID, 'BCI02')
        tidx = find([ProcessedDetectionData.Channel] == AS4EData(e).Channel);
        dt_values1(e) = ProcessedDetectionData(tidx).MeanThreshold;
    elseif strcmp(AS4EData(e).ParticipantID, 'CRS02b')
        tidx = find(crs02b_dt(:,1) == AS4EData(e).Channel);
        dt_values1(e) = crs02b_dt(tidx,2);
    elseif strcmp(AS4EData(e).ParticipantID, 'CRS07')
        tidx = find(crs07_dt(:,1) == AS4EData(e).Channel);
        dt_values1(e) = crs07_dt(tidx,2);
    end
end

mech_jnd = NaN(length(AS4mData),1);
for i = 1:length(AS4mData)
    mech_jnd(i) = abs((invsig(AS4mData(i).SigSummary.Mean,.75) - invsig(AS4mData(i).SigSummary.Mean,.25))/2);
end
as4m_test_idx = strcmp({AS4mData.Subject}, 'BCI02') & strcmp({AS4mData.Hand}, 'Contra');
as4m_control_idx = strcmp({AS4mData.Subject}, 'BCI02') & strcmp({AS4mData.Hand}, 'Ipsi');

stim_dir = fullfile(DataPath, "BioStims");
flist = dir(stim_dir);
target_amps = 40:10:80;
files = {flist.name};

stim_table = NaN(length(target_amps),2);
for t = 1:length(target_amps)
    lin_idx = find(contains(files, num2str(target_amps(t))) & contains(files, 'Lin'));
    temp = load(fullfile(stim_dir, files{lin_idx}));
    stim_table(t,1) = sum(temp(2,2:2:end));
    bio_idx = find(contains(files, num2str(target_amps(t))) & contains(files, 'Bio'));
    temp = load(fullfile(stim_dir, files{bio_idx}));
    stim_table(t,2) = sum(temp(2,2:2:end));
end
stim_table = stim_table / 100 / 5;

xq = linspace(1, 15);

% Extra detection thresholds
crs02b_dt = [2, 6.4;...
             14, 6;...
             18, 12.4;...
             44, 19.4;...
             51, 14.4]; % electrode, dt
crs07_dt = [2, 10.2;...
            12, 14.4;...
            16, 36.8;...
            46, 29;...
            48, 14.8;...
            62, 22.8;...
            17, 18.6;...
            24, 10.2;...
            39, 14]; % electrode, dt

[lin_jnd, bio_jnd, dt_values2, lin_jnd_nc, bio_jnd_nc] = deal(NaN(length(AS4bEData),1));
jnd_scatter_colors2 = zeros([length(AS4bEData), 3]);

for e = 1:length(AS4bEData)
    lin_jnd(e) = AS4bEData(e).JNDTable{1,1};
    bio_jnd(e) = AS4bEData(e).JNDTable{1,2};

    lin_jnd_nc(e) = AS4bEData(e).JNDTable{2,1};
    bio_jnd_nc(e) = AS4bEData(e).JNDTable{2,2};

    jnd_scatter_colors2(e,:) = SubjectColors(AS4bEData(e).ParticipantID);
    
    % Detection threshold
    if strcmp(AS4bEData(e).ParticipantID, 'BCI02')
        tidx = find([ProcessedDetectionData.Channel] == AS4bEData(e).Channel);
        dt_values2(e) = ProcessedDetectionData(tidx).MeanThreshold;
    elseif strcmp(AS4bEData(e).ParticipantID, 'CRS02b')
        tidx = find(crs02b_dt(:,1) == AS4bEData(e).Channel);
        dt_values2(e) = crs02b_dt(tidx,2);
    elseif strcmp(AS4bEData(e).ParticipantID, 'CRS07')
        tidx = find(crs07_dt(:,1) == AS4bEData(e).Channel);
        dt_values2(e) = crs07_dt(tidx,2);
    end
end

lin_jnd(lin_jnd > 40) = 40;
bio_jnd(bio_jnd > 40) = 40;
lin_jnd_nc(lin_jnd_nc > 1e4) = 1e4;
bio_jnd_nc(bio_jnd_nc > 1e4) = 1e4;

icms_bio_wf = bio_jnd ./ 60;
num_bio_levels = floor(log(100./dt_values2) ./ log(1 + icms_bio_wf));
icms_lin_wf = lin_jnd ./ 60;
num_lin_levels = floor(log(100./dt_values2) ./ log(1 + icms_lin_wf));

% Ratio of bio and lin intensities
rel_int = zeros(size(AS5bData));
for i = 1:length(AS5bData)
    y1 = mean(cat(1,AS5bData(i).FormattedTable.Bio{:}),2);
    y2 = mean(cat(1,AS5bData(i).FormattedTable.Lin{:}),2);

    rel_int(i) = mean(y1 ./ y2);
end

%% Supplement
if exist('fig', 'var') 
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax
end

mech_color = rgb(0, 137, 123);
control_color = [.6 .6 .6];
bio_color = rgb(66, 165, 245);
lin_color = [.6 .6 .6];

fig = figure;
set(fig, 'Units', 'inches', 'Position', [30, 5, 6.5, 4]);
h = 2; w = 3;
ax(1) = subplot(h,w,1); hold on
    Swarm(1, mech_jnd(as4m_test_idx), SubjectColors('BCI02'), "DS","box")
    Swarm(2, mech_jnd(as4m_control_idx), SubjectColors('BCI02'), "DS","box")

    set(ax(1), 'XTick', [1,2], ...
               'XTickLabels', {'Contra', 'Ipsi'}, ...
               'XLim', [.5 2.5], ...
               'YLim', [0 .15], ...
               'YTick', [0:.05:.15])

    ylabel(ax(1), 'JND (mm)')


ax(2) = subplot(h,w,2); hold on
    c = lines(size(NC3Processed,2));
    for i = 1:size(NC3Processed,2)
        x = NC3Processed(i).MechSummary.Amplitude;
        y = cat(2,NC3Processed(i).MechSummary.Responses{:});
        AlphaLine(x.*force_ratio,y',c(i,:), 'LineWidth', 2)
    end
    
    set(ax(2), 'YLim', [0 3], ...
               'YTick', [0:1:3])

    ylabel(ax(2), 'Normalized Intensity')
    xlabel(ax(2), 'Stimulus Force (N)')

ax(3) = subplot(h,w,3); hold on
    for i = 1:size(WeberData, 2)
        x = WeberData(i).Standards;
        y = [WeberData(i).ResponseSummary(1).JND, WeberData(i).ResponseSummary(2).JND];
        plot(x, y, 'Color', [.6 .6 .6], 'LineStyle',':')

        scatter(x(1), y(1), 50, SubjectColors('BCI02'), "filled", "o")
        scatter(x(2), y(2), 80, SubjectColors('BCI02'), "filled", "square")
    end
    ylabel('JND (uA)');
    xlabel('Standard (uA)')
    xlim([45, 75]); xticks([50, 70]);
    ylim([0, 30])

ax(4) = subplot(h,w,4); hold on
    [r,p] = corrcoef(dt_values1, icms_constant_jnd, 'Rows', 'complete');
    p1 = polyfit(dt_values1(~isnan(icms_constant_jnd)), icms_constant_jnd(~isnan(icms_constant_jnd)), 1);
    plot([0, max(dt_values1, [],'omitnan')], polyval(p1, [0, max(dt_values1, [],'omitnan')]), 'Color', [.6 .6 .6], 'LineStyle','--')
    scatter(dt_values1, icms_constant_jnd, 30, jnd_scatter_colors1, "filled",'MarkerEdgeColor','flat')
    leg_text = ColorText({'C1', 'P2', 'P3'},...
        [SubjectColors('BCI02'); SubjectColors('CRS02b'); SubjectColors('CRS07')]);
    text(50,30, leg_text, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
    xlabel(ax(4), sprintf('Detection Threshold (%sA)', GetUnicodeChar('mu')))
    ylabel(ax(4), sprintf('JND (%sA)', GetUnicodeChar('mu')))

    set(ax(4), 'XLim', [0 50],...
               'YLim', [0 30], ...
               'XTick', [0:20:60], ...
               'YTick', [0:10:30])

ax(5) = subplot(h,w,5); hold on
    plot([0, 3e3], [0, 3e3], 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(lin_jnd_nc, bio_jnd_nc, 50, jnd_scatter_colors2, "filled",'MarkerEdgeColor','flat',...
            'MarkerEdgeAlpha', 1)
    
    xlabel(ax(5), 'Linear JND (nC/p)')
    ylabel(ax(5), 'Biomimetic JND (nC/p)', 'Color', bio_color)
    leg_text = ColorText({'C1', 'P2', 'P3'},...
        [SubjectColors('BCI02'); SubjectColors('CRS02b'); SubjectColors('CRS07')]);
%     text(8 * 0.05, 8, leg_text, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

    set(ax(5), 'YLim', [0,8],...
           'YTick', [0:2:8],...
           'XLim', [0, 8], ...
           'XTick', [0:2:8], ...
           'XTickLabelRotation', 0, ...
           'Box', 'off')

ax(6) = subplot(h,w,6); hold on
i = 1;
    x = 40:10:80;
    y = cat(1,AS5bData(i).FormattedTable.Bio{:});
    AlphaLine(x,y,bio_color)
    y = cat(1,AS5bData(i).FormattedTable.Lin{:});
    AlphaLine(x,y,lin_color)
    
    set(ax(6), 'XLim', [38, 82], ...
               'XTick', [x])

    xlabel(ax(6), sprintf('Stimulus Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel(ax(6), 'Normalized Intensity')


AddFigureLabels(gcf, [0.05, -0.025]);

export_path = fullfile(DataPath, 'SuppFig6_Mech');
print(gcf, export_path, '-dsvg', '-r300')