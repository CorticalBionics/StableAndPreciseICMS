% Imports
addpath("HelperFunctions")
load(fullfile(DataPath, "WeberData.mat")) % Weber fraction
load(fullfile(DataPath, "AmpDiscrim.mat")) % ICMS amp discrim
load(fullfile(DataPath, "ForceEquiv_UC.mat")) % Mech amp est
load(fullfile(DataPath, "AmpDiscrimMultiChan.mat")) % Multi chan discrim
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
% Weber fractions
weber_fracs = zeros(length(WeberData),2);
for i = 1:length(WeberData)
    weber_fracs(i,:) = [WeberData(i).ResponseSummary.JND] ./ [WeberData(i).Standards]';
end

% Detection data given by Pitt at date of data collection
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

jnd_thresh = 40;
jnd_idx = icms_constant_jnd <= jnd_thresh;
fprintf('%d of %d (%0.2f%%) Channels with JND <= 40\n',...
    sum(jnd_idx), length(jnd_idx), sum(jnd_idx) / length(jnd_idx)*100)
fprintf('Median JND = %0.2f, mean JND = %0.2f\n', median(icms_constant_jnd, 'omitnan'), mean(icms_constant_jnd, 'omitnan'))

% Compute discriminable levels
icms_constant_wf = icms_constant_jnd ./ 60;
num_disc_levels = floor(log(100./dt_values1) ./ log(1 + icms_constant_wf));

% Get max force values
[dt_subset, mf_vals, mf_jnds] = deal(NaN(size(NC1Processed,2),1));
for i = 1:size(NC1Processed,2)
    tidx = find([ProcessedDetectionData.Channel] == NC1Processed(i).Channel);
    dt_subset(i) = ProcessedDetectionData(tidx).MeanThreshold;
    mf_vals(i) = NC1Processed(i).MaxIndent;
    % Find the JND
    for k = 1:size(MC_AmpDiscrim,2)
        if ismember(NC1Processed(i).Channel, MC_AmpDiscrim(k).Channels)
            k_idx = find(NC1Processed(i).Channel == MC_AmpDiscrim(k).Channels);
            mf_jnds(i) = MC_AmpDiscrim(k).SC_JND(k_idx);
            break
        end
    end
end
mf_disc_levels = (100 - dt_subset) ./ mf_jnds;

%% Bio extraction
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


%% Main
if exist('fig', 'var') 
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax
end

% Choose example electrodes
example_idx = [20,27,34];

fig = figure('Units', 'inches', 'Position', [5, 5, 4.5, 5.5]);

%%% Amplitude discrim psychometric
ax(1) = axes('Position', [.1 .75 .25 .225]); hold on
    for i = example_idx
        isi_idx = AS4EData(i).ResponseSummary.ISI == 1;
        x = AS4EData(i).ResponseSummary.CompAmp(isi_idx);
        % Mean
        y = AS4EData(i).ResponseSummary.pH(isi_idx);
        scatter(x,y,scatter_size,'MarkerFaceColor', SubjectColors(AS4EData(i).ParticipantID),...
            'MarkerEdgeColor', SubjectColors(AS4EData(i).ParticipantID))
        isi_idx = AS4EData(i).ISISigmoids.ISI == 1;
        c = AS4EData(i).ISISigmoids.SigCoeffs(isi_idx,:);
        plot(linspace(min(x), max(x)), sigfun(c, linspace(min(x), max(x))),...
            'Color', SubjectColors(AS4EData(i).ParticipantID), 'LineWidth', line_width)
    end

    leg_text = ColorText({'C1', 'P2', 'P3'},...
        [SubjectColors('BCI02'); SubjectColors('CRS02b'); SubjectColors('CRS07')]);
    text(40, 1, leg_text, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

    set(ax(1), 'YLim', [0,1],...
               'YTick', [0,1],...
               'XLim', [38, 82], ...
               'XTick', [40:20:80], ...
               'XTickLabelRotation', 0)
    xlabel(ax(1), sprintf('Comparison Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel(ax(1), 'p(Comp Higher)')

%%% Amplitude discrimination JNDs
ax(2) = axes('Position', [.45 .75 .125 .225]); hold on
    Swarm(1, icms_constant_jnd, 'SwarmColor', jnd_scatter_colors1, "DistributionStyle", "Box",...
        'CenterMethod','median', 'Parent', ax(2), 'DistributionFaceAlpha', 0.05, 'SwarmYLimits', [0 40])
    set(ax(2), 'XTick', [], ...
              'YTick', [-40:20:40],...
              'YLim', [0 40],...
              'XLim', [.6 1.4])
    ylabel(ax(2), sprintf('JND(%sA)', GetUnicodeChar('mu')));

%%% Weber fraction
ax(3) = axes('Position', [.675 .75 .25 .225]); hold on
    plot([0 0.4], [0 0.4], 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(weber_fracs(:,1), weber_fracs(:,2), scatter_size, SubjectColors('BCI02'), 'filled')
    xlabel(sprintf('%sI/I @ 50 %sA', GetUnicodeChar('Delta'), GetUnicodeChar('mu')))
    ylabel(sprintf('%sI/I @ 70 %sA', GetUnicodeChar('Delta'), GetUnicodeChar('mu')))
    set(ax(3), 'XTick', [0:0.2:0.4], ...
               'YTick', [0:0.2:0.4],...
               'YLim', [0 0.4],...
               'XLim', [0 0.4])

%%% Discriminable levels example
ax(4) = axes('Position', [.1 .4125 .3 .225]); hold on
    ylabel(ax(4),sprintf('Amplitude (%sA)', GetUnicodeChar('mu')), 'VerticalAlignment','middle')
    ii = 1;
    for i = example_idx
        dt = dt_values1(i);
        wf = icms_constant_wf(i);
        num_levels = num_disc_levels(i); % Remaining levels (assuming max amp == 100)
        % Plot the levels
        x = [ii,ii,ii+1,ii+1];
        patch(x, [0 dt, dt,0], [.6 .6 .6], 'FaceAlpha', 0.3) % Detection threshold
        c = flipud(summer(num_levels));
        yl = dt;
        yu = dt * (1+wf);
        for j = 1:num_levels
            patch(x, [yl yu, yu,yl], c(j,:), 'FaceAlpha', 1)
            yl = yu;
            yu = yu * (1+wf);
        end
        
        text(ii + 0.5, -5, ColorText({sprintf('Ch%d',AS4EData(i).Channel)}, SubjectColors(AS4EData(i).ParticipantID)),...
            'VerticalAlignment','top', 'HorizontalAlignment','center')
    
        ii = ii + 1.25;
    end
    
    set(ax(4), 'YLim', [0, 100], ...
               'XTick', [], ...
               'XLim', [0.6 5.5], ...
               'XColor', 'none')


%%% Mechanical discriminable levels
ax(5) = axes('Position', [.45 .4125 .1 .225]); hold on    
    weber = 1.0829;
    dt = 0.025 * force_ratio;
    x = [0,0,1,1];
    patch(x, [0.00 dt, dt, 0.00], [.6 .6 .6], 'LineWidth', .1)
    num_levels = 47;
    c = flipud(summer(num_levels));
    yl = dt;
    yu = dt * weber;
    for i = 1:num_levels
        patch(x, [yl yu, yu,yl], c(i,:), 'LineWidth', .1, 'FaceAlpha', 1)
        yl = yu;
        yu = yu * weber;
    end

    ylabel(ax(5),'Force (N)', 'VerticalAlignment','top')
    set(ax(5), 'YScale', 'linear', ...
               'YLim', [0, .4],...
               'YTick', [0, .4], ...
               'XTick', [], ...
               'XColor', 'none',...
               'XLim', [-.3 1.1])

    %%% Weber Fractions
ax(6) = axes('Position', [.675 .4125 .25 .225]); hold on    
    [r,p] = corrcoef(mf_vals, mf_disc_levels, 'Rows','complete');
    p1 = polyfit(mf_vals(~isnan(mf_disc_levels)), mf_disc_levels(~isnan(mf_disc_levels)),1);
    plot([0, max(mf_vals, [],'omitnan')], polyval(p1, [0, max(mf_vals, [],'omitnan')]), 'Color', [.6 .6 .6], 'LineStyle','--')
    scatter(mf_vals, mf_disc_levels, 30, 'MarkerEdgeColor',SubjectColors('BCI02'), 'MarkerFaceColor',SubjectColors('BCI02'))
    xlabel('Max Force (N)'); ylabel('# Discriminable Levels')
    set(ax(6), 'YLim', [4 9],...
               'YTick', [4:9], ...
               'XTick', [0.3:.3:0.9], ...
               'XLim', [0.3, 0.9])

amps = [80:-10:40];
lin_colors = ColorGradient(rgb(33, 33, 33), rgb(189, 189, 189), 5);
bio_colors = ColorGradient(rgb(13, 71, 161), rgb(66, 165, 245), 5);

ax(7) = axes('Position', [.025 .05 .25 .25]); hold on 
    for i = 1:5
        if i == 3
            ls = '-';
        else
            ls = ':';
        end
        % Linear
        plot([-.1 0 .1 .9 1 1.1], [0 0 amps(end-i+1) amps(end-i+1) 0 0], 'Color', lin_colors(end-i+1,:), 'Parent', ax(7), ...
            'LineWidth', line_width, 'LineStyle', ls)

        plot([-.1 0 0 .1 .1 .9 .9 1 1 1.1], ...
            [0 0 amps(end-i+1) amps(end-i+1) amps(end-i+1)/2 amps(end-i+1)/2 amps(end-i+1) amps(end-i+1) 0 0 ] - 100, ...
            'Color', bio_colors(end-i+1,:), 'Parent', ax(7),...
            'LineWidth', line_width, 'LineStyle', ls)
    end
    % Axes
    plot([-.1, -.1, .15]-.05, [5, -5, -5], 'Color', [.6 .6 .6], 'Parent', ax(1))
    text(0,-5, '0.25s', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', [.6 .6 .6],...
        'Parent', ax(1))
    text(-.175, -0.05, sprintf('10 %sA', GetUnicodeChar('mu')), 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom', 'Color', [.6 .6 .6], 'Parent', ax(1), 'Rotation', 90)
    text(.5, 20, 'Linear', 'Color', 'k', 'HorizontalAlignment','center', 'VerticalAlignment','top')
    text(.5, -105, 'Biomimetic', 'Color', rgb(66, 165, 245), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
    set(ax(7), 'XColor', 'none',...
               'YColor', 'None', ...
               'YLim', [-100 80],...
               'XLim', [-.2 1.2])

ax(8) = axes('Position', [.4 .1 .23 .2]); hold on 
    plot([0, 38], [0, 38], 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(lin_jnd, bio_jnd, 50, jnd_scatter_colors2, "filled",'MarkerEdgeColor','flat', 'MarkerEdgeAlpha', 1)
    
    xlabel(ax(8), sprintf('Linear JND (%sA)', GetUnicodeChar('mu')))
    ylabel(ax(8), sprintf('Bio JND (%sA)', GetUnicodeChar('mu')), 'Color', rgb(66, 165, 245))

    set(ax(8), 'YLim', [0,40],...
           'YTick', [0:10:40],...
           'XLim', [0, 40], ...
           'XTick', [0:10:40], ...
           'XTickLabelRotation', 0)

ax(9) = axes('Position', [.76 .1 .175 .2]); hold on 
    Swarm(1, num_disc_levels, [.6 .6 .6], 'DS', 'Box', 'SPL', 0)
    Swarm(2, num_lin_levels, [.2 .2 .2], 'DS', 'Box', 'SPL', 0)
    Swarm(3, num_bio_levels, rgb(66, 165, 245), 'DS', 'Box', 'SPL', 0)
    set(ax(9), 'XLim', [0.5 3.5], ...
               'YLim', [0 25], ...
               'YTick', [0:10:20], ...
               'XTick', [1:3], ...
               'XTickLabel', {'Flat', 'Lin', 'Bio'})
    ylabel('# Disc. Levels')

 
AddFigureLabels(gcf, [0.1, 0.02]);
shg
