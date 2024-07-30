%% Imports
addpath("HelperFunctions")
load(fullfile(DataPath, "AmpDiscrimBioLin"))
load(fullfile(DataPath, "ForceEquivMulti.mat")) % Mech amp est 
load(fullfile(DataPath, "AmpDiscrimMultiChanRaw")) % multi chan discrim
load(fullfile(DataPath, "AmpDiscrimMultiBioLin.mat")) % multi chan discrim bio vs lin
load(fullfile(DataPath, "MultiChanDTs.mat")) % multi chan discrim
load(fullfile(DataPath, "DetectionDataProcessed_BCI02.mat")); % Detection data
% Figure options
SetFont('Arial', 9)

scatter_size = 30;
line_width = 1.5;
subject_ids = {'BCI02', 'CRS02', 'CRS07'};

force_ratio = 0.37;
invsig = @(c,y) (log(1/y - 1)/-c(1));
force2amp = @(pow_coeffs, lin_coeffs, x) (pow_coeffs(1) * (x + pow_coeffs(3)).^pow_coeffs(2) - lin_coeffs(1)) ./ lin_coeffs(2);
amp2force = @(pow_coeffs, lin_coeffs, x) (((lin_coeffs(2) * x + lin_coeffs(1)) ./ pow_coeffs(1)) .^ (1/pow_coeffs(2))) - pow_coeffs(3);


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
            17, 18.6
            24, 10.2;...
            39, 14]; % electrode, dt

%% Data extraction
[lin_jnd, bio_jnd, dt_values, lin_jnd_nc, bio_jnd_nc] = deal(NaN(length(AS4bEData),1));

for e = 1:length(AS4bEData)
    lin_jnd(e) = AS4bEData(e).JNDTable{1,1};
    bio_jnd(e) = AS4bEData(e).JNDTable{1,2};

    lin_jnd_nc(e) = AS4bEData(e).JNDTable{2,1};
    bio_jnd_nc(e) = AS4bEData(e).JNDTable{2,2};

    
    % Detection threshold
    if strcmp(AS4bEData(e).ParticipantID, 'BCI02')
        tidx = find([ProcessedDetectionData.Channel] == AS4bEData(e).Channel);
        dt_values(e) = ProcessedDetectionData(tidx).MeanThreshold;
    elseif strcmp(AS4bEData(e).ParticipantID, 'CRS02b')
        tidx = find(crs02b_dt(:,1) == AS4bEData(e).Channel);
        dt_values(e) = crs02b_dt(tidx,2);
    elseif strcmp(AS4bEData(e).ParticipantID, 'CRS07')
        tidx = find(crs07_dt(:,1) == AS4bEData(e).Channel);
        dt_values(e) = crs07_dt(tidx,2);
    end
end

mc_jnds = NaN(length(AS4mcEData),1);
for i = 1:length(AS4mcEData)
    mc_jnds(i) = AS4mcEData(i).JNDTable{"uA","pCH"};
end
mc_wf = mc_jnds ./ 60;

multi_mc_jnds = mc_jnds([AS4mcEData.IsMulti]);
mc_idx = find([AS4mcEData.IsMulti]);
mc_levels = NaN(size(multi_mc_jnds));
for i = 1:length(multi_mc_jnds)
    soi = AS4mcEData(mc_idx(i)).ParticipantID;
    coi = AS4mcEData(mc_idx(i)).Channel;
    % Find detection threshold
    dt_idx = contains({MCDTStruct.SubjectID}, soi) & cellfun(@(c) all(ismember(coi, c)), {MCDTStruct.Channel});
%     mc_levels(i) = floor((100 - MCDTStruct(dt_idx).DT50) / multi_mc_jnds(i));
    y = MCDTStruct(dt_idx).DT50;
    ndl = 0;
    while y < 100
        y = y * (1+mc_wf(i));
        ndl = ndl + 1;
    end
    mc_levels(i) = ndl - 1;
end

% Compute num levels
icms_lin_wf = lin_jnd ./ 60;
icms_bio_wf = bio_jnd ./ 60;
[num_lin_levels, num_bio_levels] = deal(NaN(size(lin_jnd)));
for i = 1:length(lin_jnd)
    % Lin
    y = dt_values(i);
    ndl = 0;
    while y < 100
        y = y * (1+icms_lin_wf(i));
        ndl = ndl + 1;
    end
    num_lin_levels(i) = ndl - 1;

    % Bio
    y = dt_values(i);
    ndl = 0;
    while y < 100
        y = y * (1+icms_bio_wf(i));
        ndl = ndl + 1;
    end
    num_bio_levels(i) = ndl - 1;
end

[P,H,STATS] = ranksum(mc_jnds([AS4mcEData.IsMulti]), mc_jnds(~[AS4mcEData.IsMulti]));
fprintf('Single vs Multi-Channel JND (Ranksum test):\n Z = %0.2f\n %s\n', STATS.zval, pStr(P))
fprintf('Median single JND = %0.2f\nMedian multi JND = %0.2f\n', median(mc_jnds(~[AS4mcEData.IsMulti])), median(mc_jnds([AS4mcEData.IsMulti])))

[mc_bio_jnd, mc_lin_jnd] = deal(zeros(length(AS4_Multi_LinBio),1));
for i = 1:length(AS4_Multi_LinBio)
    mc_bio_jnd(i) = AS4_Multi_LinBio(i).Sigmoids{"BIO", "JND"};
    mc_lin_jnd(i) = AS4_Multi_LinBio(i).Sigmoids{"LIN", "JND"};
end

%% Main
if exist('fig', 'var') 
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax
end

mech_color = rgb(0, 137, 123);
multich_color = rgb(244, 67, 54);
bio_color = rgb(66, 165, 245);
control_color = [.6 .6 .6];

fig = figure;
set(fig, 'Units', 'inches', 'Position', [30, 5, 6.25, 8.5]);
ax(1) = axes('Position', [0.075 0.75 0.385 0.225]); hold on
i = 1;
    % Plot mech
    y = horzcat(NC3Processed(i).MechSummary.Responses{:})';
    AlphaLine(linspace(0,1,size(NC3Processed(i).MechSummary,1)),y, mech_color, 'LineWidth', 1.5)
    % Plot ICMS
    for j = 1:size(NC3Processed(i).ICMSSummary.SingleChResponses,2)
        y = horzcat(NC3Processed(i).ICMSSummary.SingleChResponses{:,j})';
        AlphaLine(linspace(0,1,size(NC3Processed(i).ICMSSummary,1)),y, bio_color, 'LineWidth', 1.5)
    end
    y = horzcat(NC3Processed(i).ICMSSummary.MultiChResponses{:})';
    AlphaLine(linspace(0,1,size(NC3Processed(i).ICMSSummary.MultiChResponses,1)),y, multich_color, 'LineWidth', 1.5)

    % Format
    yl = get(gca, 'YLim');
    yl = [floor(yl(1)*10)/10, ceil(yl(2)*10)/10];
    ylabel(ax(1), 'Normalized Intensity')
    set(ax(1), 'XColor', 'none',...
               'XLim', [-0.05 1.05],...
               'YLim', yl,...
               'Clipping', 'off')

    text(.025,yl(1)+range(yl), ColorText(fliplr({'Mechanical', 'Bio Multi-Channel', 'Bio Single-Channel'}), ...
        flipud([mech_color; multich_color; bio_color])), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment','top')

    % Mech x-axis
    plot([0, 1], [yl(1), yl(1)], 'Color', mech_color, 'LineWidth', 1)
    x = NC3Processed(i).MechSummary.Amplitude * force_ratio;
    xs = linspace(0,1,length(x));
    for xsi = [1:length(xs)]
        plot([xs(xsi), xs(xsi)], [yl(1), yl(1)-range(yl)*0.025], 'Color', mech_color, 'LineWidth', 1)
        if xsi == 1
            text(xs(xsi), yl(1)-range(yl)*0.025, sprintf('%0.2f', x(xsi)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', 'Color', mech_color)
        elseif xsi == length(xs)
            text(xs(xsi), yl(1)-range(yl)*0.025, sprintf('%0.2f N', x(xsi)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', 'Color', mech_color)
        end
    end

    % ICMS x-axis
    plot([0, 1], [yl(1), yl(1)]-range(yl)*0.125, 'Color', multich_color, 'LineWidth', 1)
    x = NC3Processed(i).ICMSSummary.Amplitude;
    xs = linspace(0,1,length(x));
    for xsi = [1:length(xs)]
        plot([xs(xsi), xs(xsi)], [yl(1), yl(1)-range(yl)*0.025]-range(yl)*0.125, 'Color', multich_color, 'LineWidth', 1)
        if xsi == 1
            text(xs(xsi), yl(1)-range(yl)*0.025-range(yl)*0.125, num2str(x(xsi)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', 'Color', multich_color)
        elseif xsi == length(xs)
            text(xs(xsi), yl(1)-range(yl)*0.025-range(yl)*0.125, sprintf('%d %sA', x(xsi), GetUnicodeChar('mu')), ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', 'Color', multich_color)
        end
    end

ax(2) = axes('Position', [0.56 0.74 0.385 0.235]); hold on
    c = lines(length(NC3Processed));
    % Plot single channel
    x = linspace(0, max(NC3Processed(i).MechRange(2))*2);
    for j = 1:length(NC3Processed(i).ICMSFitSingle)
        y = force2amp(coeffvalues(NC3Processed(i).MechFit.PowFun), coeffvalues(NC3Processed(i).ICMSFitSingle(j).LinFun), x);
        plot(x.* force_ratio, y, 'Color', bio_color, 'LineWidth', line_width);
    end
    % Plot multi channel
    y = force2amp(coeffvalues(NC3Processed(i).MechFit.PowFun), coeffvalues(NC3Processed(i).ICMSFitMulti.LinFun), x);
    plot(x.* force_ratio, y, 'Color', multich_color, 'LineWidth', 1.5);


    set(ax(2), 'YLim', [0 80], 'XLim', [0 1])
    ylabel(ax(2), sprintf('ICMS Amplitude (%sA)', GetUnicodeChar('mu')))
    xlabel(ax(2), 'Stimulus Force (N)')

inset_ax(1) = axes('Position', [0.7 0.775 0.245 0.085]); hold on
    [single_max, multi_max] = deal(cell(size(NC3Processed)));
    for j = 1:length(NC3Processed)
        single_max{j} = [NC3Processed(j).SingleChMaxIndent];
        multi_max{j} = NC3Processed(j).MultiChMaxIndent;
    end
    single_max = unique(horzcat(single_max{:})) .* force_ratio;
    multi_max = cell2mat(multi_max) .* force_ratio;
    eq = linspace(0,2,12);
    [counts, edges, bins]= histcounts(single_max, eq, 'Normalization','probability');
    for i = 1:length(edges)-1
        patch([edges(i), edges(i), edges(i+1), edges(i+1)], [0 counts(i), counts(i), 0], ...
            bio_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none')
        plot([edges(i), edges(i), edges(i+1), edges(i+1)], [0 counts(i), counts(i), 0], ...
            'Color', [.2 .2 .2], 'LineWidth', 0.5)
    end
    [counts, edges, bins]= histcounts(multi_max, eq, 'Normalization','probability');
    for i = 1:length(edges)-1
        patch([edges(i), edges(i), edges(i+1), edges(i+1)], [0 counts(i), counts(i), 0], ...
            multich_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none')
        plot([edges(i), edges(i), edges(i+1), edges(i+1)], [0 counts(i), counts(i), 0], ...
            'Color', [.2 .2 .2], 'LineWidth', 0.5)
    end

    text(median(single_max, 'omitnan'), 0.5, char(8964), 'Color', bio_color, 'FontSize', 20, ...
        'FontWeight', 'bold', 'FontName', 'JuliaMono', 'VerticalAlignment','bottom', 'HorizontalAlignment','center')
    text(median(multi_max, 'omitnan'), 0.5, char(8964), 'Color', multich_color, 'FontSize', 20,...
        'FontWeight', 'bold', 'FontName', 'JuliaMono', 'VerticalAlignment','bottom', 'HorizontalAlignment','center')
    set(inset_ax(1), 'XTick', [0 2], 'YTick', [0 .5], 'YTickLabel', {'0', '.5'})
    ylabel(inset_ax(1), 'Proportion')
    xlabel(inset_ax(1), 'Max Force (N)', 'VerticalAlignment','bottom')

i = 15;
ax(3) = axes('Position', [0.075 0.45 0.3 0.2]); hold on
    sigfun = GetSigmoid(2);
    xq = linspace(40,80);
    single_channels = AS4mcEData(i).Channel;
    if length(single_channels) < 4
        error('Multi-channel is not selected')
    end
    sc_idx = find(~[AS4mcEData.IsMulti]);
    sc = [AS4mcEData(sc_idx).Channel];
    % Plot each single channel
    for j = 1:length(single_channels)
        j_idx = find(sc == single_channels(j));
        scatter(AS4mcEData(sc_idx(j_idx)).ResponseSummary.Amp, AS4mcEData(sc_idx(j_idx)).ResponseSummary.pCH, 30, 'MarkerFaceColor', ...
           bio_color, 'MarkerEdgeColor', bio_color, 'MarkerFaceAlpha', 0.5)
        plot(xq, sigfun(AS4mcEData(sc_idx(j_idx)).Sigmoids{"uA", "pCH"}{1}, xq), "Color", bio_color, 'LineWidth', 1.5)
    end
    % Plot multi channel
    scatter(AS4mcEData(i).ResponseSummary.Amp, AS4mcEData(i).ResponseSummary.pCH, 30, 'MarkerFaceColor', rgb(244, 67, 54), ...
        'MarkerEdgeColor', rgb(244, 67, 54))
    plot(xq, sigfun(AS4mcEData(i).Sigmoids{"uA", "pCH"}{1}, xq), "Color", multich_color, 'LineWidth', 1.5)
    
    xlabel(ax(3), sprintf('Comparison Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel(ax(3), 'p(Comp Higher)')
    set(ax(3), 'XLim', [38, 82],...
               'YLim', [0,1],...
               'XTick', [40:10:80], ...
               'YTick', [0 1])

inset_ax(2) = axes('Position', [0.125 0.56 0.075 0.09]); hold on
    Swarm(1, mc_jnds(~[AS4mcEData.IsMulti]), bio_color, "DS", "box", 'SYL', [0 15])
    Swarm(2,mc_jnds([AS4mcEData.IsMulti]), multich_color, "DS","box")
    set(inset_ax(2), 'XTick', [], 'XLim', [.5 2.5], 'YTick', [0 15], 'YLim', [0 15])
    ylabel(inset_ax(2), sprintf('JND (%sA)', GetUnicodeChar('mu')), 'VerticalAlignment','top', 'FontSize', 8)

ax(4) = axes('Position', [0.425 0.45 0.3 0.2]); hold on
    sigfun = GetSigmoid(1);
    i = 2;
    idx = strcmp(AS4_Multi_LinBio(i).ConditionMeanTable.Type, "LIN");
    x = str2double(AS4_Multi_LinBio(i).ConditionMeanTable.CompAmp(idx));
    y = str2double(AS4_Multi_LinBio(i).ConditionMeanTable.CompHigher(idx));
    scatter(x,y,30, 'MarkerFaceColor', ...
           rgb(66, 66, 66), 'MarkerEdgeColor', rgb(66, 66, 66), 'MarkerFaceAlpha', 0.5)
    xq = linspace(x(1), x(end));
    yq = sigfun(AS4_Multi_LinBio(1).Sigmoids{"BIO", "Coeffs"}{1}, xq - 60);
    plot(xq,yq, 'Color', rgb(69, 39, 160), 'LineWidth', 1.5)
    
    idx = strcmp(AS4_Multi_LinBio(i).ConditionMeanTable.Type, "BIO");
    x = str2double(AS4_Multi_LinBio(i).ConditionMeanTable.CompAmp(idx));
    y = str2double(AS4_Multi_LinBio(i).ConditionMeanTable.CompHigher(idx));
    scatter(x,y,30, 'MarkerFaceColor', rgb(244, 67, 54), ...
        'MarkerEdgeColor', rgb(244, 67, 54))
    yq = sigfun(AS4_Multi_LinBio(1).Sigmoids{"LIN", "Coeffs"}{1}, xq - 60);
    plot(xq,yq, 'Color', multich_color, 'LineWidth', 1.5)

    xlabel(ax(4), sprintf('Comparison Amplitude (%sA)', GetUnicodeChar('mu')))
    set(ax(4), 'XLim', [38, 82],...
               'YLim', [0,1],...
               'XTick', [40:10:80], ...
               'YTick', [0 1])
    [x,y] = GetAxisPosition(ax(4), 95, 5);
    text(x,y, ColorText({'Lin Multi-Channel', 'Bio Multi-Channel'}, [rgb(66, 66, 66); rgb(244, 67, 54)]), ...
        'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'right')

inset_ax(3) = axes('Position', [0.475 0.56 0.075 0.09]); hold on
    Swarm(1, mc_lin_jnd, rgb(66, 66, 66), "DS", "box")
    Swarm(2, mc_bio_jnd, multich_color, "DS","box")
    set(inset_ax(3), 'XTick', [], 'XLim', [.5 2.5], 'YTick', [5 15], 'YLim', [5 15])
    ylabel(inset_ax(3), sprintf('JND (%sA)', GetUnicodeChar('mu')), 'VerticalAlignment','top', 'FontSize', 8)

ax(5) = axes('Position', [0.8 0.45 0.15 0.2]); hold on
    Swarm(1, num_lin_levels, [.6 .6 .6], 'DS', 'Box')
    Swarm(2, num_bio_levels, rgb(66, 165, 245), 'DS', 'Box')
    Swarm(3, mc_levels, multich_color, 'DS', 'Box')

    ylabel(ax(5), '# Discriminable Levels')
    set(ax(5), 'YLim', [0, 45], ...
               'XLim', [.5 3.5], ...
               'YTick', [0:20:40],...
               'XTick', [1:3], ...
               'XTickLabel', ColorText({'Linear', 'BioSingle', 'BioMulti'}, [[.6 .6 .6]; rgb(66, 165, 245); multich_color]))


ax(6) = axes('Position', [0.05 0.0 0.7 0.35]); hold on
imshow(".\ReferenceImages\Force_scheme.png")
    set(gca, 'XTick', [], 'YTick', [])

ax(7) = axes('Position', [0.8 0.09 0.15 0.25]); hold on
    bio_perf = [18, 19]; % Taken from test log, blocks of 20
    lin_perf = [14, 16];
    Swarm(1,  1 - lin_perf ./ 20, [.6 .6 .6], 'DS', 'Bar', 'EW', false)
    Swarm(2, 1 - bio_perf ./ 20, multich_color, 'DS', 'Bar', 'EW', false)
    set(gca, 'XTick', [1,2], ...
             'XTickLabel', {'Lin Single', 'Bio Multi'}, ...
             'XLim', [0.5 2.5], ...
             'YTick', [0:.1:.3])
    ylabel('Error Rate')

%%% Labels
char_offset = 64;
AddFigureLabels(ax(1:5), [0.05, 0.02]);
annotation("textbox", [0.025 0.32 .05 .05], 'String', char(char_offset+6), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.75 0.32 .05 .05], 'String', char(char_offset+7), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')

export_path = fullfile(DataPath, 'Figure7_MultiChannel');
print(fig, export_path, '-dsvg', '-r300')
