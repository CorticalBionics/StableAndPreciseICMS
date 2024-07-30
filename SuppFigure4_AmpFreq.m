% Comparison of freq and amp on RF size
addpath("HelperFunctions")
load(fullfile(DataPath, 'IntraDayVariability_Processed.mat'))
load(fullfile(DataPath, 'AmpFreqSurvey_Processed.mat'))

AmpFreqResults = struct();

%% Remove CRS02 from IntraDay structure
IntraDayVariability = IntraDayVariability(strcmp({IntraDayVariability.Subject}, 'BCI02'));
% Remove num trials from Intraday?
for i = 1:length(IntraDayVariability)
    IntraDayVariability(i).SensoryField = rmfield(IntraDayVariability(i).SensoryField, 'NumTrials');
end
IntraDayVariability = rmfield(IntraDayVariability, 'Stability');
%% Concatenate the structs to make analysis easy
CombinedStruct = struct();
num_ch = length(AmpFreqSurvey);
for ch = 1:num_ch
    CombinedStruct(ch).Subject = IntraDayVariability(ch).Subject;
    CombinedStruct(ch).Channel = IntraDayVariability(ch).Channel;
    temp = rmfield(AmpFreqSurvey(ch).SensoryField, 'NumTrials');
    CombinedStruct(ch).SensoryField = [IntraDayVariability(ch).SensoryField, temp];
end

ch_colors = lines(num_ch);
%% Make a table for the intensities & size and normalize by block
table_names_types = [["Channel", "double"]; ...
                     ["Session", "double"]; ...
                     ["Set", "double"]; ...
                     ["Amplitude", "double"]; ...
                     ["Frequency", "double"]; ...
                     ["Intensity", "double"]; ...
                     ["Size", "double"]; ...
                     ["DigitVec", "cell"]];
ResponseTable = table('Size',[num_ch*30,size(table_names_types,1)],... 
     'VariableNames', table_names_types(:,1),...
     'VariableTypes', table_names_types(:,2));

ii = 1;
for ch = 1:num_ch
    for i = 1:length(AmpFreqSurvey(ch).SensoryField)
        ResponseTable{ii, 'Channel'} = CombinedStruct(ch).Channel;
        ResponseTable{ii, 'Session'} = CombinedStruct(ch).SensoryField(i).Session;
        ResponseTable{ii, 'Set'} = CombinedStruct(ch).SensoryField(i).Set;
        ResponseTable{ii, 'Amplitude'} = CombinedStruct(ch).SensoryField(i).Amplitude;
        ResponseTable{ii, 'Frequency'} = CombinedStruct(ch).SensoryField(i).Frequency;
        ResponseTable{ii, 'Intensity'} = CombinedStruct(ch).SensoryField(i).Intensity;
        ResponseTable{ii, 'Size'} = sum([CombinedStruct(ch).SensoryField(i).PalmSF.Area, CombinedStruct(ch).SensoryField(i).DorsSF.Area]);
        digit_vec = false(5,1);
        for dvi = 1:5
            ds = sprintf('D%d', dvi);
            if any(contains([{CombinedStruct(ch).SensoryField(i).PalmSF.CentroidLocation},...
                             {CombinedStruct(ch).SensoryField(i).DorsSF.CentroidLocation}], ds))
                digit_vec(dvi) = true;
            end
        end
        ResponseTable{ii, 'DigitVec'} = {digit_vec'};
        ii = ii + 1;
    end
end
ResponseTable = ResponseTable(1:ii-1,:); % Remove empty rows

% Normalize by session/set/channel
norm_vec = NaN(ii-1, 1);
u_ses = table2array(unique(ResponseTable(:,"Session")));
u_ch = table2array(unique(ResponseTable(:,"Channel")));
% One block has uneven number of trials per channel so formatting needs to
% be done awkwardly
for s = 1:size(u_ses,1)
    % Organize intensities wrt block and channel
    num_trials_in_set = sum(ResponseTable{:,["Session"]} == u_ses(s));
    blocks_in_set = ceil(num_trials_in_set / num_ch);
    ch_ints = NaN(blocks_in_set, num_ch);
    for ch = 1:num_ch
        idx = sort(find(ResponseTable{:,"Session"} == u_ses(s) & ...
                        ResponseTable{:,"Channel"} == u_ch(ch)));
        ch_ints(1:length(idx),ch) = ResponseTable{idx,"Intensity"};
    end
    ch_ints_norm = ch_ints ./ mean(ch_ints,2,'omitnan');

    % Assign normalized values back to vector
    for ch = 1:num_ch
        idx = sort(find(ResponseTable{:,"Session"} == u_ses(s) & ...
                        ResponseTable{:,"Channel"} == u_ch(ch)));
        norm_vec(idx) = ch_ints_norm(1:length(idx),ch);
    end
end
ResponseTable.IntensityNorm = norm_vec;

%% Average across trials
% Because we have trials per condition we cannot do a repeated measures
% ANOVA, thus we'll take the mean for ChxAmpxFreq and do a 3-way ANOVA
table_names_types = [["Channel", "double"]; ...
                     ["Amplitude", "double"]; ...
                     ["Frequency", "double"]; ...
                     ["Intensity", "double"]
                     ["IntensityNorm", "double"];...
                     ["Size", "double"];...
                     ["SizeNorm", "double"];...
                     ["NumDigits", "double"]];
MeanResponseTable = table('Size',[num_ch*5,size(table_names_types,1)],... 
     'VariableNames', table_names_types(:,1),...
     'VariableTypes', table_names_types(:,2));

ii = 1;
for ch = 1:num_ch
    ch_idx = ResponseTable{:,"Channel"} == u_ch(ch);
    u_amp_freq = unique(ResponseTable{ch_idx,["Amplitude", "Frequency"]}, 'rows');
    for uaf = 1:size(u_amp_freq,1)
        uaf_idx = ch_idx & ...
                  ResponseTable{:,"Amplitude"} == u_amp_freq(uaf,1) & ...
                  ResponseTable{:,"Frequency"} == u_amp_freq(uaf,2);

        MeanResponseTable{ii, 'Channel'} = u_ch(ch);
        MeanResponseTable{ii, 'Amplitude'} = u_amp_freq(uaf,1);
        MeanResponseTable{ii, 'Frequency'} = u_amp_freq(uaf,2);
        MeanResponseTable{ii, 'Intensity'} = mean(ResponseTable{uaf_idx, 'IntensityNorm'});
        MeanResponseTable{ii, 'Size'} = mean(ResponseTable{uaf_idx, 'Size'});
        dv = ResponseTable{uaf_idx, 'DigitVec'};
        MeanResponseTable{ii, "NumDigits"} = mean(sum(cat(1, dv{:}),2));
        ii = ii + 1;
    end
    ch_idx = MeanResponseTable.Channel == u_ch(ch);
    MeanResponseTable.IntensityNorm(ch_idx) = MeanResponseTable.Intensity(ch_idx) ./ mean(MeanResponseTable.Intensity(ch_idx));
    MeanResponseTable.SizeNorm(ch_idx) = MeanResponseTable.Size(ch_idx) ./ mean(MeanResponseTable.Size(ch_idx));
end
% 3-way ANOVA
[p,tbl,stats] = anovan(MeanResponseTable{:,"Intensity"}, MeanResponseTable{:,["Channel", "Amplitude", "Frequency"]},...
                       'varnames', {'Channel', 'Amplitude', 'Frequency'}, 'display', 'on');
[p,tbl,stats] = anovan(MeanResponseTable{:,"Size"}, MeanResponseTable{:,["Channel", "Amplitude", "Frequency"]},...
                       'varnames', {'Channel', 'Amplitude', 'Frequency'}, 'display', 'off');

u_conds = unique(MeanResponseTable(:,["Amplitude", "Frequency"]), "rows");
u_c = unique(MeanResponseTable.Channel);
u_a = unique(MeanResponseTable.Amplitude);
u_f = unique(MeanResponseTable.Frequency);


%% Combined figure
SetFont('Arial', 9)

if exist('fig', 'var') 
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax
end

fig = figure('Units', 'inches', 'Position', [5 1 6.48 2]);
ax(1) = axes('Position', [.1 .2 .225 .75]); hold on
    y = zeros(length(u_c), length(u_a));
    for a = 1:length(u_a)
        y(:,a) = MeanResponseTable.IntensityNorm(MeanResponseTable.Amplitude == u_a(a) & MeanResponseTable.Frequency == u_f(2));
    end
    AlphaLine(u_a, y, SubjectColors('BCI02'))
    scatter(u_a, mean(y,1), 50, 'MarkerEdgeColor', SubjectColors('BCI02'),...
        'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 1, 'LineWidth', 1)
    xlabel(sprintf('Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel('Normalized Rating')
    text(40, 2*.95, sprintf('%s Magnitude\n- - Size', GetUnicodeChar('HBar')), 'VerticalAlignment','top')

    y = zeros(length(u_c), length(u_a));
    for a = 1:length(u_a)
        y(:,a) = MeanResponseTable.SizeNorm(MeanResponseTable.Amplitude == u_a(a) & MeanResponseTable.Frequency == u_f(2));
    end
    AlphaLine(u_a, y, SubjectColors('BCI02'), 'LineStyle', '--')
    scatter(u_a, mean(y,1), 50, 'MarkerEdgeColor', SubjectColors('BCI02'),...
        'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 1, 'LineWidth', 1)
    xlabel(sprintf('Amplitude (%sA)', GetUnicodeChar('mu')))
    set(gca, 'YLim', [0 2], 'YTick', [0:.5:2],'XLim', [36 84])

ax(2) = axes('Position', [.3875 .2 .225 .75]); hold on
    y = zeros(length(u_c), length(u_f));
    for f = 1:length(u_a)
        y(:,f) = MeanResponseTable.IntensityNorm(MeanResponseTable.Amplitude == u_a(2) & MeanResponseTable.Frequency == u_f(f));
    end
    AlphaLine(u_f, y, SubjectColors('BCI02'), 'LineStyle', '--')
    scatter(u_f, mean(y,1), 50, 'MarkerEdgeColor', SubjectColors('BCI02'), 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 1, 'LineWidth', 1)
    xlabel('Frequency (Hz)')
    amp_size_ratio = y(:,3) ./ y(:,1);
    fprintf('Amplitude:Size scaling: median (Q1, Q3) = %0.2f, (%0.2f, %0.2f)\n', prctile(amp_size_ratio, [50, 25, 75]))
    
    y = zeros(length(u_c), length(u_f));
    for f = 1:length(u_a)
        y(:,f) = MeanResponseTable.SizeNorm(MeanResponseTable.Amplitude == u_a(2) & MeanResponseTable.Frequency == u_f(f));
    end
    AlphaLine(u_f, y, SubjectColors('BCI02'))
    scatter(u_f, mean(y,1), 50, 'MarkerEdgeColor', SubjectColors('BCI02'), 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 1, 'LineWidth', 1)
    set(gca, 'YLim', [0 2], 'YTick', [0:.5:2], 'XLim', [35 215], 'YTickLabels', {})
    freq_size_ratio = y(:,3) ./ y(:,1);
    fprintf('Frequency:Size scaling: median (Q1, Q3) = %0.2f, (%0.2f, %0.2f)\n', prctile(freq_size_ratio, [50, 25, 75]))

ax(3) = axes('Position', [.675 .2 .225 .75]); hold on
    p1 = polyfit(MeanResponseTable.IntensityNorm, MeanResponseTable.SizeNorm, 1);
    plot([0 3], polyval(p1, [0 3]), "Color", [.6 .6 .6], 'LineStyle','--')
    scatter(MeanResponseTable.IntensityNorm, MeanResponseTable.SizeNorm, 50, SubjectColors('BCI02'), 'filled')
    [r,p] = corr(MeanResponseTable.IntensityNorm, MeanResponseTable.SizeNorm);
    set(gca, 'YTick', [0:.5:2], 'XTick', [0:.5:2], 'XLim', [0 2], 'YLim', [0 2])
    xlabel('Normalized Magnitude')
    ylabel('Normalized Size')

AddFigureLabels(ax, [0.05, 0]); shg

export_path = fullfile(DataPath, 'SuppFig4_AmpFreq');
print(gcf, export_path, '-dsvg', '-r300')


