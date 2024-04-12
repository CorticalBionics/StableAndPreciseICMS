%% Imports
addpath("HelperFunctions")
load(fullfile(DataPath, 'DigitLocalizationData'))
DigitLocalizationData = rmfield(DigitLocalizationData,...
    {'Frequency', 'Images', 'Responses', 'TestLogComments', 'Paradigm'});
load(fullfile(DataPath, "RoboticDigitLocalization.mat"))

%% Make response tables
response_table_config = [["SesID", "double"]; ...
                         ["SetID", "double"]; ...
                         ["Block", "double"]; ...
                         ["Trial", "double"]; ...
                         ["Channel", "cell"];...
                         ["TargetDigits", "string"];...
                         ["IsMultiDigit", "logical"];...
                         ["IsMultiChannel", "logical"];...                         
                         ["ResponseDigits", "string"];...
                         ["Correct", "double"]];

for i = 1:size(DigitLocalizationData, 2)
    num_sdo = size(DigitLocalizationData(i).SDO,2);
    response_table = table('Size',[num_sdo, size(response_table_config,1)],...
                           'VariableNames', response_table_config(:,1),...
                           'VariableTypes', response_table_config(:,2));
    for j = 1:num_sdo
        response_table{j, "SesID"} = DigitLocalizationData(i).SDO(j).sessionInfo.session_num;
        response_table{j, "SetID"} = DigitLocalizationData(i).SDO(j).set;
        response_table{j, "Block"} = DigitLocalizationData(i).SDO(j).block;
        response_table{j, "Trial"} = DigitLocalizationData(i).SDO(j).trialID;
        response_table{j, "Channel"} = DigitLocalizationData(i).SDO(j).channel(1);
        f_split = strsplit(DigitLocalizationData(i).SDO(j).StimulusDelivered{1}, '_');
        if contains(f_split{end}, '.')
            response_table{j, "IsMultiDigit"} = true;
            num_digits = sum(f_split{end} == '.') + 1;
        end
        if strcmp(f_split{end}, "Ch0")
            response_table{j, "TargetDigits"} = "None";
        else
            response_table{j, "TargetDigits"} = f_split(end);
        end
        if length(response_table{j, "Channel"}{1}) / num_digits > 1
            response_table{j, "IsMultiChannel"} = true;
        end
        if isnan(DigitLocalizationData(i).SDO(j).reportedData.value)
            response_table{j, "ResponseDigits"} = "Discarded";
            response_table{j, "Correct"} = NaN;
            continue
        else
            response_table{j, "ResponseDigits"} = {DigitLocalizationData(i).SDO(j).reportedData.value};
        end
        response_split = strsplit(response_table{j, "ResponseDigits"}, '.');
        for k = 1:length(response_split)
            if contains(response_table{j, "TargetDigits"}, response_split{k})
                response_table{j, "Correct"} = response_table{j, "Correct"} + 1;
            end
        end
        
    end
    DigitLocalizationData(i).ResponseTable = response_table(~isnan(response_table.Correct),:);
end

%% Analyze combined response table - OLS
ResponseTable = cat(1, DigitLocalizationData.ResponseTable);
multi_ch_idx = logical([ResponseTable.IsMultiChannel]);
multi_digit_idx = logical([ResponseTable.IsMultiDigit]);
null_idx = strcmp(ResponseTable.TargetDigits, "None");
no_percept_idx = strcmp(ResponseTable.ResponseDigits, "None");
ols_correct_idx = [ResponseTable.Correct] > 0;

num_bootstraps = 50;
[sd_sc, sd_mc, md_sc, md_mc] = deal(zeros(num_bootstraps,1));
% OLS - bootstrap
sd_sc_idx = find(~multi_digit_idx & ~multi_ch_idx & ~null_idx); % Single digit, single channel
sd_mc_idx = find(~multi_digit_idx & multi_ch_idx & ~null_idx); % Single digit, multi channel
md_sc_idx = find(multi_digit_idx & ~multi_ch_idx & ~null_idx); % Multi digit, single channel
md_mc_idx = find(multi_digit_idx & multi_ch_idx & ~null_idx); % Multi digit, multi channel

for b = 1:num_bootstraps
    temp_idx = datasample(sd_sc_idx, 10, 'Replace', false);
    sd_sc(b) = mean([ResponseTable.Correct(temp_idx)]);
    temp_idx = datasample(sd_mc_idx, 10, 'Replace', false);
    sd_mc(b) = mean([ResponseTable.Correct(temp_idx)]);
    temp_idx = datasample(md_sc_idx, 10, 'Replace', false);
    md_sc(b) = mean([ResponseTable.Correct(temp_idx)]) / 2;
    temp_idx = datasample(md_mc_idx, 10, 'Replace', false);
    md_mc(b) = mean([ResponseTable.Correct(temp_idx)]) / 2;
end

% Stats tests
[h,p] = BinomTest2(ResponseTable.Correct(sd_sc_idx), ResponseTable.Correct(sd_mc_idx));
[h,p] = BinomTest2(ResponseTable.Correct(md_sc_idx), ResponseTable.Correct(md_mc_idx));

OLS_sd_sc = ResponseTable.Correct(sd_sc_idx);
OLS_sd_mc = ResponseTable.Correct(sd_mc_idx);
OLS_md_sc = ResponseTable.Correct(md_sc_idx);
OLS_md_mc = ResponseTable.Correct(md_mc_idx);

unique_target_digits = unique(ResponseTable.TargetDigits(~null_idx));
utm_idx = contains(unique_target_digits, '.');
usd = [unique_target_digits(~utm_idx); "None"];
umd = [unique_target_digits(utm_idx); "None"];

% Single digit confusion matrix
[sd_response_matrix_single, sd_response_matrix_multi] = deal(zeros(length(usd)-1, length(usd)));
for td_idx = 1:length(usd)-1
    idx = find(strcmp(ResponseTable.TargetDigits, usd{td_idx}));
    for i = 1:length(idx) % Indices of matching target digits
        for rd_idx = 1:length(usd)-1
            if strcmp(usd{rd_idx}, ResponseTable{idx(i), "ResponseDigits"})
                if ResponseTable{idx(i), "IsMultiChannel"}
                    sd_response_matrix_multi(td_idx,rd_idx) = sd_response_matrix_multi(td_idx,rd_idx) + 1;
                else
                    sd_response_matrix_single(td_idx,rd_idx) = sd_response_matrix_single(td_idx,rd_idx) + 1;
                end
            end
        end
    end
end
sd_single_rate_matrix = sd_response_matrix_single ./ sum(sd_response_matrix_single,2);
sd_multi_rate_matrix = sd_response_matrix_multi ./ sum(sd_response_matrix_multi,2);

% Multi digit copnfusion matrix
[md_response_matrix_single, md_response_matrix_multi] = deal(zeros(length(umd)-1, length(umd)-1));
for td_idx = 1:length(umd)-1
    idx = find(strcmp(ResponseTable.TargetDigits, umd{td_idx}));
    for i = 1:length(idx) % Indices of matching target digits
        for j = 1:length(response_split)
            for rd_idx = 1:length(umd)-1
                if strcmp(umd{rd_idx}, ResponseTable{idx(i), "ResponseDigits"})
                    if ResponseTable{idx(i), "IsMultiChannel"}
                        md_response_matrix_multi(td_idx,rd_idx) = md_response_matrix_multi(td_idx,rd_idx) + 1;
                    else
                        md_response_matrix_single(td_idx,rd_idx) = md_response_matrix_single(td_idx,rd_idx) + 1;
                    end
                end
            end
        end
    end
end
md_single_rate_matrix = md_response_matrix_single ./ sum(md_response_matrix_single,2);
md_multi_rate_matrix = md_response_matrix_multi ./ sum(md_response_matrix_multi,2);

%% Combined heatmap
[comb_sc_confmat, comb_mc_confmat] = deal(zeros(length(unique_target_digits)));
for i = 1:size(ResponseTable,1)
    if strcmp(ResponseTable{i, "ResponseDigits"}, "None")
        continue
    end
    row = strcmp(ResponseTable{i, "TargetDigits"}, unique_target_digits);
    col = strcmp(ResponseTable{i, "ResponseDigits"}, unique_target_digits);
    if ResponseTable{i, "IsMultiChannel"}
        comb_mc_confmat(row, col) = comb_mc_confmat(row, col) + 1;
    else
        comb_sc_confmat(row, col) = comb_sc_confmat(row, col) + 1;
    end
end

comb_sc_confmat = comb_sc_confmat ./ sum(comb_sc_confmat, 1);
comb_mc_confmat = comb_mc_confmat ./ sum(comb_mc_confmat, 1);

%% Analyze robotic localization
rob_nostim_idx = strcmp(RoboticDigitLocalizationTable.TargetDigits, "None");
rob_noperceive_idx = strcmp(RoboticDigitLocalizationTable.ResponseDigits, "None");
rob_multichan_idx = RoboticDigitLocalizationTable.IsMultiChannel;
rob_multidigit_idx = RoboticDigitLocalizationTable.IsMultiDigit;
false_alarm_rate = mean(rob_nostim_idx & ~rob_noperceive_idx);
no_perceive_single_chan = mean(~rob_nostim_idx & rob_noperceive_idx & ~rob_multichan_idx);

% Single channel single digit
sd_sc_idx = find(~rob_multichan_idx & ~rob_multidigit_idx & ~rob_nostim_idx);
% FInd unique channels
sd_sd_channels = RoboticDigitLocalizationTable.Channel(sd_sc_idx);
sd_sd_channels = unique(cat(1,sd_sd_channels{:}));
% Find matching channel
ic = zeros(size(sd_sc_idx));
for i = 1:length(sd_sc_idx)
    ic(i) = find(ismember(sd_sd_channels, RoboticDigitLocalizationTable.Channel{sd_sc_idx(i)}));
end
% Indices of each channel
sch_idx = cell(size(sd_sd_channels));
rob_sd_sc_perf = zeros(size(sd_sd_channels));
for i = 1:length(sd_sd_channels)
    sch_idx{i} = sd_sc_idx(ic == i);
    rob_sd_sc_perf(i) = mean(RoboticDigitLocalizationTable.Correct(sd_sc_idx(ic == i)));
end

% Multi channel single digit
sd_mc_idx = find(rob_multichan_idx & ~rob_multidigit_idx & ~rob_nostim_idx);
% FInd unique channels
multi_digit_channels = RoboticDigitLocalizationTable.Channel(sd_mc_idx);
multi_digit_channels = unique(cat(1,multi_digit_channels{:}), 'Rows');
% Find matching channel
ic = zeros(size(sd_mc_idx));
for i = 1:length(sd_mc_idx)
    ic(i) = find(all(ismember(multi_digit_channels, RoboticDigitLocalizationTable.Channel{sd_mc_idx(i)}),2));
end

% Compare
[rob_sd_mc_perf, comp_sd_sc_perf] = deal(zeros(size(multi_digit_channels,1),1));
for i = 1:size(multi_digit_channels,1)
    mch_idx = sd_mc_idx(ic == i);
    rob_sd_mc_perf(i) = mean(RoboticDigitLocalizationTable.Correct(mch_idx));
    i1 = find(sd_sd_channels == multi_digit_channels(i,1));
    i2 = find(sd_sd_channels == multi_digit_channels(i,2));
    comp_sd_sc_perf(i) = max(rob_sd_sc_perf([i1,i2]));
end

% Single channel multi digit
md_sc_idx = find(~rob_multichan_idx & rob_multidigit_idx & ~rob_nostim_idx);
% Find unique channels
md_sc_channels = RoboticDigitLocalizationTable.Channel(md_sc_idx);
temp = zeros(size(md_sc_channels,1),2);
for i = 1:size(md_sc_channels,1)
    temp(i,1:length(md_sc_channels{i})) = md_sc_channels{i};
end
md_sc_channels = unique(temp, 'rows');

% Find matching channel
ic = zeros(size(md_sc_idx));
for i = 1:length(md_sc_idx)
    ic(i) = find(all(md_sc_channels == RoboticDigitLocalizationTable.Channel{md_sc_idx(i)}, 2));
end
% Indices of each channel
sch_md_idx = cell(size(md_sc_channels,1));
rob_md_sc_perf = zeros(size(md_sc_channels,1),1);
for i = 1:length(md_sc_channels)
    sch_md_idx{i} = md_sc_idx(ic == i);
    rob_md_sc_perf(i) = mean(RoboticDigitLocalizationTable.Correct(sch_md_idx{i}));
end

% Multi channel multi digit
md_mc_idx = find(rob_multichan_idx & rob_multidigit_idx & ~rob_nostim_idx);
% Find unique channels
md_mc_channels = RoboticDigitLocalizationTable.Channel(md_mc_idx);
temp = zeros(size(md_mc_channels,1),2);
for i = 1:size(md_mc_channels,1)
    temp(i,1:length(md_mc_channels{i})) = md_mc_channels{i};
end
md_mc_channels = unique(temp, 'rows');

% Find matching channel
ic = zeros(size(md_mc_idx));
for i = 1:length(md_mc_idx)
    ic(i) = find(all(md_mc_channels == RoboticDigitLocalizationTable.Channel{md_mc_idx(i)}, 2));
end
% Indices of each channel
mch_md_idx = cell(size(md_mc_channels));
[rob_md_mc_perf, comp_md_mc_perf, ] = deal(zeros(size(md_mc_channels,1),1));
for i = 1:length(md_mc_channels)
    mch_md_idx{i} = md_mc_idx(ic == i);
    rob_md_mc_perf(i) = mean(RoboticDigitLocalizationTable.Correct(mch_md_idx{i}));
    % Find matching single channels
    comp_idx = find(all(ismember(md_sc_channels, md_mc_channels(i,:)),2));
    if rob_md_sc_perf(comp_idx(1)) > rob_md_sc_perf(comp_idx(2))
        comp_md_mc_perf(i) = rob_md_sc_perf(comp_idx(1));
    else
        comp_md_mc_perf(i) = rob_md_sc_perf(comp_idx(2));
    end
end

RD_sd_sc = RoboticDigitLocalizationTable.Correct(sd_sc_idx);
RD_sd_mc = RoboticDigitLocalizationTable.Correct(sd_mc_idx);
RD_md_sc = RoboticDigitLocalizationTable.Correct(md_sc_idx);
RD_md_mc = RoboticDigitLocalizationTable.Correct(md_mc_idx);

%% Statistical tests
% Combined
[p,h,s] = ranksum([OLS_sd_sc; RD_sd_sc], [OLS_sd_mc; RD_sd_mc]);
[p,h,s] = ranksum([OLS_md_sc; RD_md_sc], [OLS_md_mc; RD_md_mc]);
% Robotic
[p,h,s] = ranksum(RD_sd_sc, RD_sd_mc);
[p,h,s] = ranksum(OLS_sd_sc, OLS_sd_mc);
[p,h,s] = ranksum(RD_md_sc, RD_md_mc);
[p,h,s] = ranksum(OLS_md_sc, OLS_md_mc);


%% Supplement Plot
sc_color = rgb(158, 158, 158);
mc_color = rgb(244, 67, 54);
clf; set(gcf, 'Units', 'Inches', 'Position', [30, 1, 6.48, 4]);
axes('Position', [.075 .25 .15 .5])
    Swarm(1, [rob_sd_sc_perf; rob_md_sc_perf; rob_sd_mc_perf; rob_md_mc_perf], 'GroupName', 'Single Electrode, Robot','Color', sc_color,...
        'DS', ds, 'SPL', spl, 'CenterMethod', 'Median', 'ShowStats', true, 'HS', '/')
    Swarm(2, [sd_sc; md_sc; sd_mc; md_mc], 'GroupName', 'Single Channel, OLS','Color', sc_color,...
        'DS', ds, 'SPL', spl, 'CenterMethod', 'Median', 'ShowStats', true)

    text(1, 0.025, {'Robot', 'Hand'}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment','center')
    text(2, 0.025, {'Computer', 'Controlled'}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment','center')

    set(gca, 'XTick', [1, 2], 'XTickLabel', {''}, 'XLim', [.5 2.5], 'YTick', [0:.2:1], 'YLim', [0 1])
    ylabel('Performance')
    

cmap = [linspace(.99,.55,255); linspace(.99,.14,255); linspace(.99,.66,255)]';
axes('Position', [.35 .1 .3 .8])
    imagesc(comb_sc_confmat)
    colormap(cmap)
    set(gca, 'XLim', [.5 10.5], 'Ylim', [.5 10.5], 'YDir', 'Reverse', 'YTick', [1:10],...
        'YTickLabel', unique_target_digits, 'XTick', [1:10], 'XTickLabel', unique_target_digits, 'CLim', [0 1],...
        'PlotBoxAspectRatio', [1,1,1])
    title('Single Electrode')
    ylabel('Reported Digit(s)')
    xlabel('Target Digit(s)')

axes('Position', [.675 .1 .3 .8])
    imagesc(comb_mc_confmat)
    colormap(cmap)
    set(gca, 'XLim', [.5 10.5], 'Ylim', [.5 10.5], 'YDir', 'Reverse', 'YTick', [1:10],...
        'XTick', [1:10], 'YTickLabel', {}, 'XTickLabel', unique_target_digits, 'CLim', [0 1],...
        'PlotBoxAspectRatio', [1,1,1])
    title('Multi Electrode')
    xlabel('Target Digit(s)')


% AddFigureLabels(gcf, [0.05 0]); shg
annotation("textbox", [.025 .8 .05 .05], 'String', 'A', ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [.325 .8 .05 .05], 'String', 'B', ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [.65 .8 .05 .05], 'String', 'C', ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')


%% Helper function
function perf_vec = block_target_performance(temp)
    [c,ia,ic] = unique(temp{:, ["Block", "TargetDigits"]}, 'Rows');
    perf_vec = zeros(size(c,1),1);
    for i = 1:length(c)
        perf_vec(i) = mean(temp{ic == i, "Correct"});
    end
end
