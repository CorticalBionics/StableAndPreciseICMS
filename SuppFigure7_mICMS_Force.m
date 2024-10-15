%% Imports
addpath("HelperFunctions")
load(fullfile(DataPath, "ForceEquivMulti.mat")) % Mech amp est 
load(fullfile(DataPath, "AmpDiscrimMultiChanRaw")) % multi chan discrim
load(fullfile(DataPath, "ForceEquiv_UC.mat")) % Mech amp est

SetFont('Arial', 9)
pow_fun = fittype('a*(x+o)^b');
lin_fun = fittype('m*x+c');
force2amp = @(pow_coeffs, lin_coeffs, x) (pow_coeffs(1) * (x + pow_coeffs(3)).^pow_coeffs(2) - lin_coeffs(1)) ./ lin_coeffs(2);
amp2force = @(pow_coeffs, lin_coeffs, x) (((lin_coeffs(2) * x + lin_coeffs(1)) ./ pow_coeffs(1)) .^ (1/pow_coeffs(2))) - pow_coeffs(3);

%% Data extraction - Variability - multi
[mech_multi_cv, icms_single_ch_cv, icms_multi_ch_cv] = deal(cell(length(NC3Processed), 1));
[cv_flat, stim_type] = deal([]); % kinda messy but whatever
ii = 1;
for i = 1:length(NC3Processed)
    % Mech
    temp = NaN(size(NC3Processed(i).MechSummary,1),1);
    for j = 1:size(NC3Processed(i).MechSummary,1)
        temp(j) = std(NC3Processed(i).MechSummary.Responses{j}) / mean(NC3Processed(i).MechSummary.Responses{j});
    end
    cv_flat(:,ii) = temp; 
    stim_type(:,ii) = ones(size(NC3Processed(i).MechSummary,1),1);
    ii = ii +1;
    mech_multi_cv{i} = temp;
    
    % ICMS Single
    temp = NaN(size(NC3Processed(i).ICMSSummary.SingleChResponses,1), size(NC3Processed(i).ICMSSummary.SingleChResponses,2));
    for j = 1:size(NC3Processed(i).ICMSSummary.SingleChResponses,2)
        for k = 1:size(NC3Processed(i).ICMSSummary.SingleChResponses,1)
            temp(k,j) = std(NC3Processed(i).ICMSSummary.SingleChResponses{k,j}) / mean(NC3Processed(i).ICMSSummary.SingleChResponses{k,j});
        end
        cv_flat(:,ii) = temp(:,j);
        stim_type(:,ii) = ones(size(NC3Processed(i).ICMSSummary.SingleChResponses,1),1) * 2;
        ii = ii +1;
    end
    icms_single_ch_cv{i} = temp;

    % ICMS Multi
    temp = NaN(size(NC3Processed(i).ICMSSummary.MultiChResponses));
    for j = 1:size(NC3Processed(i).ICMSSummary.MultiChResponses,1)
        temp(j) = std(NC3Processed(i).ICMSSummary.MultiChResponses{j}) / mean(NC3Processed(i).ICMSSummary.MultiChResponses{j});
    end
    cv_flat(:,ii) = temp;
    stim_type(:,ii) = ones(size(NC3Processed(i).ICMSSummary.MultiChResponses,1),1) * 3;
    ii = ii +1;
    icms_multi_ch_cv{i} = temp;
end

stim_int = repmat([1:8]', 1, size(cv_flat,2));

[p,tbl,stats] = anovan(cv_flat(:), {stim_type(:), stim_int(:)}, 'VarNames', {'StimType', 'Intensity'});


%% Data extraction - JNDs
mc_jnds = NaN(length(AS4mcEData),1);
for i = 1:length(AS4mcEData)
    mc_jnds(i) = AS4mcEData(i).JNDTable{"nc","pCH"} / 100 / 5;
end
mc_jnds([AS4mcEData.IsMulti]) = mc_jnds([AS4mcEData.IsMulti]) * 4;

[P,H,STATS] = ranksum(mc_jnds([AS4mcEData.IsMulti]), mc_jnds(~[AS4mcEData.IsMulti]));
fprintf('Single vs Multi-Channel JND (Ranksum test):\n Z = %0.2f\n %s\n', STATS.zval, pStr(P))

%% Data extraction - FitType
[mech_lin_r2, mech_pow_r2, single_icms_lin_r2, single_icms_pow_r2] = deal(NaN(length(NC1Processed),1));
for i = 1:length(NC1Processed)
    mech_lin_r2(i) = NC1Processed(i).MechFit.LinGOF.adjrsquare;
    mech_pow_r2(i) = NC1Processed(i).MechFit.PowGOF.adjrsquare;
    single_icms_lin_r2(i) = NC1Processed(i).ICMSFit.LinGOF.adjrsquare;
    single_icms_pow_r2(i) = NC1Processed(i).ICMSFit.PowGOF.adjrsquare;
end

[multi_icms_lin_r2, multi_icms_pow_r2] = deal(NaN(length(NC3Processed),1));
for i = 1:length(NC3Processed)
    multi_icms_lin_r2(i) = NC3Processed(i).ICMSFitMulti.LinGOF.adjrsquare;
    multi_icms_pow_r2(i) = NC3Processed(i).ICMSFitMulti.PowGOF.adjrsquare;
end

%% Plotting
multich_color = rgb(244, 67, 54);
control_color = [.6 .6 .6];
mech_color = rgb(0, 137, 123);

if exist('fig', 'var') 
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax
end

fig = figure;
set(fig, 'Units', 'inches', 'Position', [30, 5, 6.25, 2]);
t = tiledlayout(1,3, "TileSpacing", "compact", "Padding", "compact");

ax(1) = nexttile; hold on
    plot([.5 1], [.5 1], 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(mech_lin_r2, mech_pow_r2, 30, 'MarkerFaceColor', mech_color, 'MarkerEdgeColor', mech_color)
    scatter(single_icms_lin_r2, single_icms_pow_r2, 30, 'MarkerFaceColor', control_color, 'MarkerEdgeColor', control_color)
    scatter(multi_icms_lin_r2, multi_icms_pow_r2, 30, 'MarkerFaceColor', multich_color, 'MarkerEdgeColor', multich_color)

    set(ax(1), 'XLim', [.5 1], 'YLim', [.5 1], 'XTick', [.5:.25:1], 'YTick', [.5:.25:1])
    xlabel(ax(1), 'Linear Adj. R^2')
    ylabel(ax(1), 'Power Adj. R^2')

    text(.975, .525, ColorText({'Mechanical', 'Bio Multi ICMS', 'Bio Single ICMS'}, [mech_color;multich_color;control_color]), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment','right')

ax(2) = nexttile; hold on
    Swarm(1, cellfun(@(c) mean(c), mech_multi_cv), mech_color, "DS","box", 'DistributionMethod','none')
    temp = cellfun(@(c) mean(c,1), icms_single_ch_cv, 'UniformOutput', false);
    Swarm(2, cat(2, temp{:}), control_color, "DS","box")
    Swarm(3, cellfun(@(c) mean(c), icms_multi_ch_cv), multich_color, "DS","box", 'DistributionMethod','none')
    
    ylabel(ax(2), 'Coefficient of Variation')
    set(ax(2), 'XTick', [1:3], 'XTickLabel', {'Mech', 'ICMS_{BS}', 'ICMS_{BM}'}, 'XLim', [.5 3.5], 'YLim', [0 1.2], 'YTick', [0:.4:1.2])
    text(.7, 1.2*0.95, pStr(tbl{2,7}), 'VerticalAlignment', 'top', 'Color', [.6 .6 .6])
    
ax(3) = nexttile; hold on
    Swarm(1, mc_jnds(~[AS4mcEData.IsMulti]), control_color,"DS","box", 'SYL', [0 60])
    Swarm(2, mc_jnds([AS4mcEData.IsMulti]), multich_color, "DS","box", 'SYL', [0 60])

    ylabel(ax(3), 'JND (nC/p)')
    set(ax(3), 'XTick', [1,2], ...
               'XTickLabel', {'Single', 'Multi'}, ...
               'XLim', [.5 2.5], ...
               'YTick', [0:4:12], ...
               'YLim', [0 12])
    text(.6, 12, {sprintf('Z = %0.2f',STATS.zval); pStr(P)}, 'Color', [.6 .6 .6], 'VerticalAlignment','top', 'HorizontalAlignment','left')



%%% Labels
AddFigureLabels(ax, [0.075 -0.05])
export_path = fullfile(DataPath, 'SuppFig7_MultiChForce');
print(gcf, export_path, '-dpng', '-r300')
print(gcf, export_path, '-depsc', '-r300')