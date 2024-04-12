%% Imports
addpath("HelperFunctions")
load(fullfile(DataPath, 'ForceEquiv_Pitt.mat')) % Pitt mech amp
load(fullfile(DataPath, "ForceEquiv_UC.mat")) % Mech amp est
load(fullfile(DataPath, "ForceEquiv_Control.mat")) % Mech amp est

% Plotting setup
SetFont('Arial', 9)
scatter_size = 30;
line_width = 1.5;
force_ratio = 0.37;

icms_color = [.6 .6 .6];
mech_color = rgb(0, 137, 123);
multich_color = rgb(244, 67, 54);

force2amp = @(pow_coeffs, lin_coeffs, x) (pow_coeffs(1) * (x + pow_coeffs(3)).^pow_coeffs(2) - lin_coeffs(1)) ./ lin_coeffs(2);
amp2force = @(pow_coeffs, lin_coeffs, x) (((lin_coeffs(2) * x + lin_coeffs(1)) ./ pow_coeffs(1)) .^ (1/pow_coeffs(2))) - pow_coeffs(3);

%% Data extraction
%%% Mag est extraction
mag_est = NaN(length(AS5EData) + length(AS5pData),5);
mag_est_col = NaN(length(AS5EData) + length(AS5pData),3);
x_range = [40,80];
mag_est_x = [40:10:80];
mag_est_r = NaN(size(mag_est_col,1),1);
ii = 1;
% BCI02
for i = 1:length(AS5EData)
    [x,y] = get_condition_magnitudes(AS5EData(i).ResponseSummary, [0,0,1]); % Catch trials only
    mag_est(ii,:) = mean(y,1);
    mag_est_col(ii,:) = SubjectColors('BCI02');
    r = corrcoef(mag_est_x, mag_est(ii,:), 'Rows', 'complete');
    mag_est_r(i) = r(1,2);
    ii = ii + 1;
end
% CRS02 + 07
for i = 1:size(AS5pData,2)
    if i == 9
        ii = ii + 1;
        continue
    end
    x_alt = AS5pData(i).ResponseTable.Amplitude;
    y = AS5pData(i).ResponseTable.MeanNormRating;
    for j = 1:length(mag_est_x)
        xi = find(x_alt == mag_est_x(j));
        if ~isempty(xi)
            mag_est(ii,j) = y(xi);
        end
    end
    mag_est_col(ii,:) = SubjectColors(AS5pData(i).Subject);
    r = corrcoef(mag_est_x, mag_est(ii,:), 'Rows', 'complete');
    mag_est_r(ii) = r(1,2);
    ii = ii + 1;
end

% Fits
[mech_lin_r2, mech_pow_r2, single_icms_lin_r2, single_icms_pow_r2] = deal(NaN(length(NC1Processed),1));
for i = 1:length(NC1Processed)
    mech_lin_r2(i) = NC1Processed(i).MechFit.LinGOF.adjrsquare;
    mech_pow_r2(i) = NC1Processed(i).MechFit.PowGOF.adjrsquare;
    single_icms_lin_r2(i) = NC1Processed(i).ICMSFit.LinGOF.adjrsquare;
    single_icms_pow_r2(i) = NC1Processed(i).ICMSFit.PowGOF.adjrsquare;
end

pow_exps = zeros(size(NC1Processed));
for i = 1:length(NC1Processed)
    pow_exps(i) = NC1Processed(i).MechFit.PowFun.b;
end

%% Plot
if exist('fig', 'var') 
    if isgraphics(fig)
        close(fig)
    end
    clearvars fig ax
end

fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 5, 2.5, 8.75]);

ax(1) = axes('Position', [0.2 0.8 0.675 .175]); hold on
    for i = 1:size(mag_est,1)
        if ~all(isnan(mag_est(i,:)))
            plot(mag_est_x, mag_est(i,:), 'Color', mag_est_col(i,:), 'LineWidth', line_width)
        end
    end
    p = polyfit(mag_est_x, mean(mag_est,1, 'omitnan'),1);
    plot(x_range, polyval(p, x_range), 'Color', 'k', 'LineWidth', 2)
    text(80, 0.55, sprintf('r = %0.2f', median(mag_est_r, 'omitnan')), 'VerticalAlignment','bottom', ...
        'HorizontalAlignment', 'right', 'Color', [.6 .6 .6])

    leg_text = ColorText({'C1', 'P2', 'P3'},...
        [SubjectColors('BCI02'); SubjectColors('CRS02b'); SubjectColors('CRS07')]);
    text(40, 1.5, leg_text, 'VerticalAlignment','top')
    
    xlabel(ax(1), sprintf('Stimulus Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel(ax(1), 'Normalized Intensity')
    set(ax(1), 'YTick', [.5:.5:1.5],...
               'XLim', [38, 82], ...
               'XTick', [40:10:80],...
               'YLim', [.5 1.5])

ax(2) = axes('Position', [0.2 0.55 0.675 .175]); hold on
    ioi = 4;
    % Plot mech
    x = linspace(0,1,size(NC1Processed(ioi).MechSummary,1));
    y = mean(horzcat(NC1Processed(ioi).MechSummary.Responses{:}),1);
    scatter(x, y, 'MarkerFaceColor', mech_color, ...
        'MarkerEdgeColor', mech_color)
    x = [NC1Processed(ioi).MechSummary.Amplitude];
    xq = linspace(min(x), max(x));
    plot(linspace(0,1), NC1Processed(ioi).MechFit.PowFun(xq), 'Color', mech_color, 'LineWidth', 1.5)
    % Plot ICMS
    x = linspace(0,1,size(NC1Processed(ioi).ICMSSummary,1));
    y = mean(horzcat(NC1Processed(ioi).ICMSSummary.Responses{:}),1);
    scatter(x, y, 'MarkerFaceColor', icms_color, ...
        'MarkerEdgeColor', icms_color)
    p = polyfit(x,y',1);
    plot([0,1], polyval(p,[0,1]), 'Color', icms_color, 'LineWidth', 1.5)

    % Format
    yl = [.7, 1.2];
    ylabel(ax(2), 'Normalized Intensity')
    set(ax(2), 'XColor', 'none',...
             'XLim', [-0.05 1.05],...-
             'YLim', yl,...
             'YTick', [.6:.2:1.2],...
             'Clipping', 'off')
    % Legend
    text(1,yl(1)+range(yl)*0.05, ColorText({'ICMS', 'Mechanical'}, [icms_color; mech_color]), 'HorizontalAlignment', 'right', ...
        'VerticalAlignment','bottom')
    % Mech x-axis
    plot([0, 1], [yl(1), yl(1)], 'Color', mech_color, 'LineWidth', 1)
    x = NC1Processed(ioi).MechSummary.Amplitude(1:2:end) .* force_ratio;
    xs = linspace(0,1,length(x));
    for xsi = 1:length(xs)
        plot([xs(xsi), xs(xsi)], [yl(1), yl(1)-range(yl)*0.025], 'Color', mech_color, 'LineWidth', 1)
        if xsi < length(xs)
            text(xs(xsi), yl(1)-range(yl)*0.025, sprintf('%0.2f', x(xsi)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', 'Color', mech_color)
        else
            text(xs(xsi), yl(1)-range(yl)*0.025, sprintf('%0.2f N', x(xsi)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', 'Color', mech_color)
        end
    end

    % ICMS x-axis
    plot([0, 1], [yl(1), yl(1)]-range(yl)*0.2, 'Color', icms_color, 'LineWidth', 1)
    x = NC1Processed(ioi).ICMSSummary.Amplitude(1:2:end);
    xs = linspace(0,1,length(x));
    for xsi = 1:length(xs)
        plot([xs(xsi), xs(xsi)], [yl(1), yl(1)-range(yl)*0.025]-range(yl)*0.2, 'Color', icms_color, 'LineWidth', 1)
        if xsi < length(xs)
            text(xs(xsi), yl(1)-range(yl)*0.025-range(yl)*0.2, num2str(x(xsi)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', 'Color', icms_color)
        else
            text(xs(xsi), yl(1)-range(yl)*0.025-range(yl)*0.2, sprintf('%d %sA', x(xsi), GetUnicodeChar('mu')), ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', 'Color', icms_color)
        end
    end

ax(3) = axes('Position', [0.2 0.29 0.675 .175]); hold on
    c = lines(length(NC1Processed));
    for i = 1:length(NC1Processed)
        x = linspace(0, max(NC1Processed(i).MechRange(2)));
        y = force2amp(coeffvalues(NC1Processed(i).MechFit.PowFun), coeffvalues(NC1Processed(i).ICMSFit.LinFun), x);
        if i == ioi
            continue
        else
            plot(x.* force_ratio, y, 'Color', [.6 .6 .6], 'LineWidth', line_width);
        end
    end
    x = linspace(0, max(NC1Processed(ioi).MechRange(2)));
    y = force2amp(coeffvalues(NC1Processed(ioi).MechFit.PowFun), coeffvalues(NC1Processed(ioi).ICMSFit.LinFun), x);
    plot(x .* force_ratio, y, 'Color', 'k', 'LineWidth', line_width+.5);

    ylabel(ax(3), sprintf('ICMS Amplitude (%sA)', GetUnicodeChar('mu')))
    xlabel(ax(3), 'Stimulus Force (N)')
    set(ax(3), 'YTick', [0:25:100],...
               'XLim', [0 .5], ...
               'XTick', [0:.1:.5],...
               'YLim', [0, 100],...
               'XTickLabelRotation', 0)

ax(4) = axes('Position', [0.2 0.05 0.675 .175]); hold on
    plot([.5 1], [.5 1], 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(mech_lin_r2, mech_pow_r2, 30, 'MarkerFaceColor', mech_color, 'MarkerEdgeColor', mech_color)
    scatter(single_icms_lin_r2, single_icms_pow_r2, 30, 'MarkerFaceColor', [.6 .6 .6], 'MarkerEdgeColor', [.6 .6 .6])

    set(ax(4), 'XLim', [.5 1], 'YLim', [.5 1], 'XTick', [.5:.25:1], 'YTick', [.5:.25:1])
    xl = xlabel(ax(4), 'Linear Adj.R^2');
    xl.Position(2) = xl.Position(2) + .01;
    ylabel(ax(4), 'Power Adj.R^2')

    text(.975, .525, ColorText({'Single ICMS', 'Mechanical'}, [[.6 .6 .6];mech_color]), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment','right')
