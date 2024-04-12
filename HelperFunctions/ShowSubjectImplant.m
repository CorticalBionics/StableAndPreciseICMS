function [fig, ax, properties] = ShowSubjectImplant(SubjectID, varargin)
    % Get the channel map - will also error if invalid SubjectID
    ChannelMap = LoadSubjectChannelMap(SubjectID);
    % Set defaults
    ArrayColors = SubjectColors(SubjectID);
    SubjectLabelColor = SubjectColors(SubjectID);
    RegionLabelColor = 'k';
    ScaleColor = 'k';
    ShowScale = true;
    SubjectString = SubjectID;
    SubjectBackgroundColor = [1 1 1 .3];
    SplineStyle = '--';
    LineWidth = 2;
    FontSize = 18;
    Parent = [];
    ShowCorner = true;
    CornerLocation = 'Original'; % Original or PostMed
    ParseVarargin()
    
    % Check ArrayColors
    if all(size(ArrayColors) == [1,3])
        ArrayColors = repmat(ArrayColors, [4,1]);
    elseif ~all(size(ArrayColors) == [4,3])
        error('Size of ArrayColors must be [1,3] or [4,3]')
    end
    
    mpath = mfilename('fullpath');
    [parent_folder, ~] = fileparts(mpath);
    img_path = fullfile(DataPath, 'SubjectMRIPhotos');
    img_fullpath = fullfile(img_path, sprintf('%s_ImplantSiteMRI.png', SubjectID));
    
    % Assert plotting in the order (medial motor, medial sensory, lateral
    % motor, lateral sensory) - this must also be the order the colors are
    % given (if multiple)
    array_order = {'MedialMotor', 'MedialSensory', 'LateralMotor', 'LateralSensory'};
    array_idx = zeros(1,4);
    array_sizes = cellfun(@size, ChannelMap.ArrayLocations, 'UniformOutput', false);
    for i = 1:4
        array_idx(i) = find(strcmp(ChannelMap.ArrayNames, array_order{i}));
    end
    array_rotations = ChannelMap.ArrayRotations(array_idx);
    array_sizes = array_sizes(array_idx);
    
    % Declare hardcoded parameters
    switch SubjectID
        case 'BCI02'
            img = imread(img_fullpath);
            array_positions = [250, 575; ... % Medial motor
                               715, 845; ... % Medial sensory
                               315, 800; ... % Lateral motor
                               730, 1050];   % Lateral sensory
            motor_label = [170, 750];
            sensory_label = [600, 765];
            spline_x = [427, 457, 500, 484, 445, 444, 457, 459];
            spline_y = [1323, 1198, 1053, 829, 709, 583, 478, 336];
            array_scalar = 11.1; % Pixels per ch
            im_center = [475, 825];
            im_offset = 850 / 2;
            ml_ap_axis_pos = [401, 1000];

        case 'BCI03'
            img = imread(img_fullpath);
            array_positions = [215, 840; ... % Medial motor
                               777, 1250; ... % Medial sensory
                               480, 946; ... % Lateral motor
                               850, 1510]; % Lateral sensory
            motor_label = [274, 1046];
            sensory_label = [706, 1375];
            spline_x = [572, 554, 582, 643, 631, 680, 707, 628, 604, 604, 651];
            spline_y = [399, 505, 594, 675, 830, 909, 1067, 1238, 1425, 1578, 1677];
            array_scalar = 9.6; % Pixels per ch
            im_center = [550, 1160];
            im_offset = 900 / 2;
            ml_ap_axis_pos = [444, 1406];
    
        case 'CRS02'
            img = imread(img_fullpath);
            array_positions = [350, 550; ... % Medial motor
                               493, 950; ... % Medial sensory
                               325, 790; ... % Lateral motor
                               475, 1172]; % Lateral sensory
            motor_label = [170, 930];
            sensory_label = [592, 827];
            spline_x = [400, 454, 477, 511, 515, 486, 448, 434, 436, 400, 381, 400, 434];
            spline_y = [426, 565, 652, 717, 780, 834, 873, 958, 1045, 1138, 1196, 1263, 1312];
            array_scalar = 8.7; % Pixels per ch
            im_center = [425, 850];
            im_offset = 850 / 2;
            ml_ap_axis_pos = [305, 1100];
    
        case 'CRS07'
            img = imread(img_fullpath);
            array_positions = [322, 585; ... % Medial motor
                               600, 885; ... % Medial sensory
                               320, 760; ... % Lateral motor
                               605, 1050]; % Lateral sensory
            motor_label = [220, 525];
            sensory_label = [653, 661];
            spline_x = [403, 509, 481, 454, 477, 529, 582, 554];
            spline_y = [1268, 1171, 1070, 983, 875, 683, 545, 430];
            array_scalar = 11.3; % Pixels per ch
            im_center = [515, 795];
            im_offset = 750 / 2;
            ml_ap_axis_pos = [320, 920];

        case 'CRS08'
            img = imread(img_fullpath);
            array_positions = [490, 335; ... % Medial motor
                               835, 818; ... % Medial sensory
                               572, 545; ... % Lateral motor
                               795, 961]; % Lateral sensory
            motor_label = [390, 546];
            sensory_label = [774, 700];
            spline_x = [525, 570, 655, 683, 666, 689, 677, 647, 574, 574];
            spline_y = [220, 335, 458, 526, 618, 818, 983, 1062, 1150, 1221];
            array_scalar = 7.3; % Pixels per ch
            im_center = [667, 650];
            im_offset = 900 / 2;
            ml_ap_axis_pos = [561, 781];
    end

    % Store values in the properties struct in-case someone wants it returned
    properties = struct();
    properties.ArrayLabels =  array_order;
    properties.ArrayPositions = array_positions;
    properties.ArrayRotations = array_rotations;

    % Fit the spline
    yq = linspace(min(spline_y), max(spline_y));
    xq = spline(spline_y,spline_x,yq);
    
    % Plotting
    if isempty(Parent)
        fig = figure('Units', 'inches', 'Position', [1 1 5 5]);
        ax = axes('Position', [0 0 1 1], 'XColor', 'none', 'YColor', 'none'); hold on
        Parent = ax;
    end
    % Image
    imshow(img)
    % Find limits
    set(gca, 'XLim', im_center(1) + [-im_offset, im_offset], 'YLim', im_center(2) + [-im_offset, im_offset])
    % Arrays
    for a = 1:4
        ys = array_sizes{a}(1)/2;
        xs = array_sizes{a}(2)/2;
        y_scaled = [-ys ys ys -ys] .* array_scalar;
        x_scaled = [-xs -xs xs xs] .* array_scalar;
        [x_rotated,y_rotated] = rotate_array(x_scaled,y_scaled,array_rotations(a));
        x_final = x_rotated + array_positions(a,1);
        y_final = y_rotated + array_positions(a,2);
        patch(x_final, y_final, ArrayColors(a,:), 'LineWidth', LineWidth, 'EdgeColor',...
            ArrayColors(a,:), 'FaceAlpha', .2, 'Parent', Parent)
        properties.ArrayCoordinates{a,1} = x_final;
        properties.ArrayCoordinates{a,2} = y_final;

        if ShowCorner
            if strcmpi(CornerLocation, 'Original')
                y = [-ys+2 -ys -ys] .* array_scalar;
                x = [xs xs xs-2] .* array_scalar;
                [x,y] = rotate_array(x,y,array_rotations(a));
                x = x + array_positions(a,1);
                y = y + array_positions(a,2);
                plot(x, y, 'Color', 'k', 'LineWidth', LineWidth*1.1)
            elseif strcmpi(CornerLocation, 'PostMed')
                % Find the posterior-medial corner 
                theta = NaN(4,1);
                for j = 1:4
                    if x_final(j) - array_positions(a,1) > 0 % Assert posterior
                        theta(j) = atan2d(y_final(j) - array_positions(a,2),...
                                          x_final(j) - array_positions(a,1));
                    end
                end
                [~, delta_idx] = min(abs(theta + 46)); % Assert medial
                % Shift
                x_circ = circshift(x_final, -delta_idx-2);
                x_off = x_circ(1:3);
                y_circ = circshift(y_final, -delta_idx-2);
                y_off = y_circ(1:3);
                % Shorten radius to 2 * array_scaling
                for j = [1,3]
                    theta = atan2d(y_off(j) - y_off(2),x_off(j) - x_off(2));
                    xr = (2 * array_scalar * cosd(theta));
                    yr = (2 * array_scalar * sind(theta));
                    x_off(j) = x_off(2) + xr;
                    y_off(j) = y_off(2) + yr;
                end
                % Plot
                plot(x_off, y_off, 'Color', 'k', 'LineWidth', LineWidth*1.1)
            else
                error('"CornerLocation" must be "Original" or "PostMed".')
            end
        end
    end

    % Spline
    plot(xq, yq, 'k', 'LineStyle', SplineStyle, 'LineWidth', LineWidth)
    % Labels
    if FontSize > 0
        text(im_center(1) + im_offset - im_offset*0.1, im_center(2) - im_offset + im_offset*0.1, SubjectString, 'FontSize', FontSize, ...
            'VerticalAlignment','top','HorizontalAlignment','right', 'FontWeight', 'bold', 'Color', SubjectLabelColor,...
            'BackgroundColor', SubjectBackgroundColor, 'Parent', Parent)
        text(motor_label(1), motor_label(2), 'MC', 'FontSize', FontSize, 'VerticalAlignment','middle', ...
            'HorizontalAlignment','center', 'FontWeight', 'bold', 'Color', RegionLabelColor, 'Parent', Parent)
        text(sensory_label(1), sensory_label(2), 'SC', 'FontSize', FontSize, 'VerticalAlignment','middle', ...
            'HorizontalAlignment','center', 'FontWeight', 'bold', 'Color', RegionLabelColor, 'Parent', Parent)
    end
    % Legend
    if ShowScale
        offset = 10 * array_scalar;
        plot(ml_ap_axis_pos(1) - [offset 0 0], ml_ap_axis_pos(2) + [0 0 offset], 'Color', ScaleColor, ...
            'LineWidth', LineWidth, 'Parent', Parent)
        if FontSize > 0
            text(ml_ap_axis_pos(1)-offset, ml_ap_axis_pos(2), 'A', 'FontSize', FontSize,...
                'VerticalAlignment','middle','HorizontalAlignment','right', 'FontWeight','bold', 'Parent', Parent)
            text(ml_ap_axis_pos(1), ml_ap_axis_pos(2)+offset, 'L', 'FontSize', FontSize,...
                'VerticalAlignment','top','HorizontalAlignment','center', 'FontWeight','bold', 'Parent', Parent)
        end
    end
    
    function [xr,yr] = rotate_array(x,y,r)
        r = deg2rad(r);
        xr = x*cos(r) - y*sin(r);
        yr = y*cos(r) + x*sin(r);
    end
    
    function ParseVarargin()
        if ~isempty(varargin)
            nargin = ceil(length(varargin)/2);
            varargin = reshape(varargin, [2, nargin]);
            for n = 1:nargin
                if contains(varargin{1,n}, 'ArrayColor')
                    if  isnumeric(varargin{2,n})
                        ArrayColors = varargin{2,n};
                    elseif ischar(varargin{2,n}) || isstring(varargin{2,n})
                        try 
                            ArrayColors = validatecolor(varargin{2,n});
                        catch
                            error('Color %s could not be validated', varargin{2,n})
                        end
                    else
                        error('"Color" must be a numeric vector of size [1,3], or a char/string.')
                    end
                elseif strcmp(varargin{1,n}, 'SubjectLabelColor')
                    SubjectLabelColor = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'RegionLabelColor')
                    RegionLabelColor = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'SubjectString')
                    SubjectString = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'SplineStyle')
                    SplineStyle = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'LineWidth')
                    LineWidth = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'FontSize')
                    FontSize = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'ScaleColor')
                    ScaleColor = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'ShowScale')
                    ShowScale = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'SubjectBackgroundColor')
                    SubjectBackgroundColor = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'Parent')
                    Parent = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'ShowCorner')
                    ShowCorner = varargin{2,n};
                elseif strcmp(varargin{1,n}, 'CornerLocation')
                    CornerLocation = varargin{2,n};
                else
                    error('%s is an unrecognized input.', varargin{1,n})
                end
            end
        end
    end
end