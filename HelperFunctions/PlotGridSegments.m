function PlotGridSegments(grid_locations, grid_tags, segmenter, text_size, axis_handle, rotation)
    if nargin < 6
        rotation = 0;
    end
    x_offset = size(grid_locations,2)/2 + 1;
    y_offset = size(grid_locations,1)/2 + 1;
    tag_list = {segmenter.Tag};
    for i = 1:numel(grid_locations)
        [gy,gx] = ind2sub(size(grid_locations), i); % Get grid position
        x = [gx, gx, gx+1, gx+1]-x_offset;
        y = [gy, gy+1, gy+1, gy]-y_offset;
        if rotation ~= 0
            [x,y] = rotate_grid(x,y,rotation);
        end
        if isnan(grid_locations(i)) % Show empty channels as white
            patch(x, y, [1 1 1], ...
                'EdgeColor', [.6 .6 .6], 'Parent', axis_handle)
        else
            if isempty(grid_tags{i}) % Indicate no RF with gray
                patch(x, y, [.6 .6 .6], ...
                    'EdgeColor', [.6 .6 .6], 'Parent', axis_handle)
            else
                % Check the grid tag
                if isa(grid_tags{i}, 'char') || (isa(grid_tags{i}, 'cell') && length(grid_tags{i}) == 1)
                    t_idx = find(strcmp(tag_list, grid_tags{i}));
                    patch(x, y, segmenter(t_idx).Color, ...
                        'EdgeColor', [.6 .6 .6], 'Parent', axis_handle, 'FaceAlpha', 1) %#ok<*FNDSB> 
                else
                    u_tags = unique(grid_tags{i});
                    num_tags = length(u_tags);
                    x0 = 0;
                    x1 = 1/num_tags;
                    for t = 1:num_tags
                        t_idx = find(strcmp(tag_list, u_tags{t}));
                        xt = [gx+x0, gx+x0, gx+x1, gx+x1]-x_offset;
                        yt = [gy, gy+1, gy+1, gy]-y_offset;
                        if rotation ~= 0
                            [xt,yt] = rotate_grid(xt,yt,rotation);
                        end
                        patch(xt, yt, segmenter(t_idx).Color, ...
                            'EdgeColor', [.6 .6 .6], 'Parent', axis_handle, 'FaceAlpha', 1)
                        x0 = x0 + 1/num_tags;
                        x1 = x1 + 1/num_tags;
                    end
                end
            end
            if text_size > 0
                if rotation ~= 0
                    [xt,yt] = rotate_grid(gx+.5-x_offset, gy+.5-y_offset,rotation);
                end
                text(xt, yt, num2str(grid_locations(i)),...
                    'Color', 'w', 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontSize', text_size)
            end
        end
    end
    % Better bounding box
    x = [1 , size(grid_locations,2) + 1, size(grid_locations,2) + 1, 1]-x_offset;
    y = [1, 1, size(grid_locations,1) + 1, size(grid_locations,1) + 1]-y_offset;
    if rotation ~= 0
        [x,y] = rotate_grid(x,y,rotation);
    end
    patch(x, y, 'w', 'EdgeColor', 'k', 'FaceAlpha', 0, 'Parent', axis_handle, 'LineWidth', 1)
end

function [xr,yr] = rotate_grid(x,y,r)
    r = deg2rad(r);
    xr = x*cos(r) - y*sin(r); % only care about one axis
    yr = y*cos(r) + x*sin(r);
end