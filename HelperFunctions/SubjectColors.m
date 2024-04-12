function c = SubjectColors(subjectID, cb)
    if nargin == 1
        cb = false; % Colorblind alternatives
    end
    
    % Get multiple from ordered cell array
    if length(subjectID) > 1 && isa(subjectID , 'cell') 
        c = zeros(length(subjectID),3);
        for i = 1:length(subjectID)
            c(i,:) = SubjectColors(subjectID{i});
        end
        return
    end

    % Colors specific to subjects for plots
    if strcmpi(subjectID, 'BCI02') % Purple (700)
        if cb
            c = [.49 .44 .70];
        else
            c = [.48 .12 .66];
        end
        
    elseif contains(subjectID, 'CRS02','IgnoreCase',true) % Orange
        if cb
            c = [.85 .37 0];
        else
            c = [.96 .49 0];
        end
         
    elseif contains(subjectID, 'CRS07','IgnoreCase',true) % Green
        if cb
            c = [.11 .62 .47];
        else
            c =  [.22 .56 .24];
        end

    elseif contains(subjectID, 'CRS08','IgnoreCase',true) % Blue
        if cb
            c = [.22 .42 .69];
        else
            c =  [.19 .25 .62];
        end

    elseif strcmpi(subjectID, 'BCI03') % Magenta
        if cb
            c = [.91 .16 .54];
        else
            c = [.76 .09 .36];
        end
    elseif strcmpi(subjectID, 'RP1') % Lime
        if cb
            c = [.91 .16 .54];
        else
            c = [.69 .71 .17];
        end
    else % Backup
        warning('Invalid SubjectID')
        c = [.6 .6 .6];
    end
end