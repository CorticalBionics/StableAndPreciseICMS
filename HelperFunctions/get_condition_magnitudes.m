function [x,y] = get_condition_magnitudes(table, condition_vector)
    if condition_vector(3) == 0
       condition_vector = condition_vector(1:2);
    end

    cidx = all(table{:,1:length(condition_vector)} == condition_vector,2);
    table = table(cidx,:);

    x = unique([table.JudgeAmp]);
    y = cell([length(x), 1]);
    for j = 1:length(x)
        x_idx = [table.JudgeAmp] == x(j);
        y_temp = cat(1,table.IntAll{x_idx});
        y{j} = y_temp;
    end
    y = cat(2, y{:});
end