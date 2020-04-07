function d = batch_select(param, entry, value)

% return all cases with one param constrained
d = [];
for i = 1:numel(param)
    if(param(i).(entry) == value)
        d = [d;param(i)];
    end
end

end