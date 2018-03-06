function [y,x] = psr_remove_field(y,str)

% Removes field and returns removed field

try
    x = y.(str);
    y = rmfield(y,str);
catch 
    % Assumes that field does not exist
end

end