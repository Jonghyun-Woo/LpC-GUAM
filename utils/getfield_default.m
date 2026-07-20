function v = getfield_default(s, f, d)
% GETFIELD_DEFAULT Return s.(f) if present, otherwise the default d.
%   Used by config constructors to parse optional 'overrides' structs.
    if isstruct(s) && isfield(s, f)
        v = s.(f);
    else
        v = d;
    end
end
