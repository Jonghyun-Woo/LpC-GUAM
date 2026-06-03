classdef Units < handle
    properties (Constant)
        ft_to_m = 0.3048;       % feet to meters
        m_to_ft = 1/0.3048;     % meters to feet
        lb_to_kg = 0.453592;    % pounds to kilograms
        kg_to_lb = 1/0.453592;  % kilograms to pounds
        ft_s_to_m_s = 0.3048;   % feet per second to meters per second
        m_s_to_ft_s = 1/0.3048; % meters per second to feet per second
    end
end