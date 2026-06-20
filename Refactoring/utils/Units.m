classdef Units < handle
    properties (Constant)
        ft_to_m = 0.3048;       % feet to meters
        m_to_ft = 1/0.3048;     % meters to feet
        lb_to_kg = 0.453592;    % pounds to kilograms
        kg_to_lb = 1/0.453592;  % kilograms to pounds
        ft_s_to_m_s = 0.3048;   % feet per second to meters per second
        m_s_to_ft_s = 1/0.3048; % meters per second to feet per second
    end

    methods (Static)
        function m = ft2m(ft)
            % converts feet to meters
            m = ft * Units.ft_to_m;
        end
        
        function ft = m2ft(m)
            % converts meters to feet
            ft = m * Units.m_to_ft;
        end
        
        function kg = lb2kg(lb)
            % converts pounds to kilograms
            kg = lb * Units.lb_to_kg;
        end
        
        function lb = kg2lb(kg)
            % converts kilograms to pounds
            lb = kg * Units.kg_to_lb;
        end
        
        function m_s = ft_s2m_s(ft_s)
            % converts feet per second to meters per second
            m_s = ft_s * Units.ft_s_to_m_s;
        end
        
        function ft_s = m_s2ft_s(m_s)
            % converts meters per second to feet per second
            ft_s = m_s * Units.m_s_to_ft_s;
        end
    end
end