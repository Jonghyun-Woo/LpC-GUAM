classdef SurfaceDynamics < handle
    properties
        surfaces

        rate_limit
        bandwidth

        ail_max
        ail_min
        flp_max
        flp_min
        ele_max
        ele_min
        rud_max
        rud_min
    end

    methods 
        function obj = SurfaceDynamics()
            obj.ail_max = rad2deg(30);
            obj.ail_min = rad2deg(-30);
            obj.flp_max = rad2deg(30);
            obj.flp_min = rad2deg(-30);
            obj.ele_max = rad2deg(30);
            obj.ele_min = rad2deg(-30);
            obj.rud_max = rad2deg(30);
            obj.rud_min = rad2deg(-30);
        end

    end
end