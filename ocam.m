classdef ocam < kam
    %Oscillating cam
    
    properties
        m_rocker
        rocker2cam

        l_roller
        l_load

        initial_angular_displacement
        angular_displacement = 0;
    end
    
    methods
        function obj = ocam(instanceName)
            % Constructor

            obj = obj@kam(instanceName)
        end
        
        function obj = calculate_roller_position(obj)
            % 

            s2rad = obj.displacement/obj.l_load; % convert arc length to angular displacement
            s_rad_initial = deg2rad(obj.initial_angular_displacement);
            s2rad = s2rad + s_rad_initial; % angular displacement

            roller_position1 = obj.l_roller*exp(s2rad*1i);  %unregulated position, rocker is at center
            obj.roller_position = roller_position1 - obj.rocker2cam; % cam is at center, rocker axis is at (-rocker2cam,0)
        end
        
    end
end
