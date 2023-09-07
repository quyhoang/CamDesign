classdef tcam < kam
    %Translating cam
    
    properties
        rBase
    end
    
    methods
        function obj = tcam(instanceName)
            % Constructor

            obj = obj@kam(instanceName)
        end
        
        function obj = calculate_roller_position(obj)
            % 

            rPrime = obj.rBase + obj.rRoller;
            obj.roller_position = obj.displacement + rPrime; % cam is at center
        end

    end
end
