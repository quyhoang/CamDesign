classdef oskam < kam
    %Oscillating cam
    
    properties
        m_rocker
        rocker2cam
    end
    
    methods
        function obj = oskam(instanceName)
            % Constructor

            obj = obj@kam(instanceName)
        end
        
        function abc(obj)
            % Placeholder for the show method
        end
    end
end
