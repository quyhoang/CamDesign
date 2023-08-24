classdef oscillatingCam < Cam
    
    properties
        transition 
        l_roller
        l_load
        
        % Additional subclass properties
        estLoad 
        m_rocker
        rocker2cam
        
    end
    
    methods
        % Implement abstract method
        function calculate(obj)
            % Calculations
        end
        
    end
    
end