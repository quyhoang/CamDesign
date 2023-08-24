% CamDesign class
classdef cam
    properties
        transition % Transition points
        l_roller % Distance from arm rotating axis to roller center
        l_load % Distance from arm center to load
        estLoad % Estimated load
        m_roller % Roller mass
        rRoller % Roller radius
        m_rocker % Rocker arm mass
        rocker2cam % Distance between rocker arm and cam axes
        RPM % Motor velocity
        maxPressureAngle_deg % Maximum pressure angle
        kFriction % Friction coefficient
        sampleRate % Sample rate for pitch curve
        step % Step size for calculations
        
        % Results
        displacement 
        velocity
        acceleration
        angularDisplacement
        curvature
        pressureAngle
        motorTorque
        
    end
    
    methods
        function obj = cam(filename)
            if nargin > 0
                fid = fopen(filename, 'r');
                while ~feof(fid)
                    line = fgetl(fid);
                    [key, value] = strtok(line, '=');
                    key = strtrim(key);
                    value = str2double(strtrim(value(2:end)));
                    obj.(key) = value;
                end
                fclose(fid);
            end
        end
        
        function calculate(obj)
            % Perform calculations
            
            % Preliminary calculations
            theta = 0:obj.step:360;
            T = 60/obj.RPM; 
            time = linspace(0,T,length(theta));
            timeStep = T/length(time);
            
            % Displacement
            obj.displacement = obj.calcDisplacement(theta);
            
            % Velocity
            obj.velocity = obj.calcVelocity(timeStep);
            
            % Acceleration 
            obj.acceleration = obj.calcAcceleration(timeStep);
            
            % Angular displacement
            obj.angularDisplacement = obj.calcAngularDisplacement();
            
            % Curvature
            obj.curvature = obj.calcCurvature();
            
            % Pressure angle
            obj.pressureAngle = obj.calcPressureAngle();
            
            % Motor torque
            obj.motorTorque = obj.calcMotorTorque();
            
        end
                
        % Other calculation methods
        function disp = calcDisplacement(obj,theta)
           % Calculate displacement
        end
        
        function vel = calcVelocity(obj,timeStep)
            % Calculate velocity
        end
        
        function acc = calcAcceleration(obj,timeStep)
            % Calculate acceleration
        end
        
        function angDisp = calcAngularDisplacement(obj)
            % Calculate angular displacement
        end
        
        function curv = calcCurvature(obj)
            % Calculate curvature
        end
        
        function pAngle = calcPressureAngle(obj)
            % Calculate pressure angle
        end
        
        function mTorque = calcMotorTorque(obj)
            % Calculate motor torque
        end
        
    end
end

% Usage
transition = [...]; 
cam = CamDesign(transition, l_roller, l_load, ..., step);
cam.calculate();
disp = cam.displacement;