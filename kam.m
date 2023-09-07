classdef kam  < handle
    %Kam is a general class representing cam mechanism
    % there are children classes of specific cam type
    % handle class will pass instance to method as reference
    
    properties
        transition
        transition_angle
        transition_displacement

        m_roller

        m_load

        rRoller

        RPM
        kFriction

        rollerColor
        pitchColor
        camColor

        maxPressureAngle
        maxPressureAngle_deg
        theta
        T
        time
        timeStep
        sampleRate
        step

        fullStroke = 0
        displacement = 0
        velocity = 0
        acceleration = 0

        roller_position = 0 % center of roller
        pitchCurve = 0
        pitchX = 0
        pitchY = 0
        camSurfX = 0
        camSurfY = 0
    end
    

    methods
        function obj = kam(instanceName)
            % Constructor

            % Construct file paths based on instance name
            transitionFilePath = sprintf('%s_transition.txt', instanceName);
            parameterFilePath = sprintf('%s_parameter.txt', instanceName);

            % Read transition data using readmatrix
            obj.transition = readmatrix(transitionFilePath);
            obj.transition_angle = obj.transition(:, 1)';
            obj.transition_displacement = obj.transition(:, 2)';

            % Call the method to read the parameters
            obj = obj.readParameters(parameterFilePath);

            % Additional calculations
            obj.maxPressureAngle = deg2rad(obj.maxPressureAngle_deg);
            obj.theta = 0:obj.step:360;
            obj.T = 60/obj.RPM;
            obj.time = linspace(0, obj.T, length(obj.theta));
            obj.timeStep = obj.T / size(obj.time, 2);

            obj = obj.calculate_displacement();
            obj = obj.calculate_velocity();
            obj = obj.calculate_acceleration();

            obj = obj.calculate_roller_position();
            obj = obj.calculate_pitchCurve();
            obj = obj.calculate_profile();

            % Validate that all properties have been assigned values
            propertyList = properties(obj);
            for i = 1:length(propertyList)
                if isempty(obj.(propertyList{i}))
                    error('Value not specified for: %s', propertyList{i});
                end
            end
        end
        
        function obj = readParameters(obj, parameterFilePath)
            % Check whether all attributes are read in properly

            % Open the parameter file for reading
            fid = fopen(parameterFilePath, 'r');
            if fid == -1
                error('Could not open the file');
            end

            % Loop through each line in the file to read parameters
            while ~feof(fid)
                % Read and parse the current line
                tline = fgets(fid);
                eqIndex = strfind(tline, '=');
                if ~isempty(eqIndex)
                    varName = strtrim(tline(1:eqIndex-1));
                    varValue = strtrim(tline(eqIndex+1:end));
                    scIndex = strfind(varValue, ';');
                    if ~isempty(scIndex)
                        varValue = strtrim(varValue(1:scIndex-1));
                    end
                    
                    % Convert to number or matrix if possible, otherwise keep as string
                    numericValue = str2num(varValue);
                    if ~isempty(numericValue)
                        varValue = numericValue;
                    end
                    
                    % Assign the parsed value to the respective property of the object
                    obj.(varName) = varValue;
                end
            end
            fclose(fid);
        end
        
        function filtered_pairs = filterConsecutivePairs(obj)
            % manipulate cam transition data
            
            % Add 0 to the beginning and end of transition_displacement
            obj.transition_displacement = horzcat(0, obj.transition_displacement, 0);
        
            % Add 0 and 360 to the beginning and the end of transition_angle
            obj.transition_angle = horzcat(0, obj.transition_angle, 360);
        
            % Check that all elements in transition_angle are unique
            assert(all(diff(obj.transition_angle) > 0) && (length(unique(obj.transition_angle)) == length(obj.transition_angle)), 'The angles must be unique and monotonic.');
        
            % Initialize filtered_pairs as an empty array with the maximum possible size
            filtered_pairs = zeros(length(obj.transition_angle)-1, 2); % avoid dynamic allocation of memory
        
            % Initialize counter
            count = 0;
        
            % Iterate through the angles
            for i = 1:length(obj.transition_angle)-1
                % Check if the current displacement is different from the next one
                if obj.transition_displacement(i) ~= obj.transition_displacement(i+1)
                    % Increment counter
                    count = count + 1;
                    % Add the pair of consecutive angles to filtered_pairs
                    filtered_pairs(count, :) = [obj.transition_angle(i), obj.transition_angle(i+1)];
                end
            end
        
            % Trim filtered_pairs to its actual size
            filtered_pairs = filtered_pairs(1:count, :);
        end
        
        function obj = calculate_displacement(obj)
            % Calculating displacement
            obj.displacement = zeros(size(obj.theta));
            
            % Get the transition points
            filtered_pairs = obj.filterConsecutivePairs();
            
            % Iterate over the transition points
            for i = 1:size(filtered_pairs, 1)
                % Extract the start and end points of the transition
                point = filtered_pairs(i, :);
                h = obj.transition_displacement(obj.transition_angle == point(2)) - obj.transition_displacement(obj.transition_angle == point(1));
                
                if (abs(h) > obj.fullStroke)
                    obj.fullStroke = abs(h);
                end

                bRise = point(2) - point(1);
                
                % Calculate the three sections of the displacement curve
                tempTheta1 = obj.theta(obj.theta >= point(1) & obj.theta < point(1) + bRise/8) - point(1);
                sRise1 = h/(4+pi)*(pi*tempTheta1/bRise - 1/4*sin(4*pi*tempTheta1/bRise));
                
                tempTheta2 = obj.theta(obj.theta >= point(1) + bRise/8 & obj.theta < point(1) + 7*bRise/8) - point(1);
                sRise2 = h/(4+pi)*(2 + pi*tempTheta2/bRise - 9/4*sin(pi/3 + 4*pi/3*tempTheta2/bRise));
                
                tempTheta3 = obj.theta(obj.theta >= point(1) + 7*bRise/8 & obj.theta <= point(2)) - point(1);
                sRise3 = h/(4+pi)*(4 + pi*tempTheta3/bRise - 1/4*sin(4*pi*tempTheta3/bRise));
                
                % Combine all parts and add the previous displacement value
                sRise = [sRise1, sRise2, sRise3];
                
                % Update the displacement array
                obj.displacement(obj.theta >= point(1) & obj.theta <= point(2)) = obj.displacement(obj.theta >= point(1) & obj.theta <= point(2)) + sRise;
                obj.displacement(obj.theta > point(2)) = obj.displacement(obj.theta > point(2)) + sRise(end);
                % Update the cumulative displacement for the next period
            end
        end
        
        function pos(obj)
            % Plot position vs angle in cartesian coordinate

            figure;
            plot(obj.theta,obj.displacement);
            grid on;
            grid minor;
            xlim([0 360]);
            xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
            ylim([-2*obj.fullStroke + obj.fullStroke/2 2*obj.fullStroke + obj.fullStroke/2]);
            ylabel({'位置','mm'},'FontSize',15,'FontWeight','light','Color','b');
            [tit,] = title({'';'S Diagram'},{['モーター回転速度 ',num2str(obj.RPM),'rpm   ','T = ', num2str(obj.T),'s'];''},...
                'Color','blue');
            tit.FontSize = 15;
        end

        function obj = calculate_velocity(obj)
            % velocity with respect to time
            
            obj.velocity = diff(obj.displacement)/obj.timeStep;
            obj.velocity = [obj.velocity obj.displacement(1)-obj.displacement(length(obj.displacement))]; 
            %add the last element to make the length of vv and theta equal
        end

        function velo(obj)
            % Plot velocity vs angle in cartesian coordinate

            figure;
            plot(obj.theta,obj.velocity);
            grid on;
            grid minor;
            xlim([0 360]);
            xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
            ylabel({'速度','mm/s'},'FontSize',15,'FontWeight','light','Color','b');
            tempV = strcat('最大速度 ',num2str(max(obj.velocity)),' mm/s');
            [tit,] = title({'';'V Diagram';tempV},{['モーター回転速度 ',num2str(obj.RPM),'rpm   ','T = ', num2str(obj.T),'s'];''},...
                'Color','blue');
           tit.FontSize = 15;
        end

        function obj = calculate_acceleration(obj)
            % acceleration with respect to time
            
            obj.acceleration = diff(obj.velocity)/obj.timeStep;
            obj.acceleration = [obj.acceleration obj.velocity(1)-obj.velocity(length(obj.velocity))]/1000;
        end

        function accel(obj)
            % Plot acceleration vs angle in cartesian coordinate

            figure;
            plot(obj.theta,obj.acceleration);
            grid on;
            grid minor;
            xlim([0 360]);
            xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
            ylabel({'加速','m/s^2'},'FontSize',15,'FontWeight','light','Color','b');
            tempA = strcat('最大加速  ', num2str(max(obj.acceleration)),' m/s^2');
            [tit,] = title({'';'A Diagram';tempA},{['モーター回転速度 ',num2str(obj.RPM),'rpm   ','T = ', num2str(obj.T),'s'];''},...
                'Color','blue');
            tit.FontSize = 15;
        end 

        function obj = calculate_pitchCurve(obj)
            thetaRadian = deg2rad(obj.theta);
            obj.pitchCurve = obj.roller_position.*exp(thetaRadian*1i);
            obj.pitchX = real(obj.pitchCurve); 
            obj.pitchY = imag(obj.pitchCurve);
        end
   
        function pitch(obj)
            % Show pitch curve
            figure;
            plot(obj.pitchX,obj.pitchY,'color',obj.pitchColor)
            axis equal; grid on; grid minor;
        end
        
        function obj = calculate_profile(obj)
            [obj.camSurfX,obj.camSurfY] = offsetIn(obj.pitchX,obj.pitchY,obj.rRoller);
        end

        function profile(obj)
            % Show cam profile
            figure;
            plot(obj.camSurfX,obj.camSurfY,'color',obj.camColor)
            axis equal; grid on; grid minor;
        end



    end

    methods (Abstract)
        calculate_roller_position(obj) % Abstract method, to be defined by each child class
    end
end
