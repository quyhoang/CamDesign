classdef kam  < handle
    %Kam is a general class representing cam mechanism
    % there are children classes of specific cam type
    % handle class will pass instance to method as reference

    properties
        camType % ocam (oscillating cam) or tcam (translating cam)
        min_curvature = 0 % Minimum radius of curvature
        rRoller % Radius of the roller
        max_pressureAngle = 0 % Maximum pressure angle
        allowedPressureAngle_deg % Max pressure angle in degree
        RPM % Motor velocity in rounds per minutes
        transition % Input transition matrix
        transition_angle % Angle values extracted from the transition matrix
        transition_displacement % Displacement values extracted from the transition matrix

        m_roller % Roller mass in kilograms
        m_load % Load mass in kilograms


        rPrime = 0 % Prime radius (roller radius plus base radius)

        rollerColor % Color to use for the roller in graphics
        pitchColor % Color to use for the pitch in graphics
        camColor % Color to use for the cam in graphics
        angleColor % to use in graphics

        allowedPressureAngle % Max pressure angle in radian

        theta % Range of angles from 0 to 360 with step size as defined by 'step'
        T % Period of moving 360 degrees, in seconds
        time % Array of time values corresponding to each angle in 'theta'
        timeStep % Convert step in degree to step in time
        sampleRate % For animation
        step % step size of theta

        fullStroke = 0
        displacement = camdata('ローラー変位','mm') % follower (roller) displacement
        load_displacement = camdata('部品変位','mm') % load displacement

        velocity = camdata('速度', 'mm/s')
        acceleration = camdata('加速','m/s2')

        roller_position = 0 % center of roller
        pitchCurve = 0
        pitchX = 0 % x coordinate of pitch curve
        pitchY = 0 % y coordinate of pitch curve
        camSurfX = 0 % x coordinate of cam profile
        camSurfY = 0 % y coordinate of cam profile

        curvature = camdata('曲率半径','mm') % Radius of curvature
        pressureAngle = camdata('圧力角','°')
        motorTorque = camdata('モータートルク','N/m')

        initialSpringDisplacement % initial spring displacement
        minK = 0 % minimum acceptable spring hardness
        springK % actual spring hardness

        fSpring = 0 % Force exerted by spring. Positive direction is from cam center outward
        objlength = 0 % number of sampled points over 360 deg period

        fCutMax = 0 % Maximum cutting force (N)
        cutThickness = 0 %mm
        cutClearance = 0 %mm
        backClampForce = 0 %N
        fCut = 0 % Force needed to cut a workpiece

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

            obj.allowedPressureAngle = deg2rad(obj.allowedPressureAngle_deg);
            obj.theta = 0:obj.step:360;
            obj.objlength = size(obj.theta);
            obj.T = 60/obj.RPM;
            obj.time = linspace(0, obj.T, length(obj.theta));
            obj.timeStep = obj.T / size(obj.time, 2);

            if strcmpi(obj.camType, 'tcam') % Calculate rPrime is the child class has rBase attribute
                obj.rPrime = obj.rBase + obj.rRoller;
            end

            % Calculating cam dynamic properties
            obj = obj.calculate_displacement(); % Follower displacement
            obj = obj.calculate_load_displacement(); % Load displacement
            obj = obj.calculate_velocity(); 
            obj = obj.calculate_acceleration();

            obj = obj.calculate_roller_position();
            obj = obj.calculate_pitchCurve();
            obj = obj.calculate_profile();
            obj = obj.calculate_curvature();
            obj = obj.calculate_roller_position();
            obj = obj.calculate_pressure_angle();


            % Validate that all properties have been assigned values
            propertyList = properties(obj);
            for i = 1:length(propertyList)
                if isempty(obj.(propertyList{i}))
                    error('Value not specified for: %s', propertyList{i});
                end
            end

            % Check that conditions for the CAM are met
            assert(obj.max_pressureAngle < obj.allowedPressureAngle_deg, 'Max pressure angle is too large.');
            assert(obj.min_curvature > obj.rRoller, 'Minimum radius of curvature is too small.');
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

        function s(obj)
            % Plot position vs angle in cartesian coordinate

            figure;
            plot(obj.theta,obj.load_displacement.data);
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
            % velocity.data with respect to time

            obj.velocity.data = diff(obj.displacement.data)/obj.timeStep;
            obj.velocity.data = [obj.velocity.data obj.displacement.data(1)-obj.displacement.data(length(obj.displacement.data))];
            %add the last element to make the length of vv and theta equal
        end

        function v(obj)
            % Plot velocity.data vs angle in cartesian coordinate

            figure;
            plot(obj.theta,obj.velocity.data);
            grid on;
            grid minor;
            xlim([0 360]);
            xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
            ylabel({'速度','mm/s'},'FontSize',15,'FontWeight','light','Color','b');
            tempV = strcat('最大速度 ',num2str(max(obj.velocity.data)),' mm/s');
            [tit,] = title({'';'V Diagram';tempV},{['モーター回転速度 ',num2str(obj.RPM),'rpm   ','T = ', num2str(obj.T),'s'];''},...
                'Color','blue');
            tit.FontSize = 15;
        end

        function obj = calculate_acceleration(obj)
            % acceleration.data with respect to time

            obj.acceleration.data = diff(obj.velocity.data)/obj.timeStep;
            obj.acceleration.data = [obj.acceleration.data obj.velocity.data(1)-obj.velocity.data(length(obj.velocity.data))]/1000;
        end

        function a(obj)
            % Plot acceleration.data vs angle in cartesian coordinate

            figure;
            plot(obj.theta,obj.acceleration.data);
            grid on;
            grid minor;
            xlim([0 360]);
            xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
            ylabel({'加速','m/s^2'},'FontSize',15,'FontWeight','light','Color','b');
            tempA = strcat('最大加速  ', num2str(max(obj.acceleration.data)),' m/s^2');
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

        function export(obj)
            % Export profile data as txt file
            % Save both camProfile.txt and camDirection.txt files to Creo
            % working directory and then type "gcam" or click the カム作成
            % icon, the 3D model of the cam will be created.

            x_cord = transpose(obj.camSurfX);
            y_cord = transpose(obj.camSurfY);
            z_cord = zeros(length(obj.theta),1);

            camProfile = [x_cord y_cord z_cord];
            writematrix(camProfile,'camProfile.txt','Delimiter','tab');

            % Show cam motion direction
            indices = 2.^(0:log2(length(x_cord)));
            x_direction = x_cord(indices);
            y_direction = y_cord(indices);
            z_direction = ones(size(x_direction));
            camDirection = [x_direction y_direction z_direction];
            writematrix(camDirection,'camDirection.txt','Delimiter','tab');

            str = "CreoAutomation.exe がアクティブで、" + ...
                "txt データが Creo 作業ディレクトリに保存されている場合、" + ...
                "「gcam」を押すか、UI からか　カムの 3D モデルが作成されます。";
            disp(str)
        end

        function machining(obj)
            % Show pitch curve, cutting tool, and cam profile

            splRate = round(obj.sampleRate*length(obj.theta)/360);
            x_sample = transpose(obj.pitchX(1:splRate:length(obj.pitchX)));
            y_sample = transpose(obj.pitchY(1:splRate:length(obj.pitchY)));

            figure;

            plot(obj.pitchX,obj.pitchY,'color',obj.pitchColor)
            hold on
            plot(obj.camSurfX,obj.camSurfY,'color',obj.camColor)
            hold on
            axis equal; grid on;

            for k = 2:1:length(x_sample)
                viscircles([x_sample(k),y_sample(k)],obj.rRoller,'LineWidth',1,'Color',obj.rollerColor);
                xlim([min(obj.pitchX)-obj.rRoller*2 max(obj.pitchX)+obj.rRoller*2]);
                ylim([min(obj.pitchY)-obj.rRoller*2 max(obj.pitchY)+obj.rRoller*2]);
                axis equal; grid on;
                drawnow
            end
            plot(obj.camSurfX,obj.camSurfY,'color',obj.camColor)
            xlim([min(obj.pitchX)-obj.rRoller*2 max(obj.pitchX)+obj.rRoller*2]);
            ylim([min(obj.pitchY)-obj.rRoller*2 max(obj.pitchY)+obj.rRoller*2]);
            axis equal; grid on; grid minor;
        end

        function obj = calculate_curvature(obj)
            obj.curvature.data = zeros(size(obj.camSurfX));
            L = length(obj.camSurfX);
            % Boundary. Note that the first and the last points on profile curve are the
            % same
            obj.curvature.data(1) = circumscribedR([obj.camSurfX(L-1) obj.camSurfX(1) obj.camSurfX(2)],[obj.camSurfY(L-1) obj.camSurfY(1) obj.camSurfY(2)]);
            obj.curvature.data(L) = circumscribedR([obj.camSurfX(L-1) obj.camSurfX(L) obj.camSurfX(2)],[obj.camSurfY(L-1) obj.camSurfY(L) obj.camSurfY(2)]);

            for k = 1:1:L-2
                X = obj.camSurfX(k:1:k+2);
                Y = obj.camSurfY(k:1:k+2);
                obj.curvature.data(k+1) = circumscribedR(X,Y);
            end

            obj.min_curvature = min(obj.curvature.data);
        end

        function r(obj)
            % Show radius of curvature.data

            figure;
            yyaxis left
            semilogy(obj.theta, obj.curvature.data,'Color',obj.angleColor);
            ax = gca;
            ax.YColor = obj.angleColor;
            grid on;
            grid minor;
            xlim([0 360]);
            xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color',obj.angleColor);
            ylabel({'曲率半径','mm'},'FontSize',15,'FontWeight','light','Color',obj.angleColor);

            yyaxis right
            strokeColor = [0.6350 0.0780 0.1840];
            plot(obj.theta,obj.displacement.data,'Color',strokeColor);
            ax = gca;
            ax.YColor = strokeColor;
            grid on;
            grid minor;
            xlim([0 360]);
            ylabel({'位置','mm'},'FontSize',15,'FontWeight','light','Color',strokeColor);
            hold on

            tempP = strcat('最小曲率半径 ',num2str(min(obj.curvature.data)),'mm');
            title({'';'曲率半径・位置　vs　回転角度'; tempP; ''},'Color','b','FontSize',15,'FontWeight','light');
            disp(strcat('最小曲率半径: ',num2str(min(obj.curvature.data)),'mm'));
        end

        function pressure(obj)
            % Show pressure angle

            figure;
            yyaxis left
            plot(obj.theta, obj.pressureAngle.data,'Color',obj.angleColor);
            ax = gca;
            ax.YColor = obj.angleColor;
            grid on;
            grid minor;
            xlim([0 360]);
            xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color',obj.angleColor);
            ylabel({'圧角','degree'},'FontSize',15,'FontWeight','light','Color',obj.angleColor);

            yyaxis right
            strokeColor = [0.6350 0.0780 0.1840];
            plot(obj.theta,obj.displacement.data,'Color',strokeColor);
            ax = gca;
            ax.YColor = strokeColor;
            grid on;
            grid minor;
            xlim([0 360]);
            ylabel({'位置','mm'},'FontSize',15,'FontWeight','light','Color',strokeColor);
            hold on

            tempP = strcat('最大圧角 ',num2str(max(obj.pressureAngle.data)),'^o');
            %  title(temp,'Color','b','FontSize',15,'FontWeight','light');

            title({'';'圧角・位置　vs　回転角度'; tempP; ''},'Color','b','FontSize',15,'FontWeight','light');

            disp(strcat('最大圧角: ',num2str(max(obj.pressureAngle.data)),'度'));
        end

        function spring(obj)
            % Plot spring force vs angle in cartesian coordinate

            figure;
            plot(obj.theta,obj.fSpring);
            grid on;
            grid minor;
            xlim([0 360]);
            xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
            ylabel({'バネ力','N'},'FontSize',15,'FontWeight','light','Color','b');
            tempA = strcat('最大力  ', num2str(min(obj.fSpring)),' N');
            [tit,] = title({'';'バネ力 vs 回転角度';tempA},{['モーター回転速度 ',num2str(obj.RPM),'rpm   ','T = ', num2str(obj.T),'s'];''},...
                'Color','blue');
            tit.FontSize = 15;
        end

        function obj = calculate_cutting_force(obj)
            if isempty(obj.fCutMax) || (obj.fCutMax <= 0)
                % disp("切断操作なし。");
                obj.fCut = 0;
            else
                cutStartAngleIndex = find(obj.load_displacement.data > obj.cutClearance, 1) - 1;
                % find the first index at which displacement bigger than rBase,
                % which means after cutting process start, minus 1 so that the
                % index become that of the position just before cutting
                cutEndAngleIndex = find(obj.load_displacement.data >= obj.cutClearance + obj.cutThickness,1);
                % the index at which cutting process ends
                obj.fCut = zeros(obj.objlength);
                % Initialize a vector of size equal to that of displacement with all elements set to 0
                obj.fCut(cutStartAngleIndex:cutEndAngleIndex) = obj.fCutMax + obj.backClampForce;
                % Set the elements corresponding to position during cutting to max cutting force
            end
        end

        function torque(obj)
            % Show motor torque
            
            figure;
            plot(obj.theta,obj.motorTorque.data);
            xlabel('角度') ; 
            xlim([0 360]);
            ylabel('トルク (Nm)') ;
            torqueTitle = strcat('最大トルク  ', num2str(max(obj.motorTorque.data)),' Nm');
            title(torqueTitle,'Color','b','FontSize',15,'FontWeight','light');
            grid on; grid minor;
            disp(torqueTitle);
        end
    end

    methods (Abstract) % Abstract method, to be defined by each child class
        calculate_displacement(obj) % roller displacement, not load
        % load displacement is described by transition matrix

        calculate_load_displacement(obj)

        calculate_roller_position(obj)

        calculate_pressure_angle(obj)

        animation(obj)

        calculate_spring_force(obj)
        % Calculate the minimum force the spring must exert to keep the
        % cam and the follower always in contact
   
        calculate_motor_torque(obj)
    end

end














