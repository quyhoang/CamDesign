classdef ocam < kam
    %Oscillating cam

    properties
        m_rocker
        rocker2cam

        l_roller
        l_load

        initial_angular_displacement
        angular_displacement = 0;
        angular_acceleration = 0;

        s2rad = 0;
        roller_position_origin = 0 

        l_fSpring % distance from spring applied force to rocker rotation center

        initial_rocker_angle % initial angle from rocker arm to the line perpendicular to action line
        % counter clockwise is positive. This initial angle is from -180 to
        % 180 degree. Note that the cam also rotate counterclockwise.

        inertialMoment = 0
        frictionTorque = 0

    end

    methods
        function obj = ocam(instanceName)
            % Constructor

            obj = obj@kam(instanceName)
            obj.initial_rocker_angle = deg2rad(obj.initial_rocker_angle);
            obj = obj.calculate_angular_acceleration();
            obj = obj.calculate_spring_force();
            obj = obj.calculate_cutting_force();
            obj = obj.calculate_motor_torque();
        end

        function obj = regulate_transition_displacement(obj)
            % convert from displacement of load to arc displacement of follower

            initial_roller_arc_displacement = obj.l_roller * obj.initial_rocker_angle;
            % distance from initial position to projection of rocker
            % rotation center to action line
            initial_load_position = obj.l_load * sin(obj.initial_rocker_angle);
            respective_rocker_angle = asin((obj.transition_displacement+initial_load_position)/obj.l_load);
            % note that radius of roller at load side does not affect this

            obj.transition_displacement = obj.l_roller * respective_rocker_angle - initial_roller_arc_displacement;
        end
        
        function obj = calculate_displacement(obj)
            % Follower arc displacement

            obj = obj.regulate_transition_displacement();

            % Calculating displacement of oscillating cam
            obj.displacement = zeros(obj.objlength);

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
            % obj.angular_displacement = obj.displacement/obj.l_roller;
            % obj.load_displacement = obj.l_load * obj.angular_displacement;
        end        

        function obj = calculate_angular_displacement(obj)
            obj.angular_displacement = obj.displacement/obj.l_roller; %in radian
        end

        function obj = calculate_load_displacement(obj)
            obj = obj.calculate_angular_displacement();
            obj.load_displacement.data = obj.l_load * obj.angular_displacement;
        end

        function obj = calculate_angular_acceleration(obj)
            % angular acceleration with respect to time

            obj.angular_acceleration = obj.acceleration.data /(obj.l_roller/1000);
        end
            
        function obj = calculate_roller_position(obj)
            % Roller position in Cartesian coordinate

            obj.s2rad = obj.displacement/obj.l_roller + deg2rad(obj.initial_angular_displacement); % angular displacement
            obj.roller_position = obj.l_roller*exp(obj.s2rad*1i) - obj.rocker2cam; % cam is at center, rocker axis is at (-rocker2cam,0)
        end

        function obj = calculate_pressure_angle(obj)
            %

            thetaRadian = deg2rad(obj.theta);

            rollerCenterX = real(obj.roller_position);
            rollerCenterY = imag(obj.roller_position);

            rockerNormalAngle = rad2deg(angle(obj.roller_position_origin))+90;

            normalPhase = zeros(size(obj.theta));
            contactPointonCamPhase = zeros(size(obj.theta));
            contactPoint2CamDistance = zeros(size(obj.theta));

            camCenter = [0 0];

            for i = 1:length(obj.theta)
                j = thetaRadian(i);

                % Update roller center position
                tempRollerCenter = [rollerCenterX(i) rollerCenterY(i)];

                % Update cam-roller contact point
                rotatedCam = rotateCw([obj.camSurfX;obj.camSurfY],j);
                contactPoint = [rotatedCam(1,i) rotatedCam(2,i)];

                normalPhase(i) = segmentPhase(contactPoint,tempRollerCenter);
                contactPointonCamPhase(i) = segmentPhase(camCenter,contactPoint);
                contactPoint2CamDistance(i) = norm(contactPoint); %distance from cam center to contact point
            end

            obj.pressureAngle.data = (rockerNormalAngle - normalPhase);
            obj.max_pressureAngle = max(obj.pressureAngle.data);
        end

        function obj = animation(obj)
            %

            % Draw rocker arm
            armCenterX = -obj.rocker2cam;
            armCenterY = 0;


            % Draw roller
            index = linspace(0,2*pi,100);

            plot(0,0,'o','MarkerFaceColor','b');
            axis equal; grid on; grid minor;


            maxDim = obj.rocker2cam+5;
            thetaRadian = deg2rad(obj.theta);

            rollerCenterX = real(obj.roller_position);
            rollerCenterY = imag(obj.roller_position);

            for i = 1:length(obj.theta)
                j = thetaRadian(i);
                hold on
                grid on;
                grid minor;
                axis equal;
                xlim([-maxDim maxDim]);
                ylim([-maxDim maxDim]);

                plot(armCenterX,armCenterY,'o','MarkerFaceColor','r');
                plot(0,0,'o','MarkerFaceColor','r');

                temp6 = strcat('曲率半径　',num2str(obj.curvature.data(i)),' mm     ');
                temp5 = strcat('圧角　',num2str(obj.pressureAngle.data(i)),'^o     ');

                temp3 = strcat('回転角度　',num2str(obj.theta(i)),'^o   ');
                updatedTitle = {temp3; temp5; temp6};
                title(updatedTitle,'Color',[0 0.4470 0.7410],'FontSize',14);

                temp4 = strcat('変位　',num2str(obj.displacement(i)),' mm     ');
                ylabel(temp4,'Color',obj.angleColor,'FontSize',15);
                temp1 = strcat('経過時間　',num2str(obj.time(i)),' s     ');
                xlabel(temp1,'Color',obj.angleColor,'FontSize',15);

                % Update rocker arm
                armX = [armCenterX,rollerCenterX(i)];
                armY = [armCenterY,rollerCenterY(i)];
                plot(armX,armY); hold on

                % Update roller position
                xC = obj.rRoller*cos(index) + rollerCenterX(i);
                yC = obj.rRoller*sin(index) + rollerCenterY(i);
                plot(xC,yC,'color',obj.rollerColor); hold on

                % Update roller center position
                tempRollerCenterX = rollerCenterX(i);
                tempRollerCenterY = rollerCenterY(i);
                plot(tempRollerCenterX,tempRollerCenterY,'o','MarkerFaceColor',[0 0.4470 0.7410]); hold on

                % Update rotated pitch
                rotatedPitch = rotateCw([obj.pitchX;obj.pitchY],j);
                xx = rotatedPitch(1,:);
                yy = rotatedPitch(2,:);
                plot(xx,yy,'color',obj.pitchColor);hold on;

                % Update rotated cam
                rotatedCam = rotateCw([obj.camSurfX;obj.camSurfY],j);
                xx2 = rotatedCam(1,:);
                yy2 = rotatedCam(2,:);
                plot(xx2,yy2,'color',obj.camColor); hold on;

                % Update cam-roller contact point
                % rotatedCam = rotateCw([camSurfX;camSurfY],j);
                contactPointX = rotatedCam(1,i);
                contactPointY = rotatedCam(2,i);
                plot(contactPointX, contactPointY,'.','color','r'); hold on;
                pause(0.005)
                clf
            end

            hold on
            grid on;
            grid minor;
            axis equal;
            xlim([-maxDim maxDim]);
            ylim([-maxDim maxDim]);

            plot(armCenterX,armCenterY,'o','MarkerFaceColor','r');
            plot(0,0,'o','MarkerFaceColor','r'); hold on
            plot(armX,armY); hold on
            plot(xC,yC,'color',obj.rollerColor); hold on
            plot(tempRollerCenterX,tempRollerCenterY,'o','MarkerFaceColor',[0 0.4470 0.7410]); hold on
            plot(xx,yy,'color',obj.pitchColor);hold on;
            plot(xx2,yy2,'color',obj.camColor); hold on;
            plot(contactPointX, contactPointY,'.','color','r'); hold on;

        end

        function obj = calculate_spring_force(obj)
            % Find the momentum needed for follower to always remains in
            % contact with cam.
            % k is spring constant in N/mm
            % initialDisplacement is initial displacement in mm

            I_roller = obj.m_roller*(obj.l_roller/1000)^2;
            I_rocker = 1/3*obj.m_rocker*(obj.l_load/1000)^2;
            I_load = obj.m_load*(obj.l_load/1000)^2;
            obj.inertialMoment = I_roller+I_rocker+I_load;
            
            negativeAcceleration = obj.angular_acceleration;
            negativeAcceleration(obj.angular_acceleration > 0) = 0;
            minimum_spring_torque = obj.inertialMoment*negativeAcceleration - obj.frictionTorque;
            
            % approximation, assume that spring force is almost always
            % perpendicular to rocker
            minimum_spring_force = minimum_spring_torque/(obj.l_fSpring/1000);
            % figure; plot(minimum_spring_force)

            springDisplacement = obj.initialSpringDisplacement + obj.angular_displacement*obj.l_fSpring;
            % Check if all elements are positive using assert
            assert(all(springDisplacement > 0), 'The spring must always be compressed. Please change increase initial spring displacement.');

            % Find the smallest value of k so that fSpring is always
            % more negative than minimum_spring_force
            obj.minK = max(-1 * minimum_spring_force ./ springDisplacement);
            if (obj.springK < obj.minK)
                obj.springK = obj.minK;
                disp(['Spring Hook modulus is reset to minimum allowable value: ', num2str(obj.minK)]);
            end
            obj.fSpring = -1 * obj.springK * springDisplacement;
        end

        function obj = calculate_motor_torque(obj)
            angular_position = deg2rad(obj.initial_rocker_angle) + obj.angular_displacement;
            cutTorque = obj.fCut .* cos(angular_position) * obj.l_load/1000;
            inertialTorque = obj.inertialMoment * obj.angular_acceleration;
            springTorque = obj.fSpring * obj.l_fSpring/1000; % in opposite direction of motor torque
            obj.motorTorque.data = inertialTorque + cutTorque - springTorque;
        end

    end
end
