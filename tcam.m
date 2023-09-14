classdef tcam < kam
    %Translating cam

    properties
        rBase % Base radius
        fFriction % Maximum static friction force
    end

    methods
        function obj = tcam(instanceName)
            % Constructor

            obj = obj@kam(instanceName)
            obj = obj.calculate_spring_force();
            obj = obj.calculate_cutting_force();
            obj = obj.calculate_motor_torque();
        end

        function obj = calculate_displacement(obj)
            % Calculating displacement
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
        end

        function obj = calculate_load_displacement(obj)
            obj.load_displacement.data = obj.displacement;
            % load displacement and follower (roller) displacement of a
            % translating cam are the same
        end

        function obj = calculate_roller_position(obj)
            %

            obj.roller_position = obj.displacement + obj.rPrime; % cam is at center
        end

        function obj = calculate_pressure_angle(obj)
            %
            % tan a = {ds/d(theta)}/(s + rb + rr) %theta in degree

            radianStep = deg2rad(obj.step);
            d_s = [diff(obj.displacement) obj.displacement(1)-obj.displacement(length(obj.displacement))];
            v_theta = d_s/radianStep; % differentiate s with respect to theta in radian

            pitch_radius = obj.displacement + obj.rRoller + obj.rBase;
            tanPressureAngle = v_theta./pitch_radius;

            obj.pressureAngle.data = rad2deg(atan(tanPressureAngle));

            obj.max_pressureAngle = max(obj.pressureAngle.data);
        end


        function animation(obj)
            %

            angleLineFactor = 4;
            [angleEndPointX,angleEndPointY] = offsetIn(obj.pitchX,obj.pitchY,-angleLineFactor*obj.rRoller);
            theta2 = deg2rad(obj.theta);
            h = obj.fullStroke;
            plot(0,0,'o','MarkerFaceColor','r');
            % Pitch Circle
            rotatedPitch = rotateCw([obj.pitchX;obj.pitchY],-pi/2);

            % Cam Profile
            rotatedCam = rotateCw([obj.camSurfX;obj.camSurfY],-pi/2);

            % Roller
            index = linspace(0,2*pi,100);
            xC = obj.rRoller*cos(index);
            rotatedAngleEnd = rotateCw([angleEndPointX(1);angleEndPointY(1)],-pi/2);

            % Animation
            maxDim = obj.rBase + abs(h) + 2*obj.rRoller + 10;

            for i = 1:length(theta2)

                hold on
                grid on;
                grid minor;
                axis equal;

                xlim([-maxDim maxDim]);
                ylim([-maxDim maxDim]);

                j = theta2(i)-pi/2;
                rotatedPitch = rotateCw([obj.pitchX;obj.pitchY],j);
                rotatedCam = rotateCw([obj.camSurfX;obj.camSurfY],j);

                % Pressure angle
                yC = obj.rRoller*sin(index) + obj.roller_position(i);
                rollerCenterY = obj.roller_position(i);
                rollerCenterY_angle = obj.roller_position(i)+angleLineFactor*obj.rRoller;
                angle1y = [rollerCenterY rollerCenterY_angle];

                rotatedAngleEnd = rotateCw([angleEndPointX(i);angleEndPointY(i)],j);
                angle2x = [0 rotatedAngleEnd(1)];
                angle2y = [rollerCenterY rotatedAngleEnd(2)];

                temp6 = strcat('曲率半径　',num2str(obj.curvature.data(i)),' mm     ');
                temp5 = strcat('圧角　',num2str(obj.pressureAngle.data(i)),'^o     ');
                temp2 = strcat('変位　',num2str(obj.roller_position(i)-obj.rPrime),' mm     ');
                temp3 = strcat('回転角度　',num2str(obj.theta(i)),'^o   ');
                updatedTitle = {temp3; temp2; temp5; temp6};
                title(updatedTitle,'Color',[0 0.4470 0.7410],'FontSize',14);

                temp4 = strcat('位置　',num2str(obj.roller_position(i)),' mm     ');
                ylabel(temp4,'Color',obj.angleColor,'FontSize',15);
                temp1 = strcat('経過時間　',num2str(obj.time(i)),' s     ');
                xlabel(temp1,'Color',obj.angleColor,'FontSize',15);

                plot(0,0,'o','MarkerFaceColor','r');
                hold on
                plot(rotatedPitch(1,:),rotatedPitch(2,:),'color',obj.pitchColor);
                hold on
                plot(rotatedCam(1,:),rotatedCam(2,:),'color',obj.camColor);
                hold on
                plot(xC,yC,'color',obj.rollerColor);
                hold on
                plot(0,rollerCenterY,'o','MarkerFaceColor',[0 0.4470 0.7410]); % Roller center
                hold on
                plot([0 0],angle1y,'MarkerFaceColor',[0 0.4470 0.7410]);
                hold on
                plot(angle2x,angle2y,'MarkerFaceColor',[0 0.4470 0.7410]);
                pause(0.005)
                clf;
            end

            hold on
            grid on;
            grid minor;
            axis equal;
            maxDim = obj.rBase + abs(h) + 2*obj.rRoller + 10;
            xlim([-maxDim maxDim]);
            ylim([-maxDim maxDim]);
            plot(0,0,'o','MarkerFaceColor','r');
            hold on
            plot(rotatedPitch(1,:),rotatedPitch(2,:),'color',obj.pitchColor);
            hold on
            plot(rotatedCam(1,:),rotatedCam(2,:),'color',obj.camColor);
            hold on
            plot(xC,yC,'color',obj.rollerColor);
            hold on
            plot(0,rollerCenterY,'o','MarkerFaceColor',[0 0.4470 0.7410]); % Roller center
            hold on
            plot([0 0],angle1y,'MarkerFaceColor',[0 0.4470 0.7410]);
            hold on
            plot(angle2x,angle2y,'MarkerFaceColor',[0 0.4470 0.7410]);
        end

        function obj = calculate_spring_force(obj)
            % Find the force need for follower to always remains in
            % contact with cam.
            % k is spring constant in N/mm
            % initialDisplacement is initial displacement in mm

            m = obj.m_load + obj.m_roller;
            negativeAcceleration = obj.acceleration.data;
            negativeAcceleration(obj.acceleration.data > 0) = 0;
            minimum_spring_force = m*negativeAcceleration - obj.fFriction;
           

            springDisplacement = obj.initialSpringDisplacement + obj.displacement;
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
            m = obj.m_load + obj.m_roller;
            tanPressureAngle = tan(deg2rad(obj.pressureAngle.data));
            parallelForce = m*obj.acceleration.data - obj.fSpring + obj.fCut + obj.fFriction; % in N
            perpendicularForce = parallelForce.*tanPressureAngle;
            obj.motorTorque.data = obj.roller_position/1000.*perpendicularForce;
        end


    end
end