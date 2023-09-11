classdef ocam < kam
    %Oscillating cam

    properties
        m_rocker
        rocker2cam

        l_roller
        l_load

        initial_angular_displacement
        angular_displacement = 0;

        s2rad = 0;
        roller_position_origin = 0 %unregulated position, rocker is at center
    end

    methods
        function obj = ocam(instanceName)
            % Constructor
            obj = obj@kam(instanceName)


        end

        function obj = calculate_roller_position(obj)
            %

            obj.s2rad = obj.displacement/obj.l_load; % convert arc length to angular displacement
            s_rad_initial = deg2rad(obj.initial_angular_displacement);
            obj.s2rad = obj.s2rad + s_rad_initial; % angular displacement

            obj.roller_position_origin = obj.l_roller*exp(obj.s2rad*1i);  %unregulated position, rocker is at center
            obj.roller_position = obj.roller_position_origin - obj.rocker2cam; % cam is at center, rocker axis is at (-rocker2cam,0)
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

            obj.pressureAngle = (rockerNormalAngle - normalPhase);
            obj.max_pressureAngle = max(obj.pressureAngle);
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

                temp6 = strcat('曲率半径　',num2str(obj.curvature(i)),' mm     ');
                temp5 = strcat('圧角　',num2str(obj.pressureAngle(i)),'^o     ');

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
    end
end
