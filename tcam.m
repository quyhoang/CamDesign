classdef tcam < kam
    %Translating cam

    properties
        rBase % Base radius
    end

    methods
        function obj = tcam(instanceName)
            % Constructor

            obj = obj@kam(instanceName)
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

            obj.pressureAngle = rad2deg(atan(tanPressureAngle));

            obj.max_pressureAngle = max(obj.pressureAngle);
        end


        function animation(obj)
            %

            angleLineFactor = 4;
            [angleEndPointX,angleEndPointY] = offsetIn(obj.pitchX,obj.pitchY,-angleLineFactor*obj.rRoller);
            theta2 = deg2rad(obj.theta);
            h = obj.fullStroke;
            position = obj.displacement + obj.rBase + obj.rRoller;
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
                yC = obj.rRoller*sin(index) + position(i);
                rollerCenterY = position(i);
                rollerCenterY_angle = position(i)+angleLineFactor*obj.rRoller;
                angle1y = [rollerCenterY rollerCenterY_angle];

                rotatedAngleEnd = rotateCw([angleEndPointX(i);angleEndPointY(i)],j);
                angle2x = [0 rotatedAngleEnd(1)];
                angle2y = [rollerCenterY rotatedAngleEnd(2)];

                temp6 = strcat('曲率半径　',num2str(obj.curvature(i)),' mm     ');
                temp5 = strcat('圧角　',num2str(obj.pressureAngle(i)),'^o     ');
                temp2 = strcat('変位　',num2str(position(i)-obj.rPrime),' mm     ');
                temp3 = strcat('回転角度　',num2str(obj.theta(i)),'^o   ');
                updatedTitle = {temp3; temp2; temp5; temp6};
                title(updatedTitle,'Color',[0 0.4470 0.7410],'FontSize',14);

                temp4 = strcat('位置　',num2str(position(i)),' mm     ');
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
    end
end
