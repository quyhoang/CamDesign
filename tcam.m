classdef tcam < kam
    %Translating cam
    
    properties
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



        function simulation(obj)
            %
            angleColor = 'b';
            angleLineFactor = 4;
            [angleEndPointX,angleEndPointY] = offsetIn(obj.pitchX,obj.pitchY,-angleLineFactor*obj.rRoller);
            theta2 = deg2rad(obj.theta);
            h = obj.fullStroke;
            position = obj.displacement + obj.rBase + obj.rRoller;
            figure;
            hold on
            plot(0,0,'o','MarkerFaceColor','r');
            % Pitch Circle
            rotatedPitch = rotateCw([obj.pitchX;obj.pitchY],-pi/2);
            pl = plot(rotatedPitch(1,:),rotatedPitch(2,:),'color',obj.pitchColor);
            axis equal;
            maxDim = obj.rBase + abs(h) + 2*obj.rRoller + 10;
            xlim([-maxDim maxDim]);
            ylim([-maxDim maxDim]);

            hold on;
            grid on;
            grid minor;
            pl.XDataSource = 'xx';
            pl.YDataSource = 'yy';


            % Cam Profile
            rotatedCam = rotateCw([obj.camSurfX;obj.camSurfY],-pi/2);
            pl2 = plot(rotatedCam(1,:),rotatedCam(2,:),'color',obj.camColor);
            axis equal;
            maxDim = obj.rBase + abs(h) + 2*obj.rRoller + 10;
            xlim([-maxDim maxDim]);
            ylim([-maxDim maxDim]);

            hold on;
            pl2.XDataSource = 'xx2';
            pl2.YDataSource = 'yy2';

            % Roller

            index = linspace(0,2*pi,100);
            xC = obj.rRoller*cos(index);
            yC = obj.rRoller*sin(index) + position(1);
            pl3 = plot(xC,yC,'color',obj.rollerColor);
            pl3.YDataSource = 'yC';

            rollerCenterY = position(1);
            pl4 = plot(0,rollerCenterY,'o','MarkerFaceColor',[0 0.4470 0.7410]); % Roller center
            pl4.YDataSource = 'rollerCenterY';

            maxDim = obj.rBase + abs(h) + 2*obj.rRoller + 10;
            xlim([-maxDim maxDim]);
            ylim([-maxDim maxDim]);
            hold on;

            % Pressure angle
            rollerCenterY_angle = position(1)+angleLineFactor*obj.rRoller;
            angle1y = [rollerCenterY rollerCenterY_angle];
            pl5 = plot([0 0],angle1y,'MarkerFaceColor',[0 0.4470 0.7410]);
            pl5.YDataSource = 'angle1y';
            hold on;

            rotatedAngleEnd = rotateCw([angleEndPointX(1);angleEndPointY(1)],-pi/2);
            angle2x = [0 rotatedAngleEnd(1)];
            angle2y = [rollerCenterY rotatedAngleEnd(2)];
%             pl6 = plot(angle2x,angle2y,'MarkerFaceColor',[0 0.4470 0.7410]);
            pl6.XDataSource = 'angle2x';
            pl6.YDataSource = 'angle2y';

            % Animation
            for i = 1:length(theta2)
                j = theta2(i)-pi/2;

                rotatedPitch = rotateCw([obj.pitchX;obj.pitchY],j);
                xx = rotatedPitch(1,:);
                yy = rotatedPitch(2,:);

                rotatedCam = rotateCw([obj.camSurfX;obj.camSurfY],j);
                xx2 = rotatedCam(1,:);
                yy2 = rotatedCam(2,:);

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
                [titleAni,] = title(updatedTitle,'Color',[0 0.4470 0.7410],'FontSize',14);

                temp4 = strcat('位置　',num2str(position(i)),' mm     ');
                ylabel(temp4,'Color',angleColor,'FontSize',15);
                temp1 = strcat('経過時間　',num2str(obj.time(i)),' s     ');
                xlabel(temp1,'Color',angleColor,'FontSize',15);

                refreshdata
                pause(0.01)
            end
        end
    end
end
