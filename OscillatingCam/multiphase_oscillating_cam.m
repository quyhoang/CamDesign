% Multi-phase oscillating CAM
% 2023-07-12

%% INPUT 入力
%============================================
clc; close all; clear;

% All values in degree
% Input vectors
transition = [30 0; 130 3.5; 260 3.5; 350 0];
transition_angle = transition(:, 1)';
transition_displacement = transition(:, 2)';

initial_angular_displacement = 30;  % initial angular displacement in coordinate system
s_rad_initial = deg2rad(initial_angular_displacement); 

l_roller = 65; % distance from arm rotating axis to roller center
l_load = 52; % distance from arm center to load
% l_roller/l_load = displacement_roller/displacement_load
% displacement_roller = displacement_load*l_roller/l_load
% All calculation is performed on displacement of roller
%transition_displacement = transition_displacement*l_roller/l_load;

m_load = 1; % load mass in kilogram
m_roller = 0.1; % roller mass in kilogram
rRoller = 9.5;

m_rocker = 0.4; % rocker arm mass in kilogram

rocker2cam = 94.2; % distance between rocker arm axis and cam axis
RPM = 200; % motor velocity in rounds per minutes
% m = 1; % follower mass in kg

%============================================

% recommended values
maxPressureAngle_deg = 20; % in degree
kFriction = 0.7;
sampleRate = 5; % for showing roller on pitch curve with distance in degree
step = .5; % for caculation, the smaller the more accurate, sampling rate in degree

% graphic color
rollerColor = [0.4660 0.6740 0.1880];
pitchColor = [0.8500 0.3250 0.0980];
camColor = 'b';
%============================================
% PRELIMINARY CALCULATION
%============================================

maxPressureAngle = deg2rad(maxPressureAngle_deg);

theta = 0:step:360;
T = 60/RPM; % period of moving 360 degree, in second
time = linspace(0,T,length(theta));
timeStep = T/size(time,2); % convert step in degree to step in time


%% DISPLACEMENT
%============================================

% Initialize the displacement array
displacement = zeros(size(theta));

% Get the transition points
filtered_pairs = filterConsecutivePairs(transition_angle, transition_displacement);

% Iterate over the transition points
for i = 1:size(filtered_pairs, 1)
    % Extract the start and end points of the transition
    point = filtered_pairs(i, :);
    h = transition_displacement(transition_angle == point(2)) - transition_displacement(transition_angle == point(1));
    bRise = point(2) - point(1);
    
    % Calculate the three sections of the displacement curve
    tempTheta1 = theta(theta >= point(1) & theta < point(1) + bRise/8) - point(1);
    sRise1 = h/(4+pi)*(pi*tempTheta1/bRise - 1/4*sin(4*pi*tempTheta1/bRise));
    
    tempTheta2 = theta(theta >= point(1) + bRise/8 & theta < point(1) + 7*bRise/8) - point(1);
    sRise2 = h/(4+pi)*(2 + pi*tempTheta2/bRise - 9/4*sin(pi/3 + 4*pi/3*tempTheta2/bRise));
    
    tempTheta3 = theta(theta >= point(1) + 7*bRise/8 & theta <= point(2)) - point(1);
    sRise3 = h/(4+pi)*(4 + pi*tempTheta3/bRise - 1/4*sin(4*pi*tempTheta3/bRise));
    
    % Combine all parts and add the previous displacement value
    sRise = [sRise1, sRise2, sRise3];
    
    % Update the displacement array
    displacement(theta >= point(1) & theta <= point(2)) = displacement(theta >= point(1) & theta <= point(2)) + sRise;
    displacement(theta > point(2)) = displacement(theta > point(2)) + sRise(end);
    % Update the cumulative displacement for the next period
end


% Plot position vs angle in cartesian coordinate
figure;
plot(theta, displacement);
xlim([0 360]);
xlabel({'degree'},'FontSize',15,'FontWeight','light','Color','b');
ylim([-2*abs(h)+h/2 2*abs(h)+h/2]);
ylabel({'位置','mm'},'FontSize',15,'FontWeight','light','Color','b');

figure;
plot(theta,displacement);
grid on;
grid minor;
xlim([0 360]);
xlabel({'回転角度'},'FontSize',15,'FontWeight','light','Color','b');
ylim([-2*abs(h)+h/2 2*abs(h)+h/2]);
ylabel({'位置','mm'},'FontSize',15,'FontWeight','light','Color','b');
% legend("位置 s");
[tit,] = title({'';'S V A Diagram'},{['モーター回転速度 ',num2str(RPM),'rpm   ','T = ', num2str(T),'s'];''},...
    'Color','blue');
tit.FontSize = 15;

%% VELOCITY
%============================================
% velocity with respect to time
velocity = diff(displacement)/timeStep;
velocity = [velocity displacement(1)-displacement(length(displacement))]; %add the last element to make the length of vv and theta equal
figure;
plot(theta,velocity);
grid on;
grid minor;
xlim([0 360]);
xlabel({'回転角度'},'FontSize',15,'FontWeight','light','Color','b');
ylabel({'速度','mm/s'},'FontSize',15,'FontWeight','light','Color','b');
%title({'';'速度　vs　時間';''},'Color','b','FontSize',15,'FontWeight','light');
tempV = strcat('最大速度 ',num2str(max(velocity)),' mm/s');
title(tempV ,'Color','b','FontSize',15,'FontWeight','light');

%% ACCELERATION
%============================================
% acceleration with respect to time
acceleration = diff(velocity)/timeStep;
acceleration = [acceleration velocity(1)-velocity(length(velocity))]/1000;
figure;
plot(theta,acceleration);
grid on; 
grid minor;
xlim([0 360]);
xlabel({'回転角度'},'FontSize',15,'FontWeight','light','Color','b');
ylabel({'加速','m/s^2'},'FontSize',15,'FontWeight','light','Color','b');

%title({'';'加速　vs　時間';''},'Color','b','FontSize',15,'FontWeight','light');
tempA = strcat('最大加速  ', num2str(max(acceleration)),' m/s^2');
title(tempA,'Color','b','FontSize',15,'FontWeight','light');

disp(tempV)
disp(tempA)

%% ANGULAR DISPLACEMENT
%============================================

s2rad = displacement/l_load; % convert arc length to angular displacement
s2rad = s2rad + s_rad_initial; % angular displacement
thetaRadian = deg2rad(theta);

roller_position1 = l_roller*exp(s2rad*1i);  %unregulated position, rocker is at center
roller_position = roller_position1 - rocker2cam; % cam is at center, rocker axis is at (-rocker2cam,0)
pitchCurve = roller_position.*exp(thetaRadian*1i); 
% move points counterclockwise, this means the cam rotates clockwise

rollerCenterX = real(roller_position);
rollerCenterY = imag(roller_position);

pitchX = real(pitchCurve); pitchY = imag(pitchCurve);
[camSurfX,camSurfY] = offsetIn(pitchX,pitchY,rRoller);

%% RADIUS OF CURVATURE
%============================================

curvature = zeros(size(camSurfX));
L = length(camSurfX);
% Boundary. Note that the first and the last points on profile curve are the
% same
curvature(1) = circumscribedR([camSurfX(L-1) camSurfX(1) camSurfX(2)],[camSurfY(L-1) camSurfY(1) camSurfY(2)]); 
curvature(L) = circumscribedR([camSurfX(L-1) camSurfX(L) camSurfX(2)],[camSurfY(L-1) camSurfY(L) camSurfY(2)]); 

for k = 1:1:L-2
X = camSurfX(k:1:k+2);
Y = camSurfY(k:1:k+2);
curvature(k+1) = circumscribedR(X,Y);
end

figure;
yyaxis left
angleColor = 'b';
semilogy(theta, curvature,'Color',angleColor);
ax = gca;
ax.YColor = angleColor;
grid on;
grid minor;
xlim([0 360]);
xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color',angleColor);
ylabel({'曲率半径','mm'},'FontSize',15,'FontWeight','light','Color',angleColor);

yyaxis right
strokeColor = [0.6350 0.0780 0.1840];
plot(theta,displacement,'Color',strokeColor);
ax = gca;
ax.YColor = strokeColor;
grid on;
grid minor;
xlim([0 360]);
% ylim([rPrime-2*abs(h)+h/2 rPrime+2*abs(h)+h/2]);
% xlabel({'角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
ylabel({'位置','mm'},'FontSize',15,'FontWeight','light','Color',strokeColor);
hold on

tempP = strcat('最小曲率半径 ',num2str(min(curvature)),'mm');
%  title(temp,'Color','b','FontSize',15,'FontWeight','light');

title({'';'曲率半径・位置　vs　回転角度'; tempP; ''},'Color','b','FontSize',15,'FontWeight','light');

disp(strcat('最小曲率半径: ',num2str(min(curvature)),'mm'));
% input(['\n','何かキーを押すとカムのピッチ円を表示します。','\n'])


%%
%============================================
% PRESSURE ANGLE
%============================================
rockerNormalAngle = rad2deg(angle(roller_position1))+90;

pressureAngle = zeros(size(theta));
normalPhase = zeros(size(theta));
contactPointonCamPhase = zeros(size(theta));
contactPoint2CamDistance = zeros(size(theta));

camCenter = [0 0];

for i = 1:length(theta)
j = thetaRadian(i);

% Update roller center position
tempRollerCenter = [rollerCenterX(i) rollerCenterY(i)];

% Update cam-roller contact point
rotatedCam = rotateCw([camSurfX;camSurfY],j);
contactPoint = [rotatedCam(1,i) rotatedCam(2,i)];

normalPhase(i) = segmentPhase(contactPoint,tempRollerCenter);
contactPointonCamPhase(i) = segmentPhase(camCenter,contactPoint);
contactPoint2CamDistance(i) = norm(contactPoint); %distance from cam center to contact point
end

pressureAngle = (rockerNormalAngle - normalPhase);

figure;
angleColor = 'b';
plot(theta, pressureAngle,'Color',angleColor);
xlim([0 360]);
grid on;
grid minor;

xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color',angleColor);
ylabel({'圧角','degree'},'FontSize',15,'FontWeight','light','Color',angleColor);

tempP = strcat('最大圧角 ',num2str(max(abs(pressureAngle))),'^o');
title({'';'圧角・位置　vs　回転角度'; tempP; ''},'Color','b','FontSize',15,'FontWeight','light');

%%
%============================================
% MOTOR TORQUE 
%============================================

angularAcceleration = acceleration/(l_roller/1000);
inertialMoment1 = m_roller*(l_roller/1000)^2;
inertialMoment2 = 1/3*m_rocker*(l_load/1000)^2;
inertialMoment = inertialMoment1 + inertialMoment2;
combinedTorque = inertialMoment*angularAcceleration;
regulatedMinimumSpringTorque = min(combinedTorque)*ones(size(combinedTorque));

estimatedLoadTorque = -l_load/1000*(m_load*10);

motorInducedTorque = combinedTorque - regulatedMinimumSpringTorque - estimatedLoadTorque;

% Calculate motor force
rockerAngle = rad2deg(angle(roller_position1));
rocker2forceAngle = deg2rad(normalPhase - rockerAngle);
forceModuli = 1/l_roller*motorInducedTorque./sin(rocker2forceAngle);
camForceArm2ForceAngle = deg2rad(contactPointonCamPhase-normalPhase);
motorTorque = contactPoint2CamDistance.*forceModuli.*sin(camForceArm2ForceAngle);
figure; 

plot(theta,motorTorque);
xlim([0 360]);
xlabel('角度') ; % misumi nomenclature
ylabel('トルク (Nm)') ;
torqueTitle = strcat('最大トルク  ', num2str(max(motorTorque)),' Nm');
title(torqueTitle,'Color','b','FontSize',15,'FontWeight','light');
grid on; grid minor;

% save torque for later use. (e.g. Combine with other torques to determine total
% motor torque.)
inputMessage =  (['トルクを保存します。ファイル名とトルク名を入力してください。', ...
    '\nファイル名とトルク名は同じになります。']);
fileName = input(inputMessage,"s");

if isempty(fileName)
   disp('ファイル名の入力がないので、トルクを保存しません。');
else
   torqueName = fileName;   
   eval([torqueName, '=motorTorque;']); % we cannot directly assign value to the new variable name
   % this is dynamic field referencing syntax to assign the value of the original variable to the new variable name.
   save(fileName, torqueName); %save the variable into the .mat file of the same name.
   disp('トルクを保存しました。');
   disp('トルク価値を使う前に、次の構文使ってください：');
   disp(['load ',fileName]);
end




%%
%============================================
% ANIMATION 
%============================================


figure;
% Draw rocker arm
armCenterX = -rocker2cam;
armCenterY = 0;
% rollerCenterX = real(roller_position);
% rollerCenterY = imag(roller_position);

armX = [armCenterX,rollerCenterX(1)];
armY = [armCenterY,rollerCenterY(1)];
pl2 = plot(armX,armY); hold on
pl2.XDataSource = 'armX';
pl2.YDataSource = 'armY';

plot(armCenterX,armCenterY,'o','MarkerFaceColor','r');

% Draw roller
index = linspace(0,2*pi,100);
xC = rRoller*cos(index) + rollerCenterX(1);
yC = rRoller*sin(index) + rollerCenterY(1);
pl3 = plot(xC,yC,'color',rollerColor); hold on
pl3.XDataSource = 'xC';
pl3.YDataSource = 'yC';

% Draw roller center
tempRollerCenterX = rollerCenterX(1);
tempRollerCenterY = rollerCenterY(1);
pl4 = plot(tempRollerCenterX,tempRollerCenterY,'o','MarkerFaceColor',[0 0.4470 0.7410]); hold on
pl4.XDataSource = 'tempRollerCenterX';
pl4.YDataSource = 'tempRollerCenterY';


% Draw pitch curve
% pitchX = real(pitchCurve); pitchY = imag(pitchCurve);
xx = pitchX;
yy = pitchY;
pl = plot(xx,yy,'color',pitchColor);hold on; 
pl.XDataSource = 'xx';
pl.YDataSource = 'yy';

% Draw cam profile
% [camSurfX,camSurfY] = offsetIn(pitchX,pitchY,rRoller);
xx2 = camSurfX;
yy2 = camSurfY;
pl2 = plot(xx2,yy2,'color',camColor); hold on;
pl2.XDataSource = 'xx2';
pl2.YDataSource = 'yy2';

plot(0,0,'o','MarkerFaceColor','b');

% Draw cam-roller contact point
contactPointX = camSurfX(1);
contactPointY = camSurfY(1);
pl5 = plot(contactPointX, contactPointY,'.','color','r'); hold on;
pl5.XDataSource = 'contactPointX';
pl5.YDataSource = 'contactPointY';
axis equal; grid on; grid minor;

initialContactPoint = [contactPointX contactPointY];
writematrix(initialContactPoint,'initialContactPoint.txt','Delimiter','tab');


maxDim = rocker2cam+5;
xlim([-maxDim maxDim]);
ylim([-maxDim maxDim]);

%============================================
% UPDATE FRAME 
%============================================

% Cam motion simulation
prompt = "Show cam motion? カムの動きのシミュレーションを表示しますか? y/n [n]: ";
txt = input(prompt,"s");

if (txt == 'y')

for loopNumber = 1:1

for i = 1:length(theta)
j = thetaRadian(i);

% Update rocker arm
armX = [armCenterX,rollerCenterX(i)];
armY = [armCenterY,rollerCenterY(i)];

% Update roller position
xC = rRoller*cos(index) + rollerCenterX(i);
yC = rRoller*sin(index) + rollerCenterY(i);

% Update roller center position
tempRollerCenterX = rollerCenterX(i);
tempRollerCenterY = rollerCenterY(i);

% Update rotated pitch
rotatedPitch = rotateCw([pitchX;pitchY],j);
xx = rotatedPitch(1,:);
yy = rotatedPitch(2,:);

% Update rotated cam
rotatedCam = rotateCw([camSurfX;camSurfY],j);
xx2 = rotatedCam(1,:);
yy2 = rotatedCam(2,:);

% Update cam-roller contact point
% rotatedCam = rotateCw([camSurfX;camSurfY],j);
contactPointX = rotatedCam(1,i);
contactPointY = rotatedCam(2,i);



refreshdata
pause(0.01)
end 

end 

end

% Export data for using in 3D CAD package
prompt1 = (['CreoAutomation.exe がアクティブで、' ...
    '\ntxt データが Creo 作業ディレクトリに保存されている場合、' ...
    '\n「gcam」を押して Enter を押すと、' ...
    '\nカムの 3D モデルが作成されます。']);
prompt = prompt1 + "\n\nExport data (.txt)? データ (.txt) をエクスポートしますか? y/n [n] : ";
txt = input(prompt,"s");

% Export if user inputs 'y'
% otherwise (input 'n' or just Enter), skip
if (txt == 'y')
    x_cord = transpose(camSurfX);
    y_cord = transpose(camSurfY);
    z_cord = zeros(size(x_cord));
    
    camProfile = [x_cord y_cord z_cord];
%     writematrix(camProfile,'camProfile.xlsx');
    writematrix(camProfile,'camProfile.txt','Delimiter','tab');

    % take element to show cam direction
    indices = 2.^(0:log2(length(x_cord)));
    x_direction = x_cord(indices);
    y_direction = y_cord(indices);
    z_direction = ones(size(x_direction));
    camDirection = [x_direction y_direction z_direction];
    writematrix(camDirection,'camDirection.txt','Delimiter','tab');
end


%% 
%============================================
% Save all data 
%============================================
camName = input("カムの名を入力してください。","s");
if ~isempty(camName)
    save(camName)
    disp(['データは全て', camName,'.mat に保存されています。']);
    disp('このデータを使う際に、次の構文使ってください：');
   disp([camName, ' = load(''''', camName, '.mat'''');']);
   disp('それから、例えば変位データを使う際に、次の構文使ってください：');
   disp([camName, '.displacement']);
else
    disp('入力がないので、データは保存されていません。');
end

%% 
%============================================
% FUNCTIONS 
%============================================

function [xOffset,yOffset] = offsetIn(x,y,R)
% Generate 2 vectors holding coordinates of offset curve by R from a closed
% curve. 
% If R is positive, the offset curve will be inside the original curve. 
% If R is negative, the offset curve will be outside the original curve.

L = length(x); % the length of both vectors
% Check whether the input is closed
%if (x(1) ~= x(L)) %the same point could be shifted because of calculation
%precision. Although the difference is extremely small, this condition may
%not hold true, hence the condition in the next line.
if (abs(x(1) - x(L)) > 10^-10)
    disp("The input is not a closed curve. The function will terminate.")
    return
end

xOffset = zeros(size(x));
yOffset = zeros(size(x));

% Boundary. Note that the first and the last points on pitch curve are the
% same
[xOffset(1),yOffset(1)] = normalIn([x(L-1) x(1) x(2)],[y(length(y)-1) y(1) y(2)],R); 
[xOffset(L),yOffset(L)] = normalIn([x(L-1) x(L) x(2)],[y(L-1) y(L) y(2)],R);

for k = 1:1:L-2
X = x(k:1:k+2);
Y = y(k:1:k+2);
[xOffset(k+1),yOffset(k+1)] = normalIn(X,Y,R);
end

end

function [xo,yo] = normalIn(x,y,R)
% x and y are row vector of length 3 (longer vectors don't cause problem,
% but only the first 3 elements will be used. 

% Call the three point represented by x and y A, B and C
% This function return the coordinates of point D such that
% * DB is perpendicular to AC
% * DB has length R
% * D is on the right hand side when moving on the curve ABC from A to C
% if R is positive and left hand side if R is negative

% calculate normal vector <a,b>
a = y(1)-y(3);
b = x(3)-x(1);
k = R/sqrt(a^2+b^2); 

% temporary factor
xo = k*a + x(2);
yo = k*b + y(2);
end

function Y = rotateCw(X,theta)
% rotate clockwise
rotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
Y = rotMat * X;
end

function segmentPhase = segmentPhase(P1,P2)
% return the angle of the vector from point P1 to point P2
complexNum = P2(1) - P1(1) + (P2(2) - P1(2))*1i;
segmentPhase = rad2deg(angle(complexNum));
end

function R = circumscribedR(X,Y)
% Calculate radius of circumscribed circle through three point in vectors X
% and Y of size 1*3

A = [X(1) Y(1)]; B = [X(2) Y(2)]; C = [X(3) Y(3)];

a = norm(B-C);
b = norm(C-A);
c = norm(A-B);

s = (a+b+c)/2;
temp = sqrt(s*(s-a)*(s-b)*(s-c));

R = a*b*c/4/temp;
end


function filtered_pairs = filterConsecutivePairs(transition_angle, transition_displacement)
    % Add 0 to the beginning and end of transition_displacement
    transition_displacement = horzcat(0, transition_displacement, 0);

    % Add 0 and 360 to the beginning and the end of transition_angle
    transition_angle = horzcat(0, transition_angle, 360);

    % Check that all elements in transition_angle are unique
    assert(all(diff(transition_angle) > 0) && (length(unique(transition_angle)) == length(transition_angle)), 'The angles must be unique and monotonic.');

    % Initialize filtered_pairs as an empty array with the maximum possible size
    filtered_pairs = zeros(length(transition_angle)-1, 2); % avoid dynamic allocation of memory

    % Initialize counter
    count = 0;

    % Iterate through the angles
    for i = 1:length(transition_angle)-1
        % Check if the current displacement is different from the next one
        if transition_displacement(i) ~= transition_displacement(i+1)
            % Increment counter
            count = count + 1;
            % Add the pair of consecutive angles to filtered_pairs
            filtered_pairs(count, :) = [transition_angle(i), transition_angle(i+1)];
        end
    end

    % Trim filtered_pairs to its actual size
    filtered_pairs = filtered_pairs(1:count, :);
end
