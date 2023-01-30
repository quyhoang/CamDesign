% CAM Design Assistant
% Dwell - Rise - Dwell - Return oscillating CAM
% 2023-01-30

%%
clc; close all; 
% clear;
%============================================
% INPUT 入力
%============================================
% All values in degree

%TEMPORARY - DEBUGGING only
rPrime = 0;

% Job-specific values
% eventAngle = [rise start - rise end - return start- return end]
eventAngle = [80 120 190 230]; % degree at which the rise/return starts/ends
h = 10; % stroke in mm
l_roller = 60; % distance from arm rotating axis to roller center
l_load = 80; % distance from arm center to load

rRoller = 8;

rocker2cam = 80; % distance between rocker arm and cam axes
RPM = 200; % motor velocity in rounds per minutes
m = 1; % follower mass in kg

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
bRise = eventAngle(2) - eventAngle(1) ; %rise period
bReturn = eventAngle(4) - eventAngle(3) ; %return period

% point in time with acceleration change
% points of events = [1-rise, 2-rise +1/8, 3-rise +7/8, 4-rise end, 5-return, 6-return +1/8, 7-return +7/8, 8-return end]
point = [eventAngle(1) eventAngle(1)+bRise/8 eventAngle(1)+7*bRise/8 eventAngle(2) eventAngle(3) eventAngle(3)+bReturn/8 eventAngle(3)+7*bReturn/8 eventAngle(4)];

theta = 0:step:360;
T = 60/RPM; % period of moving 360 degree, in second
time = linspace(0,T,length(theta));
timeStep = T/size(time,2); % convert step in degree to step in time


%% 
%============================================
% DISPLACEMENT
%============================================

% Rise
temp = theta(theta<point(1));
sDwe1 = zeros(size(temp));
tempTheta = theta(theta >= point(1) & theta < point(2))-point(1);
sRise1 = h/(4+pi)*(pi*tempTheta/bRise - 1/4*sin(4*pi*tempTheta/bRise));
tempTheta = theta(theta >= point(2) & theta < point(3))-point(1);
sRise2 = h/(4+pi)*(2+pi*tempTheta/bRise-9/4*sin(pi/3+4*pi/3*tempTheta/bRise));
tempTheta = theta(theta >= point(3) & theta <= point(4))-point(1);
sRise3 = h/(4+pi)*(4+ pi*tempTheta/bRise - 1/4*sin(4*pi*tempTheta/bRise));

% Dwell
temp = theta(theta > point(4) & theta < point(5));
sDwe2 = zeros(size(temp)) + h;

% Return
tempTheta = theta(theta >= point(5) & theta < point(6))-point(5);
sReturn1 = h/(4+pi)*(4 + pi - pi*tempTheta/bReturn + 1/4*sin(4*pi*tempTheta/bReturn));
tempTheta = theta(theta >= point(6) & theta < point(7))-point(5);
sReturn2 = h/(4+pi)*(2+ pi - pi*tempTheta/bReturn  + 9/4*sin(pi/3+4*pi/3*tempTheta/bReturn));
tempTheta = theta(theta >= point(7) & theta <= point(8))-point(5);
sReturn3 = h/(4+pi)*(pi - pi*tempTheta/bReturn + 1/4*sin(4*pi*tempTheta/bReturn));

% Dwell
temp = theta(theta > point(8) & theta <= 360);
sDwe3 = zeros(size(temp));

% Entire trajectory
s = [sDwe1 sRise1 sRise2 sRise3 sDwe2 sReturn1 sReturn2 sReturn3 sDwe3] + rPrime;

% Plot position vs angle in cartesian coordinate
figure;
subplot(3,1,1);
plot(time,s);
grid on;
grid minor;
xlim([0 T]);
xlabel({'t(s)'},'FontSize',15,'FontWeight','light','Color','b');
ylim([rPrime-2*abs(h)+h/2 rPrime+2*abs(h)+h/2]);
ylabel({'位置','mm'},'FontSize',15,'FontWeight','light','Color','b');
% legend("位置 s");
[tit,] = title({'';'S V A Diagram'},{['モーター回転速度 ',num2str(RPM),'rpm   ','T = ', num2str(T),'s'];''},...
    'Color','blue');
tit.FontSize = 15;

%============================================
% VELOCITY
%============================================
% velocity with respect to time
vv = diff(s)/timeStep;
vv = [vv s(1)-s(length(s))]; %add the last element to make the length of vv and theta equal
subplot(3,1,2);
plot(time,vv);
grid on;
grid minor;
xlim([0 T]);
xlabel({'t(s)'},'FontSize',15,'FontWeight','light','Color','b');
ylabel({'速度','mm/s'},'FontSize',15,'FontWeight','light','Color','b');
%title({'';'速度　vs　時間';''},'Color','b','FontSize',15,'FontWeight','light');
tempV = strcat('最大速度 ',num2str(max(vv)),' mm/s');
title(tempV ,'Color','b','FontSize',15,'FontWeight','light');

%============================================
% ACCELERATION
%============================================
% acceleration with respect to time
aa = diff(vv)/timeStep;
aa = [aa vv(1)-vv(length(vv))]/1000;
subplot(3,1,3);
plot(time,aa);
grid on; grid minor;

xlim([0 T]);
xlabel({'t(s)'},'FontSize',15,'FontWeight','light','Color','b');
ylabel({'加速','m/s^2'},'FontSize',15,'FontWeight','light','Color','b');

%title({'';'加速　vs　時間';''},'Color','b','FontSize',15,'FontWeight','light');
tempA = strcat('最大加速  ', num2str(max(aa)),' m/s^2');
title(tempA,'Color','b','FontSize',15,'FontWeight','light');

disp(tempV)
disp(tempA)

%%
%============================================
% ANGULAR DISPLACEMENT
%============================================

s2rad = s/l_load; % convert arc length to angular displacement
s_rad_initial = pi/6; % initial angular displacement in coordinate system
s2rad = s2rad + s_rad_initial;
thetaRadian = deg2rad(theta);

roller_position1 = l_roller*exp(s2rad*1i);  %unregulated position, rocker is at center
roller_position = roller_position1 - rocker2cam; % cam is at center, rocker axis is at (-rocker2cam,0)
pitchCurve = roller_position.*exp(thetaRadian*1i); 
% move points counterclockwise, this means the cam rotates clockwise 

%%
%============================================
% ANIMATION 
%============================================

figure

% Draw rocker arm
armCenterX = -rocker2cam;
armCenterY = 0;
rollerCenterX = real(roller_position);
rollerCenterY = imag(roller_position);

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
pitchX = real(pitchCurve); pitchY = imag(pitchCurve);
xx = pitchX;
yy = pitchY;
pl = plot(xx,yy,'color',pitchColor);hold on; 
pl.XDataSource = 'xx';
pl.YDataSource = 'yy';

% Draw cam profile
[camSurfX,camSurfY] = offsetIn(pitchX,pitchY,rRoller);
xx2 = camSurfX;
yy2 = camSurfY;
pl2 = plot(xx2,yy2,'color',camColor); hold on;
pl2.XDataSource = 'xx2';
pl2.YDataSource = 'yy2';

plot(0,0,'o','MarkerFaceColor','b');
axis equal; grid on; grid minor;




%============================================
% UPDATE FRAME 
%============================================

for loopNumber = 1:5

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

maxDim = rocker2cam+5;
xlim([-maxDim maxDim]);
ylim([-maxDim maxDim]);

refreshdata
pause(0.01)
end 

end 

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
%if (x(1) ~= x(L))
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

