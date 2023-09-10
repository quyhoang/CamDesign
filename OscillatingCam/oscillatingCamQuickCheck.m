% CAM Design Assistant
% Dwell - Rise - Dwell - Return oscillating CAM
% 2023-02-01

%%
clc; close all; clear;


% clear;
%============================================
% INPUT 入力
%============================================
% All values in degree

% Job-specific values
% eventAngle = [rise start - rise end - return start- return end]
eventAngle = [0 30 190 220]; % degree at which the rise/return starts/ends
h = 3; % stroke in mm
l_roller = 75; % distance from arm rotating axis to roller center
estLoad = 5; % estimated load in Newton
l_load = 80; % distance from arm center to load
m_roller = 0.1; % roller mass in kilogram
rRoller = 8;

m_rocker = 0.4; % rocker arm mass in kilogram

rocker2cam = 80; % distance between rocker arm and cam axes
s_initial = 40; % initial angle between rocker and cam center - rocker axis in degree
RPM = 200; % motor velocity in rounds per minutes
% m = 1; % follower mass in kg







% recommended values
maxPressureAngle_deg = 30; % in degree
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
% s = [sDwe1 sRise1 sRise2 sRise3 sDwe2 sReturn1 sReturn2 sReturn3 sDwe3] + rPrime;
displacement = [sDwe1 sRise1 sRise2 sRise3 sDwe2 sReturn1 sReturn2 sReturn3 sDwe3];



%============================================
% VELOCITY
%============================================
% velocity with respect to time
velocity = diff(displacement)/timeStep;
velocity = [velocity displacement(1)-displacement(length(displacement))]; %add the last element to make the length of vv and theta equal
tempV = strcat('最大速度 ',num2str(max(velocity)),' mm/s');

%============================================
% ACCELERATION
%============================================
% acceleration with respect to time
acceleration = diff(velocity)/timeStep;
acceleration = [acceleration velocity(1)-velocity(length(velocity))]/1000;

tempA = strcat('最大加速  ', num2str(max(acceleration)),' m/s^2');


disp(tempV)
disp(tempA)

%%
%============================================
% ANGULAR DISPLACEMENT
%============================================

s2rad = displacement/l_load; % convert arc length to angular displacement
s_rad_initial = deg2rad(s_initial); % initial angular displacement in coordinate system
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

%%
%============================================
% RADIUS OF CURVATURE
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

disp(strcat('最小曲率半径 ',num2str(min(curvature)),'mm'));
if (min(curvature) < rRoller)
disp('NG')
else
disp('OK')
end
disp('---------')

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

disp(strcat('最大圧角 ',num2str(max(pressureAngle)),'度'));
if (max(pressureAngle) > maxPressureAngle_deg)
disp('NG')
else
disp('OK')
end
disp('---------')



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

estimatedLoadTorque = -l_load/1000*estLoad;

motorInducedTorque = combinedTorque - regulatedMinimumSpringTorque - estimatedLoadTorque;

% Calculate motor force
rockerAngle = rad2deg(angle(roller_position1));
rocker2forceAngle = deg2rad(normalPhase - rockerAngle);
forceModuli = 1/l_roller*motorInducedTorque./sin(rocker2forceAngle);
camForceArm2ForceAngle = deg2rad(contactPointonCamPhase-normalPhase);
motorTorque = contactPoint2CamDistance.*forceModuli.*sin(camForceArm2ForceAngle);

torqueTitle = strcat('最大トルク  ', num2str(max(motorTorque)),' Nm');
disp(torqueTitle);
disp('---------')




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
