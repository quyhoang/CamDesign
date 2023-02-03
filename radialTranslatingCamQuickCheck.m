% CAM Design Quick Check
% Dwell - Rise - Dwell - Return translating CAM
% 2023-02-03
clc; close all;


%%
%============================================
% INPUT 入力
%============================================
% eventAngle = [rise start - rise end - return start- return end]
eventAngle = [0 30 190 220]; % degree at which the rise/return starts/ends
h = 15; % stroke in mm
RPM = 200; % motor velocity in rounds per minutes
m = 1; % follower mass in kg
rRoller = 6;
rBase = 40;

%%




















% recommended values
maxPressureAngle_deg = 30; % in degree
kFriction = 0.7;
sampleRate = 5; % for showing roller on pitch curve with distance in degree
step = .5; % for caculation, the smaller the more accurate, sampling rate in degree

%============================================
% PRELIMINARY CALCULATION
%============================================

maxPressureAngle = deg2rad(maxPressureAngle_deg);
bRise = eventAngle(2) - eventAngle(1) ; %rise period
bReturn = eventAngle(4) - eventAngle(3) ; %return period

% point in time with acceleration change
% points of events = [1-rise, 2-rise +1/8, 3-rise +7/8, 4-rise end, 5-return, 6-return +1/8, 7-return +7/8, 8-return end]
point = [eventAngle(1) eventAngle(1)+bRise/8 eventAngle(1)+7*bRise/8 eventAngle(2) eventAngle(3) eventAngle(3)+bReturn/8 eventAngle(3)+7*bReturn/8 eventAngle(4)];

%============================================
% CHOOSING BASE RADIUS AND ROLLER RADIUS
%============================================
Cv = 1.7596; % Modified sinusoidal
if h >= 0
    B = deg2rad(bRise);
elseif h < 0
    B = deg2rad(bReturn);
end

list_rRoller = linspace(5,40,200); % roller radius in mm
primeRadius = Cv*abs(h)/B/tan(maxPressureAngle)-abs(h)/2;
list_rBase = primeRadius - list_rRoller;
list_rBase = list_rBase .* (list_rBase>0);

rPrime = rBase + rRoller; %mm - Pitch circle prime radius

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
displacement = [sDwe1 sRise1 sRise2 sRise3 sDwe2 sReturn1 sReturn2 sReturn3 sDwe3] + rPrime;

%============================================
% VELOCITY
%============================================
% velocity with respect to time
velocity = diff(displacement)/timeStep;
velocity = [velocity displacement(1)-displacement(length(displacement))]; %add the last element to make the length of vv and theta equal

tempV = strcat('最大速度 ',num2str(max(velocity)),' mm/s');disp(tempV)
disp('---------')
%============================================
% ACCELERATION
%============================================
% acceleration with respect to time
acceleration = diff(velocity)/timeStep;
acceleration = [acceleration velocity(1)-velocity(length(velocity))]/1000;
tempA = strcat('最大加速  ', num2str(max(acceleration)),' m/s^2'); disp(tempA)
disp('---------')

%%
%============================================
% PRESSURE ANGLE 圧角
%============================================
% tan a = {ds/d(theta)}/(s + rb + rr) %theta in degree

radianStep = deg2rad(step);
d_s = [diff(displacement) displacement(1)-displacement(length(displacement))];
v_theta = d_s/radianStep; % differentiate s with respect to theta in radian

pitch_radius = displacement + rRoller;
tanPressureAngle = v_theta./pitch_radius;

pressureAngle = rad2deg(atan(tanPressureAngle));
disp(strcat('最大圧角 ',num2str(max(pressureAngle)),'度'));
if (max(pressureAngle) > maxPressureAngle_deg)
disp('NG')
else
disp('OK')
end
disp('---------')


%%
%============================================
% POSITION IN POLAR COORDINATES
%============================================

theta2 = deg2rad(theta);

% Converting Polar to Cartesian Coordinate System
[pitchX,pitchY] = pol2cart(theta2,displacement);

%%
%============================================
% CAM PROFILE 
%============================================
%
% Plot cam profile
[camSurfX, camSurfY] = offsetIn(pitchX,pitchY,rRoller);

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
% MOTOR TORQUE 
%============================================
% suppose that the spring force generate an acceleration at least equal to
% that by motor, and friction coefficient is 0.7
% kFriction = 1;
fFriction = kFriction * m * 9.8; 
% Friction force is in opposite direction of velocity
temp1 = -(velocity>0);
temp2 = (velocity<0);
fFrictionDirection = temp1 + temp2;
fFriction = fFriction.*fFrictionDirection;

deltaXoInput = 1;
deltaXo = deltaXoInput/1000; % convert to m
strokeT = (displacement-rPrime)/1000; % m
deltaX = strokeT + deltaXo; % spring displacement 

minFspring = m*acceleration - fFriction; % positive direction is upward, measured in N
minKspring = -minFspring./deltaX;
minAcceptableKspring = max(minKspring); 

maxDeltaX = h+deltaXoInput:1:3*h;
maxFspring = minAcceptableKspring/1000.*maxDeltaX;
strminKspring = num2str(minAcceptableKspring/1000);
Fmax = h+deltaXoInput+1;
Nmax = Fmax*minAcceptableKspring/1000+1;

selectedKspring = Nmax/Fmax*1000;

fSpring = -selectedKspring.*deltaX;

parallelForce = m*acceleration - fSpring - fFriction; % in N
perpendicularForce = parallelForce.*tanPressureAngle;
motorTorque = displacement/1000.*perpendicularForce;
torqueTitle = strcat('最大トルク  ', num2str(max(motorTorque)),' Nm');
disp(torqueTitle);
disp('---------')

%============================================
% FUNCTIONS 
%============================================


function Y = rotateCw(X,theta)
% rotate clockwise
rotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
Y = rotMat * X;
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

function [xOffset,yOffset] = offsetIn(x,y,R)
% Generate 2 vectors holding coordinates of offset curve by R from a closed
% curve. 
% If R is positive, the offset curve will be inside the original curve. 
% If R is negative, the offset curve will be outside the original curve.

L = length(x); % the length of both vectors
% Check whether the input is closed
if (x(1) ~= x(L))
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


function [radius,radArg] = radCurv(f,arg)

% Numerically computing radius of curvature
% arg is used only to define step. In this program arg is array of angles
    step = arg(2)-arg(1);
    f1 = diff(f)/step;
    f2 = diff(f1)/step;
    f1 = regulate(f1,f2);
    ff = regulate(f,f2);
    radArg = regulate(arg,f2);
    
    k1 = 2*sqrVec(f1) + sqrVec(ff) - ff.*f2;
    k2 = (sqrVec(f1) + sqrVec(ff)).^(3/2);
    radius = k2./k1;
end

function sqrVec = sqrVec(x)
    sqrVec = x.*x;
end

function regulated = regulate(x,ref)
    regulated = x(1:length(ref));
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
