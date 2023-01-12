% CAM Design Assistant
% Dwell - Rise - Dwell - Return oscillating CAM
% 2022-12-07

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
l_roller = 40; % distance from arm rotating axis to roller center
l_load = 80; % distance from arm center to load

rocker2cam = 80; % distance between rocker arm and cam axes
RPM = 200; % motor velocity in rounds per minutes
m = 1; % follower mass in kg

% recommended values
maxPressureAngle_deg = 20; % in degree
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
camProfile = roller_position.*exp(thetaRadian*1i); 
% move points counterclockwise, this means the cam rotates clockwise 
figure
plot(camProfile)
hold on
plot(roller_position,'o')
axis equal; grid on;



