% CAM Design Assistant for Radial Translating Cam
% Dwell - Rise - Dwell - Return 
% Last edited 2023-02-07

%%
clc; close all; clear;
% clear;
%============================================
% Load Input
%============================================

% % Data from input file
% % eventAngle = [rise start - rise end - return start- return end]
% eventAngle = [80 120 190 230]; % degree at which the rise/return starts/ends
% h = 15.2; % stroke in mm
% RPM = 200; % motor velocity in rounds per minutes
% m = 1; % follower mass in kg
% rRoller = 8;

inputFileName = input("設定データ名を入力してください。","s");
load(inputFileName);

% recommended values
maxPressureAngle_deg = 25; % in degree
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

% figure('Name','カムベース円の半径　ｖｓ　ローラー半径');
figure;
plot(list_rRoller,list_rBase);

tempMin = min(list_rRoller);
tempMax = max(list_rRoller);
xScaleDiv = tempMin:1:tempMax;
xticks(xScaleDiv);

tempMin = min(list_rBase);
tempMax = max(list_rBase);
yScaleDiv = round(min(list_rBase)):round((max(list_rBase)-min(list_rBase))/20):max(list_rBase);
yticks(yScaleDiv);

grid on
grid minor

xlabel({'ローラー半径(mm)'},'FontWeight','light','Color','b');
ylabel({'カムベース円の半径(mm)'},'FontWeight','light','Color','b');

temp = strcat('最大圧角  ', num2str(maxPressureAngle_deg),'^o');
temp = {'ローラー半径 vs カムベース円の半径';temp};
title(temp,'Color','b','FontSize',15,'FontWeight','light');


% Ask for user input (user can use the previous figure as reference)
rRollerInputMessage =  (['Please refer Figure 1 and input roller radius: ', ...
    '\n図1を参考して, ローラー半径(mm) を入力してください: \n']);
rBaseInputMessage = (['Please refer Figure 1 and input cam base radius: ', ...
    '\n図1を参考して, カムベース円の半径(mm) を入力してください: \n']);

rRoller = input(rRollerInputMessage);
rBase = input(rBaseInputMessage);

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

% Plot position vs angle in cartesian coordinate
figure;
subplot(3,1,1);
plot(time,displacement);
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
velocity = diff(displacement)/timeStep;
velocity = [velocity displacement(1)-displacement(length(displacement))]; %add the last element to make the length of vv and theta equal
subplot(3,1,2);
plot(time,velocity);
grid on;
grid minor;
xlim([0 T]);
xlabel({'t(s)'},'FontSize',15,'FontWeight','light','Color','b');
ylabel({'速度','mm/s'},'FontSize',15,'FontWeight','light','Color','b');
%title({'';'速度　vs　時間';''},'Color','b','FontSize',15,'FontWeight','light');
tempV = strcat('最大速度 ',num2str(max(velocity)),' mm/s');
title(tempV ,'Color','b','FontSize',15,'FontWeight','light');

%============================================
% ACCELERATION
%============================================
% acceleration with respect to time
acceleration = diff(velocity)/timeStep;
acceleration = [acceleration velocity(1)-velocity(length(velocity))]/1000;
subplot(3,1,3);
plot(time,acceleration);
grid on; grid minor;

xlim([0 T]);
xlabel({'t(s)'},'FontSize',15,'FontWeight','light','Color','b');
ylabel({'加速','m/s^2'},'FontSize',15,'FontWeight','light','Color','b');

%title({'';'加速　vs　時間';''},'Color','b','FontSize',15,'FontWeight','light');
tempA = strcat('最大加速  ', num2str(max(acceleration)),' m/s^2');
title(tempA,'Color','b','FontSize',15,'FontWeight','light');

disp(tempV)
disp(tempA)

input(['\n','何かキーを押すと圧角を表示します。\n'])
% Wait for user response to proceed
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

% figure('Name','圧角・位置 vs 回転角度');
figure;
yyaxis left
angleColor = 'b';
plot(theta, pressureAngle,'Color',angleColor);
ax = gca;
ax.YColor = angleColor;
grid on;
grid minor;
xlim([0 360]);
xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color',angleColor);
ylabel({'圧角','degree'},'FontSize',15,'FontWeight','light','Color',angleColor);

yyaxis right
strokeColor = [0.6350 0.0780 0.1840];
plot(theta,displacement,'Color',strokeColor);
ax = gca;
ax.YColor = strokeColor;
grid on;
grid minor;
xlim([0 360]);
ylim([rPrime-2*abs(h)+h/2 rPrime+2*abs(h)+h/2]);
% xlabel({'角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
ylabel({'位置','mm'},'FontSize',15,'FontWeight','light','Color',strokeColor);
hold on

tempP = strcat('最大圧角 ',num2str(max(pressureAngle)),'^o');
%  title(temp,'Color','b','FontSize',15,'FontWeight','light');

title({'';'圧角・位置　vs　回転角度'; tempP; ''},'Color','b','FontSize',15,'FontWeight','light');

disp(strcat('最大圧角: ',num2str(max(pressureAngle)),'度'));

input(['\n','何かキーを押すとカムのピッチ円を表示します。','\n'])
% Wait for user respond to proceed
%%
%============================================
% POSITION IN POLAR COORDINATES
%============================================

pitchColor = [0.8500 0.3250 0.0980];
camColor = 'b'; % cam profile color
% Plot position in polar coordinate
% figure('Name','極座標のカム ピッチ カーブ');
figure;
theta2 = deg2rad(theta);
polarplot(theta2,displacement,'Color',pitchColor);
grid on;
[title1,] = title({'';'位置　vs　角度';''},'Color','b','FontSize',15,'FontWeight','light');

% Converting Polar to Cartesian Coordinate System
[pitchX,pitchY] = pol2cart(theta2,displacement);

input(['何かキーを押すとカムの輪郭を表示します。','\n'])
% Wait for user response to proceed
%%
%============================================
% CAM PROFILE 
%============================================
% Cam Machining process

% Draw roller around cam curve and on pitch curve, sample rate is defined in input region
% figure('Name','カムの輪郭加工');
sampleRate = round(sampleRate*length(theta2)/360);
x_sample = transpose(pitchX(1:sampleRate:length(pitchX)));
y_sample = transpose(pitchY(1:sampleRate:length(pitchY)));

centers = [x_sample y_sample];
radii = rRoller*ones(length(y_sample),1);

plot(pitchX,pitchY,'Color',pitchColor);
hold on;
rollerColor = [0.4660 0.6740 0.1880];
p = viscircles([x_sample(1),y_sample(1)],rRoller,'LineWidth',1,'Color',rollerColor);
hold off; axis equal; grid on;

% Machining process
prompt = "Show machining process? 加工工程を表示しますか? y/n [n]: ";
txt = input(prompt,"s");
% Animate machining process if user input 'y'
% otherwise (input 'n' or just Enter), skip
if (txt == 'y')
    for k = 2:1:length(x_sample)
        p = viscircles([x_sample(k),y_sample(k)],rRoller,'LineWidth',1,'Color',rollerColor);
        drawnow
    end
end

% figure; % figure('Name','カムの輪郭');
viscircles(centers,radii,'LineWidth',1,'color',rollerColor);
axis equal; hold on; grid on; grid minor;
plot(pitchX,pitchY,'color','r');

% Plot cam profile
[camSurfX, camSurfY] = offsetIn(pitchX,pitchY,rRoller);
plot(camSurfX,camSurfY,'color','b')

% Export data for using in 3D CAD package
prompt1 = (['CreoAutomation.exe がアクティブで、' ...
    '\ntxt データが Creo 作業ディレクトリに保存されている場合、' ...
    '\n「gcam」を押して Enter を押すと、' ...
    '\nカムの 3D モデルが作成されます。']);
prompt = prompt1 + "\n\nExport data (.txt)? データ (.txt) をエクスポートしますか? y/n [n] : ";
txt = input(prompt,"s");

% Animate machining process if user input 'y'
% otherwise (input 'n' or just Enter), skip
if (txt == 'y')
    x_cord = transpose(camSurfX);
    y_cord = transpose(camSurfY);
    z_cord = zeros(length(theta2),1);
    
    camProfile = [x_cord y_cord z_cord];
%     writematrix(camProfile,'camProfile.xlsx');
    writematrix(camProfile,'camProfile.txt','Delimiter','tab');
end



input(['\n','何かキーを押すと曲率半径を表示します。','\n'])
% Wait for user respond to proceed
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
% p = {[(rb + s)^2 + (ds/dtheta)^2]^(3/2)}/[(rb + s)^2 + 2*(ds/dtheta)^2 - (rb+s)*d^2s/dtheta^2]


input(['何かキーを押すとモータートルクを表示します。','\n'])


% Wait for user response to proceed
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

deltaXoInput = input("最初ばね変位 (mm): ");
deltaXo = deltaXoInput/1000; % convert to m
% deltaXo = 0.001; %m
strokeT = (displacement-rPrime)/1000; % m
deltaX = strokeT + deltaXo; % spring displacement 

minFspring = m*acceleration - fFriction; % positive direction is upward, measured in N
minKspring = -minFspring./deltaX;
minAcceptableKspring = max(minKspring); 

figure;
maxDeltaX = h+deltaXoInput:1:3*h;
maxFspring = minAcceptableKspring/1000.*maxDeltaX;
plot(maxDeltaX,maxFspring);
xlabel('Fmax (mm)') ; % misumi nomenclature
xlim([h+deltaXoInput 3*h]);
ylabel('Nmax (N)') ;
grid on; grid minor;

strminKspring = num2str(minAcceptableKspring/1000);
disp(strcat('最小　K＝ ',strminKspring,'N/m'));

prompt = strcat('Fmax は ',num2str(h+deltaXoInput),'mm 以上です。');
disp('図を参照して Fmax を選んでください：');
Fmax = input(prompt);

prompt = strcat('Nmax は ',num2str(Fmax*minAcceptableKspring/1000),'N 以上です。');
disp('図を参照して　Nmax を選んでください：');
Nmax = input(prompt);

selectedKspring = Nmax/Fmax*1000;

fSpring = -selectedKspring.*deltaX;
% plot(fSpring); % N

fCutMax = input("最大切断力を入力してください:");
if isempty(fCutMax) || (fCutMax == 0)
    disp("切断操作なし。");
    fCut = 0;
else
    thickness = input("切断厚みを入力してください:");
    if h < 0
        disp("エラー：ストロークがネガティブです。");
        return
    else
        cutStartAngleIndex = find(displacement>rBase+rRoller, 1 )-1;
        % find the first index at which displacement bigger than rBase,
        % which means after cutting process start, minus 1 so that the
        % index become that of the position just before cutting
        cutEndAngleIndex = find(displacement > rBase+rRoller + thickness,1);
        % the index at which cutting process ends
        fCut = zeros(1,length(displacement));
        % Initialize a vector of size equal to that of displacement with all elements set to 0
        fCut(cutStartAngleIndex:cutEndAngleIndex) = fCutMax;  
        % Set the elements corresponding to position during cutting to max cutting force
    end
end

parallelForce = m*acceleration - fSpring - fFriction + fCut; % in N
perpendicularForce = parallelForce.*tanPressureAngle;
motorTorque = displacement/1000.*perpendicularForce;
figure;
plot(theta,motorTorque);
xlabel('角度') ; % misumi nomenclature
ylabel('トルク (Nm)') ;
torqueTitle = strcat('最大トルク  ', num2str(max(motorTorque)),' Nm');
title(torqueTitle,'Color','b','FontSize',15,'FontWeight','light');
grid on; grid minor;
disp(torqueTitle);

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
% ====================================
% ANIMATING CAM ROTATION
% ====================================

angleLineFactor = 4;
[angleEndPointX,angleEndPointY] = offsetIn(pitchX,pitchY,-angleLineFactor*rRoller);


% Cam motion simulation
prompt = "Show cam motion? カムの動きのシミュレーションを表示しますか? y/n [n]: ";
txt = input(prompt,"s");
if (txt == 'y')
    prompt = "Show cam motion in separate windows with parameters displayed? y/n [n]: ";
    separate = input(prompt,"s");
end
% Animate machining process if user input 'y'
% otherwise (input 'n' or just Enter), skip
if (txt == 'y')
    if (separate == 'y')
        figure('Name','カムの動きのシミュレーション');
    else
        figure;
    end
for loopNumber = 1:1
    
    hold on
    plot(0,0,'o','MarkerFaceColor','r');
    % Pitch Circle
    rotatedPitch = rotateCw([pitchX;pitchY],-pi/2);
    pl = plot(rotatedPitch(1,:),rotatedPitch(2,:),'color',pitchColor);
    axis equal;
    maxDim = rBase + abs(h) + 2*rRoller + 10;
    xlim([-maxDim maxDim]);
    ylim([-maxDim maxDim]);
    
    hold on;
    grid on;
    grid minor;
    pl.XDataSource = 'xx';
    pl.YDataSource = 'yy';
    
    
    % Cam Profile
    rotatedCam = rotateCw([camSurfX;camSurfY],-pi/2);
    pl2 = plot(rotatedCam(1,:),rotatedCam(2,:),'color',camColor);
    axis equal;
    maxDim = rBase + abs(h) + 2*rRoller + 10;
    xlim([-maxDim maxDim]);
    ylim([-maxDim maxDim]);
    
    hold on;
    pl2.XDataSource = 'xx2';
    pl2.YDataSource = 'yy2';
    
    % Roller
    
    index = linspace(0,2*pi,100);
    xC = rRoller*cos(index);
    yC = rRoller*sin(index) + displacement(1);
    pl3 = plot(xC,yC,'color',rollerColor);
    pl3.YDataSource = 'yC';
    
    rollerCenterY = displacement(1);
    pl4 = plot(0,rollerCenterY,'o','MarkerFaceColor',[0 0.4470 0.7410]); % Roller center
    pl4.YDataSource = 'rollerCenterY';
    
    maxDim = rBase + abs(h) + 2*rRoller + 10;
    xlim([-maxDim maxDim]);
    ylim([-maxDim maxDim]);
    hold on;
    
    % Pressure angle
    rollerCenterY_angle = displacement(1)+angleLineFactor*rRoller;
    angle1y = [rollerCenterY rollerCenterY_angle];
    pl5 = plot([0 0],angle1y,'MarkerFaceColor',[0 0.4470 0.7410]); 
    pl5.YDataSource = 'angle1y';
    hold on;
    
    rotatedAngleEnd = rotateCw([angleEndPointX(1);angleEndPointY(1)],-pi/2);
    angle2x = [0 rotatedAngleEnd(1)];
    angle2y = [rollerCenterY rotatedAngleEnd(2)];
    pl6 = plot(angle2x,angle2y,'MarkerFaceColor',[0 0.4470 0.7410]); 
    pl6.XDataSource = 'angle2x';
    pl6.YDataSource = 'angle2y';
    
    % Animation 
    for i = 1:length(theta2)
    j = theta2(i)-pi/2;
    
    rotatedPitch = rotateCw([pitchX;pitchY],j);
    xx = rotatedPitch(1,:);
    yy = rotatedPitch(2,:);
    
    rotatedCam = rotateCw([camSurfX;camSurfY],j);
    xx2 = rotatedCam(1,:);
    yy2 = rotatedCam(2,:);
    
    yC = rRoller*sin(index) + displacement(i);
    rollerCenterY = displacement(i);
    rollerCenterY_angle = displacement(i)+angleLineFactor*rRoller;
    angle1y = [rollerCenterY rollerCenterY_angle];
    
    rotatedAngleEnd = rotateCw([angleEndPointX(i);angleEndPointY(i)],j);
    angle2x = [0 rotatedAngleEnd(1)];
    angle2y = [rollerCenterY rotatedAngleEnd(2)];
    
    if (separate == 'y')
    temp5 = strcat('圧角　',num2str(pressureAngle(i)),'^o     ');
    temp2 = strcat('変位　',num2str(displacement(i)-rPrime),' mm     ');
    temp3 = strcat('回転角度　',num2str(theta(i)),'^o   ');
    updatedTitle = {temp3; temp2; temp5};
    [titleAni,] = title(updatedTitle,'Color',[0 0.4470 0.7410],'FontSize',14);
    
    temp4 = strcat('位置　',num2str(displacement(i)),' mm     ');
    ylabel(temp4,'Color',angleColor,'FontSize',15);
    temp1 = strcat('経過時間　',num2str(time(i)),' s     ');
    xlabel(temp1,'Color',angleColor,'FontSize',15);
    end

    refreshdata
    pause(0.001)
    end
end
end

%% 
%============================================
% Save all data 
%============================================
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
