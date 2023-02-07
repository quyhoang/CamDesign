clc; close all; clear;
%============================================
% INPUT 入力
%============================================
% All values in degree

% Job-specific values
% eventAngle = [rise start - rise end - return start- return end]
eventAngle = [80 120 190 230]; % degree at which the rise/return starts/ends
h = 15.2; % stroke in mm
RPM = 200; % motor velocity in rounds per minutes
m = 1; % follower mass in kg
rRoller = 8;

camName = input("設定データ名を入力してください。","s");
if ~isempty(camName)
    save(camName)
    disp(['設定データは全て', camName,'.mat に保存されています。']);
else
    disp('入力がないので、データは保存されていません。');
end
