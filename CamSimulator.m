classdef CamSimulator
    properties
        camType % Oscillating or translating cam
        transition % Input transition matrix
        transition_angle % Angle values extracted from the transition matrix
        transition_displacement % Displacement values extracted from the transition matrix
        s_rad_initial % Initial angular displacement in radian
        l_roller % Distance from arm rotating axis to roller center
        l_load % Distance from arm center to load
        estLoad % Estimated load in Newton
        m_roller % Roller mass in kilograms
        rRoller % Radius of the roller
        m_rocker % Rocker arm mass in kilograms
        rocker2cam % Distance between rocker arm axis and cam axis
        RPM % Motor velocity in rounds per minutes
        maxPressureAngle_deg % Max pressure angle in degree
        kFriction % Friction coefficient
        sampleRate % Rate at which the roller on pitch curve is sampled
        step % Sampling rate in degrees for calculation
        rollerColor % Color to use for the roller in graphics
        pitchColor % Color to use for the pitch in graphics
        camColor % Color to use for the cam in graphics
        maxPressureAngle % Max pressure angle in radian
        theta % Range of angles from 0 to 360 with step size as defined by 'step'
        T % Period of moving 360 degrees, in seconds
        time % Array of time values corresponding to each angle in 'theta'
        timeStep % Convert step in degree to step in time
    end
    
    methods
        function obj = CamSimulator() % Constructor method
            obj.transition = [30 0; 130 3.5; 260 3.5; 350 0];
            obj.transition_angle = obj.transition(:, 1)';
            obj.transition_displacement = obj.transition(:, 2)';
            obj.s_rad_initial = deg2rad(30);
            obj.l_roller = 65;
            obj.l_load = 52;
            obj.estLoad = 10;
            obj.m_roller = 0.1;
            obj.rRoller = 9.5;
            obj.m_rocker = 0.4;
            obj.rocker2cam = 94.2;
            obj.RPM = 200;
            obj.maxPressureAngle_deg = 20;
            obj.kFriction = 0.7;
            obj.sampleRate = 5;
            obj.step = .5;
            obj.rollerColor = [0.4660 0.6740 0.1880];
            obj.pitchColor = [0.8500 0.3250 0.0980];
            obj.camColor = 'b';
            obj.maxPressureAngle = deg2rad(obj.maxPressureAngle_deg);
            obj.theta = 0:obj.step:360;
            obj.T = 60/obj.RPM;
            obj.time = linspace(0,obj.T,length(obj.theta));
            obj.timeStep = obj.T/size(obj.time,2);
        end
        
        % Additional methods can be defined here for the operations you want to perform
    end
end
