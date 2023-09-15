function combinedtorque(varargin)
    % Calculate sum and torques over the period of 360 degree
    
    % Initialize the sum to 0
    
    lengthindex = length(varargin{1}.motorTorque.data);
    index = linspace(0,360,lengthindex);
    total = 0;
    
    figure; hold on;
    legendLabels = cell(1, nargin+1);
    % Loop through each input argument and add it to the sum
    for k = 1:nargin
        total = total + varargin{k}.motorTorque.data;
        plot(index,varargin{k}.motorTorque.data);
        legendLabels{k} = inputname(k);
    end   
    legendLabels{nargin+1} = 'Combined torque';
    
    % Plot result
    plot(index,total,'LineWidth',2.2);
    xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
    ylabel({'トルク','Nm'},'FontSize',15,'FontWeight','light','Color','b');
    xlim([0 360]);
    plottitle = strcat('最大モータートルク ', num2str(max(total)),' Nm');
    title(plottitle,'Color','b','FontSize',15,'FontWeight','light');
    grid on; grid minor;
    legend(legendLabels);
    disp(plottitle);
end