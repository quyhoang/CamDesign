function combination = combine(varargin)
    % Calculate sum and plot over the period of 360 degree
    
    % Initialize the sum to 0
    total = 0;
    
    index = linspace(0,360,length(varargin{1}));
    combination = zeros(nargin+2,length(index));
    combination(1,:) = index;

    figure; hold on;
    % Loop through each input argument and add it to the sum
    for k = 1:nargin
        combination(k+2,:) = varargin{k};
        total = total + varargin{k};
        plot(index,varargin{k}); 
    end

    combination(2,:) = total;

    % Plot result
    plot(index,total,'LineWidth',2.5);
    xlabel('角度') ;
    xlim([0 360]);
    plottitle = strcat('最大 ', num2str(max(total)));
    title(plottitle,'Color','b','FontSize',15,'FontWeight','light');
    grid on; grid minor;
    legend;
    disp(plottitle);

    combination = combination';
end