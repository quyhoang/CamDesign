function combination = writecsv(varargin)
    % Combine and export to csv
    % First input is the name of the csv file

    filename = [varargin{1}, '.csv'];
    
    % Initialize
    index = linspace(0,360,length(varargin{2}));
    combination = zeros(nargin,length(index));
    combination(1,:) = index;

    % Loop through each input argument
    for k = 2:nargin 
        combination(k,:) = varargin{k};
    end

    % Write
    writematrix(combination', filename);
    disp(['Data has been written to ', filename]);
end
