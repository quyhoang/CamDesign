function output = excel(varargin)
    % Combine and export to xlsx
    % First input is the name of the xlsx file
    % Other inputs are camdata objects
    
    % Initialize
    index = linspace(0,360,length(varargin{2}.data));
    output = zeros(nargin,length(index));
    output(1,:) = index;

    % Create column header
    headers = cell(1,nargin);
    headers{1} = '回転角度';
    
    % Loop through each input argument
    for k = 2:nargin 
        output(k,:) = varargin{k}.data;
        headers{k} = varargin{k}.name;
    end
    output = output';

    filename = [varargin{1}, '.xlsx'];

    % Write to xlsx file
    T = array2table(output, 'VariableNames', headers);
    writetable(T, filename);
    disp(['Data has been written to ', filename]);
end