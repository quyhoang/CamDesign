classdef camdata

    properties
        name = 'uninitilized'
        unit = 'uninitilized'
        data = 0
    end

    methods
        function obj = camdata(name, unit, data)
            if nargin > 0
                obj.name = name;
            end
            if nargin > 1
                obj.unit = unit;
            end
            if nargin > 2
                obj.data = data;
            end
        end

        function show(obj)
            % Plot data

            figure;
            plot(linspace(0,360,length(obj.data)),obj.data);
            grid on; grid minor;
            xlim([0 360]);
            xlabel({'回転角度','degree'},'FontSize',15,'FontWeight','light','Color','b');
            ylabel({obj.name, obj.unit},'FontSize',15,'FontWeight','light','Color','b');
            temp = strcat('最大',obj.name,' :  ','',num2str(max(obj.data)),' ',obj.unit);
            [t,s] = title({'',obj.name},{temp,''},'Color','blue');
            t.FontSize = 20;
            s.FontSize = 16;
            s.FontAngle = 'italic';
        end

        function export(obj, filename)
            index = linspace(0,360,length(obj.data));
            exportdata = zeros(2,length(index));
            exportdata(1,:) = index;
            exportdata(2,:) = obj.data;
            exportdata = exportdata';

            columnHeaders = {'回転角度', obj.name};
            T = array2table(exportdata, 'VariableNames', columnHeaders);
            filename = [filename, '.xlsx'];
             writetable(T, filename)
%             writetable(T, filename, 'Encoding', 'UTF-8');
            disp(['Data has been written to ', filename]);
        end
    end
end