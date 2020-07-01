function stop = plot_xy(x, ~, state, ~)
% function stop = plot_xy(varargin)

fprintf(1, 'x(1)=%f x(2)=%f\n', x(1), x(2));

stop = false;

switch state
    
    case 'init'
        cla;
        ax = gca;
        ax.XLabel.String = 'Heat extraction rate [kW]';
        ax.YLabel.String = 'Flow rate [L/s]';
        ax.UserData = [];
        
    case 'iter'
        ax = gca;
        if isempty(ax.UserData)
            x_data = [x(1)];
            y_data = [x(2)];
        else
            x_data = [ax.UserData(1,:), x(1)];
            y_data = [ax.UserData(2,:), x(2)];
        end
        plot(x_data/1000, y_data, 'ro-');
        ax.UserData = [x_data; y_data];
        
end
