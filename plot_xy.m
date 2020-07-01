function stop = plot_xy(x, ~, state)

stop = false;

switch state
    
    case 'init'
        cla;
        ax = gca;
        ax.XLabel.String = 'Heat extraction rate [W]';
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
        plot(x_data, y_data, 'ro-');
        ax.UserData = [x_data; y_data];
        
end
