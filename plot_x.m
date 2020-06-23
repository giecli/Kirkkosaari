function stop = plot_x(x, optim_values, state)

stop = false;

switch state
    
    case 'init'
        cla;
        ax = gca;
        ax.XLabel.String = 'Iteration';
        ax.YLabel.String = 'Heat extraction rate [W]';
        ax.UserData = [];
        
    case 'iter'
        ax = gca;
        if isempty(ax.UserData)
            x_data = [optim_values.iteration];
            y_data = [x];
        else
            x_data = [ax.UserData(1,:), optim_values.iteration];
            y_data = [ax.UserData(2,:), x];
        end
        plot(x_data, y_data, 'ro-');
        ax.UserData = [x_data; y_data];
        
end
