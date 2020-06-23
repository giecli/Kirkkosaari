function stop = plot_fval(x, optim_values, state)

stop = false;

switch state
    
    case 'init'
        cla;
        ax = gca;
        ax.XLabel.String = 'Iteration';
        ax.YLabel.String = 'Average borehole wall temperature [\circC]';
        ax.UserData = [];
        
    case 'iter'
        ax = gca;
        if isempty(ax.UserData)
            x_data = [optim_values.iteration];
            y_data = [optim_values.fval];
        else
            x_data = [ax.UserData(1,:), optim_values.iteration];
            y_data = [ax.UserData(2,:), optim_values.fval];
        end
        semilogy(x_data, y_data, 'ro-');
        ax.UserData = [x_data; y_data];
        
end
