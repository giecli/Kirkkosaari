function stop = print_x(x, ~, state)
stop = false;
switch state
    case 'init'
        fprintf(1, 'init: x=%s\n', x);
    case 'iter'
        fprintf(1, 'iter: x=%s\n', x);
    case 'done'
        fprintf(1, 'done: x=%s\n', x);
end
