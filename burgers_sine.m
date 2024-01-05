function [u] = burgers_sine(x, t)
    f = @(y) x - sin(y) * t - y;
    options = optimset('Display', 'off');
    u = fsolve(f, 0, options);
end
