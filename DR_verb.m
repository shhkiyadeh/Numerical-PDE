% Sanaz Hami
% Douglas-Rashford Splitting Function including Verbose Argument  

function [u_mod] = DR(m, M, u, u_guess, epsilon, max_iter, lambda, verbose)
    % Initialization
    N = length(u); % Number of elements in u
    b = sum(u); % Sum of u, which should be preserved
    c = b / N; % Average of the elements in u
    u_mod = u_guess; % Initial guess for the solution

    for iter = 1:max_iter
        % Douglas-Rachford iteration
        y = u_mod + lambda * (c - sum(u_mod) / N);

        % Apply the min/max limiter to y to get z
        z = max(min(y, M), m);

        % Update y
        y_next = u_mod + lambda * (z - u_mod);

        % Update u_mod
        u_mod = y_next;

        % Check for convergence
        if norm(u_mod - u_guess, 2) <= epsilon
            if verbose
                disp(['Converged in ', num2str(iter), ' iterations.']);
            end
            break;
        end
    end

    if iter == max_iter && verbose
        disp('Maximum iterations reached without convergence.');
    end

    % Enforce the sum constraint
    u_mod = u_mod * (b / sum(u_mod));
end
