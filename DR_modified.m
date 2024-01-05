function [u_mod] = DR_modified(u_bar, m, M, u_guess, epsilon, max_iter, lambda, verbose)
    % Initialization
    N = length(u_bar); % Number of elements in u_bar
    b = sum(u_bar); % Sum of u_bar, which should be preserved
    u_mod = u_guess; % Initial guess for the solution

    % DR iteration based on algorithm
    for iter = 1:max_iter
        % Step I: Projection onto the set [m, M]
        y = min(max(u_mod, m), M); 

        % Step II: Reflection step
        y_reflected = 2 * y - u_mod;

        % Step III: Projection onto the hyperplane defined by the conservation constraint 
        y_next = y_reflected - (sum(y_reflected) - b) / N;

        % Relaxation step
        u_mod = u_mod + lambda * (y_next - u_mod);

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
end
