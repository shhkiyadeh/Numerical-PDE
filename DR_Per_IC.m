% Sanaz Hami
% Douglas-Rashford Function for 1D Linear Advection with Periodic IC

function [u_mod] = DR(m, M, u, u_guess, epsilon, max_iter, lambda)
    N = length(u); % Number of elements in u
    b = sum(u); % Sum of u, which should be preserved
    c = b / N; % Average of the elements in u

    % Initialize the modified solution with the initial guess
    u_mod = u_guess;
    
    % Auxiliary variables for the DRS method
    y = u_mod; % Initialization of the auxiliary variable y
    z = zeros(size(u_mod)); % Initialization of the auxiliary variable z

    for iter = 1:max_iter
        % Step 1: Reflection
        y_reflected = 2 * u_mod - y;
        
        % Projection onto the set [m, M]
        z = min(max(y_reflected, m), M);
        
        % Step 2: Reflection
        z_reflected = 2 * z - y_reflected;
        
        % Projection onto the hyperplane {x | sum(x) = b}
        y = z_reflected + (b - sum(z_reflected)) / N;
        
        % Check for convergence: if the norm is less than epsilon, stop
        if norm(y - u_mod, 2) <= epsilon
            disp(['DRS converged in ', num2str(iter), ' iterations.']);
            break;
        end
        
        % Update for the next iteration
        u_mod = (1 - lambda) * u_mod + lambda * y;
    end

    if iter == max_iter
        warning('DRS did not converge within the maximum number of iterations.');
        u_mod = u_guess; % If not converged, return the initial guess
    end

    % The last step ensures that the sum of u_mod equals the sum of the original u
    
    % Adjust u_mod to satisfy the conservation property
    u_mod = u_mod + (b - sum(u_mod)) / N;
end
