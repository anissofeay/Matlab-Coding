# Matlab-Coding
This is the full matlab coding for this project

% BVP Solver for given ODEs using bvp4c
clc;
clear;

% Adjustable Parameters
Pr = 1;              % Prandtl number (adjust as needed)
gamma = 1;           % Parameter for boundary condition (adjust as needed)
etaMax = 10;         % Maximum value of eta to approximate infinity
etaMesh = linspace(0, etaMax, 100); % Mesh points for eta

% Define a range of epsilon values
epsilonValues = [2.0, 1.0, 0.5, 0]; % Adjust as needed

% Prepare results table
results = table();

for i = 1:length(epsilonValues)
    epsilon = epsilonValues(i); % Current epsilon value

    % Initial guess for the solution
    solinit = bvpinit(etaMesh, @guess);

    % Solve the BVP
    sol = bvp4c(@(eta, y) odefun(eta, y, Pr), ...
                @(ya, yb) bcfun(ya, yb, epsilon, gamma), ...
                solinit);

    % Evaluate f''(0) from the solution
    fDoublePrimeAtZero = sol.y(3, 1); % f''(0) = y(3) at eta = 0

    % Round f''(0) to 6 decimal places
    fDoublePrimeAtZero = round(fDoublePrimeAtZero, 6);

    % Append epsilon and f''(0) to results table
    results = [results; table(epsilon, fDoublePrimeAtZero)];
end

% Display the results with 6 decimal places using fprintf
fprintf('Results Table:\n');
fprintf('%10s %20s\n', 'epsilon', 'fDoublePrimeAtZero');
for i = 1:height(results)
    fprintf('%10.2f %20.6f\n', results.epsilon(i), results.fDoublePrimeAtZero(i));
end

% --- Standalone Functions ---

% ODE Function
function dydeta = odefun(eta, y, Pr)
    dydeta = zeros(5, 1);
    % y(1) = f, y(2) = f', y(3) = f''
    % y(4) = theta, y(5) = theta'
    dydeta(1) = y(2);                         % f' = y(2)
    dydeta(2) = y(3);                         % f'' = y(3)
    dydeta(3) = -y(1) * y(3) - 1 + y(2)^2;    % f'''
    dydeta(4) = y(5);                         % theta' = y(5)
    dydeta(5) = -Pr * y(1) * y(5);            % theta'' equation
end

% Boundary Condition Function
function res = bcfun(ya, yb, epsilon, gamma)
    res = zeros(5, 1);
    res(1) = ya(1);                            % f(0) = 0
    res(2) = ya(2) - epsilon;                  % f'(0) = epsilon
    res(3) = yb(2) - 1;                        % f'(infinity) = 1
    res(4) = ya(5) + gamma * (1 - ya(4));      % theta'(0) = -gamma * (1 - theta(0))
    res(5) = yb(4);                            % theta(infinity) = 0
end

% Initial Guess Function (Standalone)
function yinit = guess(eta)
    % Provide an initial guess for the solution
    yinit = [eta; 1; 0; exp(-eta); -exp(-eta)];
end
