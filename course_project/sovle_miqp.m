function [ objective_value, allocation_vector, quantity_vector ] = sovle_miqp(mean_vector, covariance_matrix, price_vector, real_budget, scale_factor, min_stocks, max_stocks, min_fraction, max_fraction, min_fraction_invested, lambda, max_sane_return )

p = scale_factor*price_vector';
N = size(mean_vector,2);
r = mean_vector';
Q = covariance_matrix;

total_budget = real_budget*scale_factor;
M = max_stocks;
m = min_stocks;
fmin = min_fraction*total_budget;
fmax = max_fraction*total_budget;

xvars = 1:N; %The price of each share
vvars = N+1:2*N; %Indicator whether share bought or not
mvars = 2*N+1:3*N; %The number of each share
zvar = 3*N+1;





% The lower bounds of all the 2N+1 variables in the problem are zero. The upper bounds of the first 2N variables are one, and the last variable has no upper bound.
lb = zeros(3*N+1,1);
ub = [total_budget*ones(N,1);ones(N,1);Inf*ones(N,1)];
ub(zvar) = Inf;

% lb
% ub
% 
% Set the number of assets in the solution to be between 100 and 150. Incorporate this constraint into the problem in the form, namely
% 
% by writing two linear constraints of the form :



A = zeros(1,3*N+1); % Allocate A matrix
A(vvars) = 1; % A*x represents the sum of the v(i)
A = [A;-A];
b = zeros(2,1); % Allocate b vector
b(1) = M;
b(2) = -m;

% Include semicontinuous constraints. Take the minimal nonzero fraction of assets to be 0.001 for each asset type, and the maximal fraction to be 0.05.


% Include the inequalities  and  as linear inequalities.
Atemp = eye(N);
Amax = horzcat(Atemp,-Atemp*fmax,zeros(N,N+1));


A = [A;Amax];
b = [b;zeros(N,1)];
Amin = horzcat(-Atemp,Atemp*fmin,zeros(N,N+1));

% Amin
A = [A;Amin];
b = [b;zeros(N,1)];

% 
A = [A; [ones(1,N),zeros(1,2*N+1)]];
A = [A; [-ones(1,N),zeros(1,2*N+1)]];
b = [b; total_budget; -min_fraction_invested*total_budget];




  

% Include the constraint that the portfolio is 100% invested, meaning .
Aeq = zeros(N,3*N+1);
for k = 1:N
    Aeq(k,k) = 1;
    Aeq(k,2*N+k) = -p(k);
end
% Aeq

beq = zeros(N,1);


% Set the risk-aversion coefficient  to 100.


% Define the objective function  as a vector. Include zeros for the multipliers of the  variables.
% r
% zeros(2*N,1)
% lambda

f = [-r;zeros(2*N,1);lambda];

% f

% Solve the Problem
% To solve the problem iteratively, begin by solving the problem with the current constraints, which do not yet reflect any linearization. The integer constraints are in the vvars vector.
options = optimoptions(@intlinprog,'Display','off'); % Suppress iterative display
% size(Aeq)
% 
% size(A)
% 
% size(f)
[xLinInt,fval,exitFlagInt,output] = intlinprog(f,[vvars;mvars],A,b,Aeq,beq,lb,ub,options);


% Prepare a stopping condition for the iterations: stop when the slack variable  is within 0.01% of the true quadratic value.
thediff = 1e-2;
iter = 1; % iteration counter
assets = xLinInt(xvars); % the x variables

truequadratic = assets'*Q*assets;
zslack = xLinInt(zvar); % slack variable value
xLinInt(1:N)'*Q*xLinInt(1:N);

% Keep a history of the computed true quadratic and slack variables for plotting.
history = [truequadratic,zslack];

% Compute the quadratic and slack values. If they differ, then add another linear constraint and solve again.
% In toolbox syntax, each new linear constraint  comes from the linear approximation
% 
% You see that the new row of  and the new element in  , with the  term represented by a -1 coefficient in .
% After you find a new solution, use a linear constraint halfway between the old and new solutions. This heuristic way of including linear constraints can be faster than simply taking the new solution. To use the solution instead of the halfway heuristic, comment the "Midway" line below, and uncomment the following one.
numIter = 0;
maxIter = 100;
while abs((zslack - truequadratic)/truequadratic) > thediff  && numIter < maxIter% relative error
%     numIter
%     zslack
%     truequadratic
    err = abs((zslack - truequadratic)/truequadratic);
    newArow = horzcat(2*assets'*Q,zeros(1,2*N),-1); % Linearized constraint
  
    
    A = [A;newArow];
    b = [b;truequadratic];
    % Solve the problem with the new constraints
    [xLinInt,fval,exitFlagInt,output] = intlinprog(f,[vvars;mvars],A,b,Aeq,beq,lb,ub,options);
%     assets = (0.1*assets+0.9*xLinInt(xvars)); % Midway from the previous to the current
    assets = (assets+xLinInt(xvars))/2; % Midway from the previous to the current

%     assets = xLinInt(xvars); % Use the previous line or this one
    truequadratic = assets'*Q*assets;
    zslack = xLinInt(zvar);
    history = [history;truequadratic,zslack];
    iter = iter + 1;
    numIter = numIter + 1;
end
% 
% Examine the Solution and Convergence Rate
% Plot the history of the slack variable and the quadratic part of the objective function to see how they converged.
% % disp("History")
% % disp(history)
% h =figure; plot(history)
% legend('Quadratic','Slack', 'location', 'NorthEast')
% xlabel('Iteration number')
% title('Quadratic and linear approximation (slack)')
% saveas(h, 'C:/Users/vignesh.subramanian/Desktop/stockPredictor/plots/slackconv', 'png')

% What is the quality of the MILP solution? The output structure contains that information. Examine the absolute gap between the internally-calculated bounds on the objective at the solution.
% disp(output.absolutegap)
% The absolute gap is zero, indicating that the MILP solution is accurate.

% Plot the optimal allocation. Use xLinInt(xvars), not assets, because assets might not satisfy the constraints when using the midway update.
% % figure; bar(xLinInt(xvars)/total_budget)
% % grid on
% % xlabel('Asset index')
% % ylabel('Proportion of investment')
% % title('Optimal asset allocation')

% You can easily see that all nonzero asset allocations are between the semicontinuous bounds  and .
% How many nonzero assets are there? The constraint is that there are between 100 and 150 nonzero assets.
% sum(xLinInt(vvars))

% xLinInt(mvars)


% xLinInt(xvars)


quantity_vector = xLinInt(mvars);
allocation_vector = xLinInt(xvars)/(total_budget*scale_factor);
objective_value = -fval;

% What is the expected return for this allocation, and the value of the risk-adjusted return?
fprintf('The expected return is %g, and the risk-adjusted return is %g.\n',...
    r'*xLinInt(xvars),-fval)
% More elaborate analyses are possible by using features specifically designed for portfolio optimization in Financial Toolbox™.
% Copyright 2015 The MathWorks, Inc.


end

