function [alt_objective_value, alt_allocation_vector, alt_quantity_vector] = solve_qp(mean_vector, covariance_matrix, price_vector, current_budget, min_stocks, max_stocks, min_fraction, max_fraction, lambda, max_sane_return ) 

n = size(mean_vector,2);
cvx_begin
    variable x(n);
    minimize( -mean_vector*x + lambda*quad_form(x, covariance_matrix)); 
    
    subject to
        0 <= x <= max_fraction*current_budget
        sum(x)<= current_budget
cvx_end

x_cleaned = [];
for i = 1:size(x,1)
    if x(i) <= min_fraction*current_budget
        x_cleaned = [x_cleaned; 0];
    else
        x_cleaned = [x_cleaned; x(i)];
    end
end

%Fix based on integer stock quantity

% quantity_vector = floor(x_cleaned./price_vector');
% 
% x_realistic = quantity_vector.*price_vector';


%Sort and only keep max number of stocks
[sorted_x_realistic, I] = sort(x_cleaned, 'descend');

x_realistic = x_cleaned;
x_realistic(I(max_stocks+1:end)) = 0;

if sum(x_realistic) > 0
    x_realistic = min(current_budget*x_realistic/sum(x_realistic), max_fraction*current_budget);
end

quantity_vector = floor(x_realistic./price_vector');

x_realistic = quantity_vector.*price_vector';



alt_quantity_vector = x_realistic./price_vector';
alt_allocation_vector = x_realistic/current_budget*100;
alt_objective_value = cvx_optval;
% figure; bar(alt_allocation_vector)



end

