function [ returns, sharpe_ratio ] = get_stats( budget_vector, days_in_period )

overall_returns = (budget_vector(end)/budget_vector(1))
total_days = size(budget_vector,1)*days_in_period;
number_years = total_days/365;

annualized_return = 100*(overall_returns^(1/number_years) - 1);
returns = annualized_return;

returns_vector = [];
for k = 2:size(budget_vector,1)
    period_return = budget_vector(k)/budget_vector(k-1) - 1;
    returns_vector = [returns_vector, period_return];
end
mean_return = 100*mean(returns_vector)
std_return = 100*std(returns_vector)
sharpe_ratio = mean_return/std_return;
end

