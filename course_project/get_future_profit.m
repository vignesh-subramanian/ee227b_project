function [future_profit] = get_future_profit(allocation_vector, current_budget, cur_week_index, history_weeks, stock_price_matrix, max_sane_return)
overall_profit = 0;
for i = 1:size(stock_price_matrix,2)
    future_week_index = cur_week_index + history_weeks;
    cur_stock_price = stock_price_matrix(cur_week_index, i);
    future_stock_price = stock_price_matrix(future_week_index, i);
    if cur_stock_price > 0
        future_stock_return = (future_stock_price/cur_stock_price) - 1;
    else
        future_stock_return = 0;
    end
    if isnan(future_stock_return)
        future_stock_return = 0;
    end
    if abs(future_stock_return) > max_sane_return
        future_stock_return = future_stock_return/abs(future_stock_return)*max_sane_return;
    end
    future_stock_profit = allocation_vector(i)*current_budget*future_stock_return/100;
    
    overall_profit = overall_profit + future_stock_profit;
    
end

future_profit = floor(overall_profit);
end

