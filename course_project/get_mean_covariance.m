function [ mean_vector, covariance_matrix ] = get_mean_covariance( returns_matrix )
%Input: returns_matrix (M by N)
%M: Number of data entries for an asset
%N: Number of assets

mean_vector = mean(returns_matrix);
covariance_matrix = cov(returns_matrix);

end

