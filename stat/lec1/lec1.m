%% Sigma 
format shortG

sigma = [[4 1 2];
         [1 9 -3];
         [2 -3 25]]

% Perform the eigenvalue decomposition
[V, D] = eig(sigma)
% V*D*V'

% Verify that the determinant is the sum of eig-vals
disp("prod of eig: " + string(prod(diag(D))))
disp("det: " + string(det(sigma)))
% thus they are equal to eachother

% Verify that the trace is the sum of eig-vals
disp("sum of sigma: " + string(trace(sigma)))
disp("trace of sigma: " + string(trace(sigma)))

%% Multivariate visualization and descriptive statistics

% We shall calculate the descriptive statistics for a datset

load('dataset_problem_1_2.dat');
x = dataset_problem_1_2;
x_bar = mean(x) 
S = cov(x)
R = corrcoef(x)

% x_bar - dependents on unit of time and gives the average time on each
% axis.
% S - describes the covariance between the different axes's, it doesn't
% really make intuitively sense, since no normalization is made.
% R - describes the correlation covariance, and is basically a normalized
% variance, here  we can see that correlation between a marathon and a 100
% meters is low.

% https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Box_Plots.pdf

figure(1)
boxplot(x, 'Whisker', 1) % three outliers on boxplot / plotmatrix
% we have the 25% percentile
% we have the 75% percentile
% the iqr is the measurement between 75% - 25%
% 25/75 percentile + pm 1.5 * (iqr)

figure(2)
plotmatrix(x)

% Calculate the generalized sample variances |S| and |R| and verify their
% expression |S| = prod of s_ii |R|

disp("|S| " + string(det(S)))
disp("|R| " + string(det(R)))

% |R| is more "saying", but the number is still not intuitive at it is the
% volume of the variances, or "they span a volume with coordinates axis's".


%%

