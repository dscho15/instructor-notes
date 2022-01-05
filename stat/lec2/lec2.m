%% lec2
clc; clear;
format shortG

%% task 1
load('dataset_problem_2_1.mat');

% the data is X
p = 5
x_bar = mean(X) % the mean axes are about the same
S = cov(X)
R = corrcoef(X) % we have negative correlation between most of the data

figure(1)
boxplot(X)
figure(2)
plotmatrix(X) 

% the assumption of normality seems to hold, but we're still gonna perform
% a Mahalanobis QQ-plot
d = diag((X - x_bar) * inv(S) * (X - x_bar)')
rng('default') % For reproducibility
y = chi2rnd(p, 10000, 1);
figure(3)
qqplot(d, y)
% very normal distributed


%% Task 2
% Generation of general multivariate normal distribution from 
% standard univariate normal distribution
clc; clear;

mu = [1; 2; 0]
sigma = [[4 1 2]; [1 9 -3]; [2 -3 5]]
[P, V] = eig(sigma) % P V P'

% now we want to show that we can generate a normal distribution
% from (mu, sigma)

Y  = randn(1000, 3);
X = (mu + sqrtm(sigma) * Y')'
sigma_X = cov(X)
mu_X = mean(X)

figure(1)
plotmatrix(X)

% we have shown that we are able to generate normal distributed dat
p = 3
d = diag((X - mu_X) * inv(sigma_X) * (X - mu_X)')
y = chi2rnd(p, 10000, 1);
figure(2)
qqplot(d, y)

%% 
cr_ellipsis_2d((x_bar_female - x_bar_male)', S_pool, alpha, n1, n2)

% simultaneous
sim_bounds = x_bar' + [-sim_cor sim_cor]
hold on;
plot([sim_bounds(1, 1) sim_bounds(1, 1) sim_bounds(1, 2) sim_bounds(1, 2) sim_bounds(1, 1)], ...
     [sim_bounds(2, 1) sim_bounds(2, 2) sim_bounds(2, 2) sim_bounds(2, 1) sim_bounds(2, 1)], ...
     'DisplayName','Simultaneous')

% marginal 
marg_cor = tinv(1-alpha/2, n-1) * sqrt(diag(S)/n)
marg_bounds = x_bar' + [-marg_cor marg_cor]
hold on;
plot([marg_bounds(1, 1) marg_bounds(1, 1) marg_bounds(1, 2) marg_bounds(1, 2) marg_bounds(1, 1)], ...
     [marg_bounds(2, 1) marg_bounds(2, 2) marg_bounds(2, 2) marg_bounds(2, 1) marg_bounds(2, 1)], ...
     'DisplayName','Marginal')

% bonf 
marg_cor = tinv(1-alpha/(2*p), n-1) * sqrt(diag(S)/n)
marg_bounds = x_bar' + [-marg_cor marg_cor]
hold on;
plot([marg_bounds(1, 1) marg_bounds(1, 1) marg_bounds(1, 2) marg_bounds(1, 2) marg_bounds(1, 1)], ...
     [marg_bounds(2, 1) marg_bounds(2, 2) marg_bounds(2, 2) marg_bounds(2, 1) marg_bounds(2, 1)], ...
     'DisplayName','Bonferroni')

legend()
disp("This is confidence region for the mean")


