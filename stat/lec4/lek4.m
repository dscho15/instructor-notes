%% Task 1
clc; clear;
format short
load('dataset_problem_4_1.mat')

% we consider a multivariate normal distribution
% and iid observations

%% Task 2
clc; clear;
format short
load('dataset_problem_4_2.mat');

% All measurements are in units [mm] and columns correspods to tail len
% and wing length

% It is desired to test for equality of mean vectors for the population of
% male and female hook-billd kites

% descriptive statistics
x_bar_male = mean(X_male);
x_bar_female = mean(X_female);
S_male = cov(X_male);
S_female = cov(X_female);

% perform a bartlett test, to investigate whether we can assume
% homoscedatiscity 

alpha = 0.05;
p = 2;
n1 = 45;
n2 = 45;
C = 1 - (2*p^2 + 3*p-1)/(6*(p+1)) * (1/(n1-1) + 1/(n2-1) - (1/(n1+n2-2)));
S_pool = ((n1-1) * S_male + (n2-1) * S_female) / (n1+n2-2);

test_statistics = C*((n1+n2-2)*log(det(S_pool)) - (n1-1)*log(det(S_male)) - (n2-1)*log(det(S_female)));
critical_value = chi2inv(1-alpha, p*(p+1)/2);

disp("Reject [y/n]: " + string(test_statistics > critical_value))
disp("p-value: " + string(1 - chi2cdf(test_statistics, p*(p+1)/2)))

% seems like we can assume normality
figure(1)
boxplot([X_male, X_female]);

rng('default') % For reproducibility
y = chi2rnd(p, 10000, 1);

% mahalanobis
figure(2)
subplot(2, 1, 1);
d = diag((X_male - x_bar_male) * inv(S_pool) * (X_male - x_bar_male)');
qqplot(d, y);
pbaspect([1 1 1]);
subplot(2, 1, 2);
d = diag((X_female - x_bar_female) * inv(S_pool) * (X_female - x_bar_female)');
qqplot(d, y);
pbaspect([1 1 1]);
% seems like it is "okay" normally distributed

% Test for equality of the mean vectors for female and male populations
% Perform an exact test of the hypothesis H_0. Use level of significance alpha = 0.05 but 
% calculate also the p-value for the test:
test_statistics = (x_bar_female - x_bar_male) * inv( (1/n1 + 1/n2)  * S_pool) * (x_bar_female - x_bar_male)';
critical_value = (n1+n2-2)*p/(n1+n2-p-1)* finv(1-alpha, p, n1+n2-p-1);
p_value = 1 - fcdf((n1+n2-p-1)/((n1+n2-2)*p) * test_statistics , p, n1+n2-p-1);
disp("exact test")
disp("p_value " + string(p_value))
critical_value = chi2inv(1-alpha, p);
p_value = 1 - chi2cdf(test_statistics, p);
disp("large-scale test")
disp("p_value " + string(p_value))

%% Confidence regions and intervals

figure(3)

cr_ellipsis_2d((x_bar_female - x_bar_male)', S_pool, alpha, n1, n2)
pbaspect([1 1 1]);

% simultaneous
sim_cor = sqrt((n1+n2-2)*p/(n1+n2-p-1) * finv(1-alpha, p, n1+n2-p-1)) * sqrt(diag(S_pool)*(1/n1 + 1/n2))
sim_bounds = (x_bar_female - x_bar_male)' + [-sim_cor sim_cor]
hold on;
plot([sim_bounds(1, 1) sim_bounds(1, 1) sim_bounds(1, 2) sim_bounds(1, 2) sim_bounds(1, 1)], ...
     [sim_bounds(2, 1) sim_bounds(2, 2) sim_bounds(2, 2) sim_bounds(2, 1) sim_bounds(2, 1)], ...
     'DisplayName','Simultaneous')

% marginal 
marg_cor = tinv(1-alpha/2, n1+n2-2) * sqrt(diag(S_pool)*(1/n1 + 1/n2))
marg_bounds = (x_bar_female - x_bar_male)' + [-marg_cor marg_cor]
hold on;
plot([marg_bounds(1, 1) marg_bounds(1, 1) marg_bounds(1, 2) marg_bounds(1, 2) marg_bounds(1, 1)], ...
     [marg_bounds(2, 1) marg_bounds(2, 2) marg_bounds(2, 2) marg_bounds(2, 1) marg_bounds(2, 1)], ...
     'DisplayName','Marginal')

% bonf 
marg_cor = tinv(1-alpha/(2*p), n1+n2-2) * sqrt(diag(S_pool)*(1/n1 + 1/n2))
marg_bounds = (x_bar_female - x_bar_male)' + [-marg_cor marg_cor]
hold on;
plot([marg_bounds(1, 1) marg_bounds(1, 1) marg_bounds(1, 2) marg_bounds(1, 2) marg_bounds(1, 1)], ...
     [marg_bounds(2, 1) marg_bounds(2, 2) marg_bounds(2, 2) marg_bounds(2, 1) marg_bounds(2, 1)], ...
     'DisplayName','Bonferroni')

legend()
disp("This is confidence region for the mean")



