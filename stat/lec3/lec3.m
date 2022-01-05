clc; clear;
format shortG

%%
clc;
load('dataset_problem_3_1.mat')

X = obs;

% The first column (x1) contains measurements of stiffness (modulus of 
% elasticity) and the second column (x2) contains the corresponding 
% measurements of bending strength for the lumber pieces,
% both values in units proportional to [N/m2].

figure(1)
plotmatrix(X)
figure(2)
boxplot(X)

% looks fairly normal distributed

x_bar = mean(X);
S = cov(X);
r = corrcoef(X);

% we want to perform a test, whether the normal distribution is 
% [1750, 950], we want to use an exact test.

alpha = 0.05;
n = 30;
p = 2;
mu_0 = [1750, 950];
test_statistics = (x_bar - mu_0) * inv(S / n) * (x_bar - mu_0)';
%
critical_value = p*(n-1)/(n-p) * finv(1 - alpha, p, n-p);
H_o_rejection = test_statistics > critical_value;
P_value = 1 - fcdf((n-p)/(p*(n-1)) * test_statistics, p, n-p);
disp("Exact test")
disp("A binary rejection: " + string(H_o_rejection))
disp("The p-value: " + string(P_value))
disp("")
%
critical_value = chi2inv(1 - alpha, p);
H_o_rejection = test_statistics > critical_value;
P_value = 1 - chi2cdf(test_statistics, p);
disp("Large sample approximation")
disp("A binary rejection: " + string(H_o_rejection))
disp("The p-value: " + string(P_value))

%% Turning to confidence region and intervals using alpha = 0.05
figure(3)
cr_ellipsis(x_bar', S, alpha, n)
hold on;
scatter(mu_0(1), mu_0(2), 'DisplayName','H_o')
xlim([1600 2100])
ylim([700 1000])

sim_cor = sqrt(p*(n-1)/(n-p) * finv(1-alpha, p, n-p)) * sqrt(diag(S)/n)

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

% What can we say about the confidence regions?
% - simultaneous is pessimistic
% - marignal is optimistic
% - bonferroni is inbetween bonf and simu






