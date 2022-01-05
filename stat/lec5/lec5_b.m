%% Task 1
clc; clear;
load('dataset_problem_5_2.mat')

%
mask_f1_l1 = (X(:, 1) == 1);
mask_f1_l2 = (X(:, 1) == 2);
mask_f2_l1 = (X(:, 2) == 1);
mask_f2_l2 = (X(:, 2) == 2);
mask_f2_l3 = (X(:, 2) == 3);

% x_bar
x_bar_total = mean(X(:, 3:end))';
% row-wise
x_bar_f1_l1 = mean(X(mask_f1_l1, 3:end))';
x_bar_f1_l2 = mean(X(mask_f1_l2, 3:end))';
% col-wise means
x_bar_f2_l1 = mean(X(mask_f2_l1, 3:end))';
x_bar_f2_l2 = mean(X(mask_f2_l2, 3:end))';
x_bar_f2_l3 = mean(X(mask_f2_l3, 3:end))';
% interaction
x_bar_int_l1_l1 = mean(X(mask_f1_l1 & mask_f2_l1, 3:end))';
x_bar_int_l1_l2 = mean(X(mask_f1_l1 & mask_f2_l2, 3:end))';
x_bar_int_l1_l3 = mean(X(mask_f1_l1 & mask_f2_l3, 3:end))';
x_bar_int_l2_l1 = mean(X(mask_f1_l2 & mask_f2_l1, 3:end))';
x_bar_int_l2_l2 = mean(X(mask_f1_l2 & mask_f2_l2, 3:end))';
x_bar_int_l2_l3 = mean(X(mask_f1_l2 & mask_f2_l3, 3:end))';

%
g = 2;
b = 3;
n = 2;

% total
SST = cov(X(:,3:end)) * (g*b*n-1);
% between 1
SSB1 = b*n * (x_bar_f1_l1 - x_bar_total) * (x_bar_f1_l1 - x_bar_total)' ...
       + b*n * (x_bar_f1_l2 - x_bar_total) * (x_bar_f1_l2 - x_bar_total)';
% between 2
SSB2 = g*n * (x_bar_f2_l1 - x_bar_total) * (x_bar_f2_l1 - x_bar_total)' ...
       + g*n * (x_bar_f2_l2 - x_bar_total) * (x_bar_f2_l2 - x_bar_total)' ...
       + g*n * (x_bar_f2_l3 - x_bar_total) * (x_bar_f2_l3 - x_bar_total)';
% interaction
SSINT =  n * (x_bar_int_l1_l1 - x_bar_f1_l1 - x_bar_f2_l1 + x_bar_total) * (x_bar_int_l1_l1 - x_bar_f1_l1 - x_bar_f2_l1 + x_bar_total)' ...
       + n * (x_bar_int_l1_l2 - x_bar_f1_l1 - x_bar_f2_l2 + x_bar_total) * (x_bar_int_l1_l2 - x_bar_f1_l1 - x_bar_f2_l2 + x_bar_total)' ...
       + n * (x_bar_int_l1_l3 - x_bar_f1_l1 - x_bar_f2_l3 + x_bar_total) * (x_bar_int_l1_l3 - x_bar_f1_l1 - x_bar_f2_l3 + x_bar_total)' ...
       + n * (x_bar_int_l2_l1 - x_bar_f1_l2 - x_bar_f2_l1 + x_bar_total) * (x_bar_int_l2_l1 - x_bar_f1_l2 - x_bar_f2_l1 + x_bar_total)' ...
       + n * (x_bar_int_l2_l2 - x_bar_f1_l2 - x_bar_f2_l2 + x_bar_total) * (x_bar_int_l2_l2 - x_bar_f1_l2 - x_bar_f2_l2 + x_bar_total)' ...
       + n * (x_bar_int_l2_l3 - x_bar_f1_l2 - x_bar_f2_l3 + x_bar_total) * (x_bar_int_l2_l3 - x_bar_f1_l2 - x_bar_f2_l3 + x_bar_total)';
% within
SSW = SST - SSB1 - SSB2 - SSINT;

% Test for systematic interaction between location and variety
% a likelihood ratio test
p = 3;
alpha = 0.01;
LAMBDA = det(SSW) / det(SSW + SSINT);
test_statistics = -(g*b*(n-1) - (p+1-(g-1)*(b-1))/2)*log(LAMBDA);
critical_value = chi2inv(1-alpha, (g-1)*(b-1)*p);
p_value = 1-chi2cdf(test_statistics, (g-1)*(b-1)*p);
disp("-")
disp("Interaction is significant?")
disp("p_value is " + string(p_value))
disp("We conclude that there is no interaction")
disp("")

% Test for systematic factor 1 effect  (location)
% a likelihood ratio test
LAMBDA = det(SSW) / det(SSW + SSB1);
test_statistics = -(g*b*(n-1) - (p+1-(g-1))/2)*log(LAMBDA);
critical_value = chi2inv(1-alpha, (g-1)*p);
p_value = 1-chi2cdf(test_statistics, (g-1)*p);
disp("-")
disp("Factor 1 is significant?")
disp("p_value is " + string(p_value))
disp("We conclude that there is not a factor1 effect")
disp("-")

% Test for systematic factor 2 effect  (variety)
% a likelihood ratio test
LAMBDA = det(SSW) / det(SSW + SSB2);
test_statistics = -(g*b*(n-1) - (p+1-(b-1))/2)*log(LAMBDA);
critical_value = chi2inv(1-alpha, (b-1)*p);
p_value = 1-chi2cdf(test_statistics, (b-1)*p);
disp("Factor 2 is significant?")
disp("p_value is " + string(p_value))
disp("We conclude that there is a factor2 effect, i.e., we reject the null hypothesis")

disp("The location has no effect, but the variety has.")



