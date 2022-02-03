% 
clc; clear;
format short

% The file %dataset_problem_5_1.mat  contains three 30x4 data matrices containing size 
% measurements for Egyptian skull samples. Each of the three samples is from a specific 
% time period in history (4000 B.C., 3300 B.C. and 1850 B.C., respectively) and contains 30 
% skull observations for that period

load('dataset_problem_5_1.mat')

alpha = 0.05;
p = 4;
n1 = 30;
n2 = 30;
n3 = 30;

S_time1 = cov(X_time1);
S_time2 = cov(X_time2);
S_time3 = cov(X_time3);
S_pool = ((n1-1)*S_time1 + (n2-1)*S_time2 + (n3-1)*S_time3) / (n1+n2+n3-3);

x_bar_time1 = mean(X_time1);
x_bar_time2 = mean(X_time2);
x_bar_time3 = mean(X_time3);

C = 1 - (2*p^2 + 3*p-1)/(6*(p+1)) * (1/(n1-1) + 1/(n2-1) + 1/(n3-1) - 1/(n1+n2+n3-3));
test_statistics = C*((n1+n2+n3-3)*log(det(S_pool))-(n1-1)*log(det(S_time1))-(n2-1)*log(det(S_time2)) - (n3-1)*log(det(S_time3)));
critical_value = chi2inv(1-alpha, p*(p+1));

disp("Reject [y/n]: " + string(test_statistics > critical_value))
disp("p-value: " + string(1 - chi2cdf(test_statistics, p*(p+1))))
disp("We have passed homoscedasticity test")

%%
figure(1)
subplot(2,3,1)
boxplot(X_time1)
pbaspect([1 1 1])
subplot(2,3,2)
boxplot(X_time2)
pbaspect([1 1 1])
subplot(2,3,3)
boxplot(X_time3)
pbaspect([1 1 1])
rng('default') % For reproducibility
y = chi2rnd(p, 10000, 1);
subplot(2, 3, 4);
d = diag((X_time1 - x_bar_time1) * inv(S_pool) * (X_time1 - x_bar_time1)');
qqplot(d, y);
pbaspect([1 1 1]);
subplot(2, 3, 5);
d = diag((X_time2 - x_bar_time2) * inv(S_pool) * (X_time2 - x_bar_time2)');
qqplot(d, y);
pbaspect([1 1 1]);
subplot(2, 3, 6);
d = diag((X_time3 - x_bar_time3) * inv(S_pool) * (X_time3 - x_bar_time3)');
qqplot(d, y);
pbaspect([1 1 1]);
% seems like it is "okay" normally distributed

%% Test for equality of skull size mean vectors from the three time periods (MANOVA1)
% We want to use one factor Manova with three factor levels (3 populations)
g = 3;

% compose
x_bar_total = mean([X_time1; X_time2; X_time3]);
SST = cov([X_time1; X_time2; X_time3]) * (n1+n2+n3-1)
SSB = (n1*(x_bar_time1 - x_bar_total)' * (x_bar_time1 - x_bar_total) +  ... 
       n2*(x_bar_time2 - x_bar_total)' * (x_bar_time2 - x_bar_total) +  ... 
       n3*(x_bar_time3 - x_bar_total)' * (x_bar_time3 - x_bar_total));
SSW = SST - SSB;

% if they are never needed
ST = SST/(n1+n2+n3-1);
SB = SSB/(g-1);

% Wilks_Lambda
LAMBDA = det(SSW)/det(SST)
test_statistics = (n1+n2+n3-p-2)/p * (1-sqrt(LAMBDA))/(sqrt(LAMBDA))
critical_value = finv(1-alpha, 2*p, 2*(n1+n2+n3-p-2))

disp("Reject [y/n]: " + string(test_statistics > critical_value))
disp("p-value: " + string(1 - fcdf(test_statistics, 2*p, 2*(n1+n2+n3-p-2))))

% Large sample approximation
test_statistics = -(n1+n2+n3-1-(p+g)/2) * log(LAMBDA);
critical_value = chi2inv(1-alpha, p*(g-1));

disp("Reject [y/n]: " + string(test_statistics > critical_value))
disp("p-value: " + string(1 - chi2cdf(test_statistics, p*(g-1))))

%
X_time = [X_time1; X_time2; X_time3];
time_index = [1*ones(n1,1); 2*ones(n2,1); 3*ones(n3,1)];
[dim,p_values,stats] = manova1(X_time,time_index);
SSB = stats.B;
SSW = stats.W;
df_B = stats.dfB;
df_W = stats.dfW;
chisq = stats.chisq(1);
p_value = p_values(1);
disp("Build in")
disp("p-value: " + string(p_value))

% make the bonferroni plots later







