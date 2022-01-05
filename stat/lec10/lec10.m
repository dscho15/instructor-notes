% Load data
clc; clear; format long

load('dataset_problem_11_1.mat')


%%
clc;
X_1 = X(X(:, 1) == 1, 2:3);
X_2 = X(X(:, 1) == 2, 2:3);
X_3 = X(X(:, 1) == 3, 2:3);
n1 = 50;
n2 = 50;
n3 = 50;

%% Descriptive stat

S_1 = cov(X_1);
S_2 = cov(X_2);
S_3 = cov(X_3);
x_1_bar = mean(X_1);
x_2_bar = mean(X_2);
x_3_bar = mean(X_3);

%% Scatter plot

scatter(X_1(:, 1), X_1(:, 2))
hold on;
scatter(X_2(:, 1), X_2(:, 2))
hold on;
scatter(X_3(:, 1), X_3(:, 2))
legend('1', '2', '3')

%% Model Check

figure(2)
subplot(2, 3, 1)
p = 2
d = diag((X_1 - x_1_bar) / S_1 * (X_1 - x_1_bar)')
y = chi2rnd(p, 10000, 1);
qqplot(d, y)
pbaspect([1 1 1])

subplot(2, 3, 2)
p = 2
d = diag((X_2 - x_2_bar) / S_2 * (X_2 - x_2_bar)')
y = chi2rnd(p, 10000, 1);
figure(2)
qqplot(d, y)
pbaspect([1 1 1])

subplot(2, 3, 3)
p = 2
d = diag((X_3 - x_3_bar) / S_2 * (X_3 - x_3_bar)')
y = chi2rnd(p, 10000, 1);
figure(2)
qqplot(d, y)
pbaspect([1 1 1])

figure(2)
subplot(2, 3, 4)
plotmatrix(X_1)
pbaspect([1 1 1])

subplot(2, 3, 5)
plotmatrix(X_2)
pbaspect([1 1 1])

subplot(2, 3, 6)
plotmatrix(X_3)
pbaspect([1 1 1])

%% Bartlett
p = 2
alpha = 0.05

% the correction factor
C = 1 - (2*p^2 + 3*p-1)/(6*(p+1)) * (1/(n1-1) + 1/(n2-1) + 1/(n3-1) - (1/(n1+n2+n3-3)));

% pooled variance
S_pool = ((n1-1)*S_1 + (n2-1)*S_2 + (n3-1)*S_3) / (n1+n2+n3-3);

% test-statistics for the bartlett test
test_statistics = C*((n1+n2+n3-3)*log(det(S_pool)) - (n1-1)*log(det(S_1)) - (n2-1)*log(det(S_2)) - (n3-1)*log(det(S_3)));

% critical value
critical_value = chi2inv(1-alpha, p*(p+1));

% p-value
p_value = 1-chi2cdf(test_statistics, p*(p+1))

% We reject the null hypothesis, i.e., we have to resolve to QDA

%% However, we still use Manova1 to compare for equal means

[d, p, stats]  = manova1(X(:, 2:3), X(:, 1))

% we reject the null hypothesis of equal means at dimension 2
% however, we still utilize the pooled variance here xD

%%
% cf1, cf2, cf3
cf1 = @(x,y) (-0.5 * [x; y]' * inv(S_1) * [x; y] + x_1_bar * inv(S_1) * [x; y] - 0.5 * log(det(S_1)) - 0.5 * x_1_bar * inv(S_1) * x_1_bar');
cf2 = @(x,y) (-0.5 * [x; y]' * inv(S_2) * [x; y] + x_2_bar * inv(S_2) * [x; y] - 0.5 * log(det(S_2)) - 0.5 * x_2_bar * inv(S_2) * x_2_bar');
cf3 = @(x,y) (-0.5 * [x; y]' * inv(S_3) * [x; y] + x_3_bar * inv(S_3) * [x; y] - 0.5 * log(det(S_3)) - 0.5 * x_3_bar * inv(S_3) * x_3_bar');

figure(3)
decisions = zeros(150, 1);
[Xx,Yy] = meshgrid(1.5:0.05:4.5, 0:0.05:3.0);
colors = ["blue", "green", "red"]

% the descisions of the data
for i = 1:size(X, 1)
    d1 = cf1(X(i, 2), X(i, 3));
    d2 = cf2(X(i, 2), X(i, 3));
    d3 = cf3(X(i, 2), X(i, 3));
    [~, I] = max([d1, d2, d3]);
    decisions(i) = I;
    scatter(X(i, 2), X(i, 3), 20, colors(X(i, 1)));
    hold on;
end

% to plot the decision regions
for i = 1:size(Xx,1)
    for j = 1:size(Xx,2)
        d1 = cf1(Xx(i, j), Yy(i, j));
        d2 = cf2(Xx(i, j), Yy(i, j));
        d3 = cf3(Xx(i, j), Yy(i, j));
        [~, I] = max([d1, d2, d3]);
        scatter(Xx(i, j), Yy(i,j), 1, colors(I));
    end
end

%%
clc;
% calculate the APER
APER = sum(decisions == X(:, 1))/150

% calculate the confusion matrix
cf11 = sum(decisions(1:50) == 1);
cf22 = sum(decisions(51:100) == 2);
cf33 = sum(decisions(101:end) == 3);

cf_mat = diag([cf11 cf22 cf33]);
cf_mat(2, 3) = 4;
cf_mat(3, 2) = 3;
cf_mat

% calculate the AER
for i = 1:size(X_new, 1)
    X_new = X
    X_new(i, :) = []

    X_1 = X_new(X_new(:, 1) == 1, 2:3);
    X_2 = X_new(X_new(:, 1) == 2, 2:3);
    X_3 = X_new(X_new(:, 1) == 3, 2:3);

    S_1 = cov(X_1);
    S_2 = cov(X_2);
    S_3 = cov(X_3);

    x_1_bar = mean(X_1);
    x_2_bar = mean(X_2);
    x_3_bar = mean(X_3);

    cf1 = @(x,y) (-0.5 * [x; y]' * inv(S_1) * [x; y] + x_1_bar * inv(S_1) * [x; y] - 0.5 * log(det(S_1)) - 0.5 * x_1_bar * inv(S_1) * x_1_bar');
    cf2 = @(x,y) (-0.5 * [x; y]' * inv(S_2) * [x; y] + x_2_bar * inv(S_2) * [x; y] - 0.5 * log(det(S_2)) - 0.5 * x_2_bar * inv(S_2) * x_2_bar');
    cf3 = @(x,y) (-0.5 * [x; y]' * inv(S_3) * [x; y] + x_3_bar * inv(S_3) * [x; y] - 0.5 * log(det(S_3)) - 0.5 * x_3_bar * inv(S_3) * x_3_bar');

    d1 = cf1(X_new(i, 2), X_new(i, 3)); 
    d2 = cf2(X_new(i, 2), X_new(i, 3));
    d3 = cf3(X_new(i, 2), X_new(i, 3));

    [~, I] = max([d1, d2, d3]);
    
    decisions(i) = I;
    
    hold on;
end

