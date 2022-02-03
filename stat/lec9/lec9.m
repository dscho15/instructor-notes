% classification10 - two-populations linear and quadratic discriminant analysis
clc; clear; format short
close all;
load('dataset_problem_10_1.mat')

X_m = X(X(:, 1) == 1, :)
X_f = X(X(:, 1) == 2, :)

% Scatter plots
figure(1)
scatter(X_m(:, 2), X_m(:, 3))
hold on;
scatter(X_f(:, 2), X_f(:, 3))
legend('male', 'female')
grid on;
xlabel('tail length[mm]')
ylabel('snout to vent length[mm]')

% Calculate descriptive statistics for both populations, i.e. sample mean vectors and
% sample covariance matrices
x_m_bar = mean(X_m(:, 2:end));
S_m = cov(X_m(:, 2:end));
R_m = corrcoef(X_m(:,2:end));

x_f_bar = mean(X_f(:,2:end));
S_f = cov(X_f(:,2:end));
R_m = corrcoef(X_m(:,2:end));

% Model check
figure(2)
subplot(2, 2, 1)
p = 2
d = diag((X_m(:, 2:end) - x_m_bar) / S_m * (X_m(:, 2:end) - x_m_bar)')
y = chi2rnd(p, 10000, 1);
qqplot(d, y)
pbaspect([1 1 1])
subplot(2, 2, 2)
p = 2
d = diag((X_f(:, 2:end) - x_f_bar) / S_f * (X_f(:, 2:end) - x_f_bar)')
y = chi2rnd(p, 10000, 1);
figure(2)
qqplot(d, y)
pbaspect([1 1 1])
figure(2)
subplot(2, 2, 3)
plotmatrix(X_m(:, 2:end))
pbaspect([1 1 1])
hold on;
subplot(2, 2, 4)
plotmatrix(X_f(:, 2:end))
pbaspect([1 1 1])

%% Bartlett test
clc; close all;
n1 = 37;
n2 = 29;
p = 2;
alpha = 0.05;
% the correction factor
C = 1 - (2*p^2 + 3*p-1)/(6*(p+1)) * (1/(n1-1) + 1/(n2-1) - (1/(n1+n2-2)));
% pooled variance
S_pool = ((n1-1)*S_m + (n2-1)*S_f) / (n1+n2-2);
% test-statistics for the bartlett test
test_statistics = C*((n1+n2-2)*log(det(S_pool)) - (n1-1)*log(det(S_m)) - (n2-1)*log(det(S_f)));
% critical value
critical_value = chi2inv(1-alpha, p*(p+1)/2);
% p-value
p_value = 1-chi2cdf(test_statistics, p*(p+1)/2)

% The low p-value for this test should indicate that the covariance matrices most
% probably are not equal, and thus QDA should be the most appropriate classification
% approach

%% LDA
clc;
S_pool = ((n1-1)*S_m + (n2-1)*S_f) / (n1+n2-2);
% we want to test for equal mean vectors, if they are equal, then it
% doesn't make sense to move on, as we assume equal covariance

test_statistics = (x_m_bar-x_f_bar) / (S_pool * (1/n1 + 1/n2)) * (x_m_bar-x_f_bar)'
critical_value = (n1+n2-2)*p/(n1+n2-p-1) * finv(1-alpha, p, n1+n2-p-1)
p_value = 1 - fcdf((n1+n2-p-1)/((n1+n2-2)*p) * test_statistics , p, n1+n2-p-1)
% we reject the null hypothesis of equal means

clc;
% calculate distriminant functions
figure(1)
scatter(X_m(:, 2), X_m(:, 3))
hold on;
scatter(X_f(:, 2), X_f(:, 3))
legend('male', 'female')
grid on;
fimplicit(@(x,y) (x_m_bar - x_f_bar) / S_pool * ([x; y] - 1/2 * (x_m_bar + x_f_bar)'), [100 220 350 700], 'black', 'LineWidth', 2);
xlabel('tail length[mm]')
ylabel('snout to vent length[mm]')

% Performance evalulation of LDA classification
classification = @(x,y)(x_m_bar - x_f_bar) / S_pool * ([x; y] - 1/2 * (x_m_bar + x_f_bar)');
male = zeros(66, 1);
female = zeros(66, 1);
for i = 1:66
    male(i, classification(X(i, 2), X(i, 3)) > 0) = 1;
    female(i, classification(X(i, 2), X(i, 3)) < 0) = 1;
end
correct_male = sum(male(1:37));
false_male = 37-correct_male;
correct_female = sum(female(38:66));
false_female = 29-correct_female;
confusion_matrix = [[correct_male, false_male];
                    [false_female, correct_female]]

%% AER
classifications = zeros(66, 1)
for i = 1:66

    X_new = X;
    X_new(i, :) = [];

    X_m = X_new(X_new(:, 1) == 1, :);
    X_f = X_new(X_new(:, 1) == 2, :);
    
    n1 = size(X_m, 1)
    n2 = size(X_f, 1)

    x_m_bar = mean(X_m(:, 2:end));
    S_m = cov(X_m(:, 2:end));
    x_f_bar = mean(X_f(:,2:end));
    S_f = cov(X_f(:,2:end));
    S_pool = ((n1-1)*S_m + (n2-1)*S_f) / (n1+n2-2);

    classification = @(x,y)(x_m_bar - x_f_bar) / S_pool * ([x; y] - 1/2 * (x_m_bar + x_f_bar)');
    
    classifications(i, classification(X(i, 2), X(i, 3)) > 0) = 1;
end
classifications(classifications(:) == 0) = 2
actual_error_rate = sum(X(:, 1) == classifications)

%% Classifications of new observations using LDA
clc;
X_m = X(X(:, 1) == 1, :);
X_f = X(X(:, 1) == 2, :);
x_m_bar = mean(X_m(:, 2:end));
S_m = cov(X_m(:, 2:end));
x_f_bar = mean(X_f(:,2:end));
S_f = cov(X_f(:,2:end));
classification = @(x,y)(x_m_bar - x_f_bar) / S_pool * ([x; y] - 1/2 * (x_m_bar + x_f_bar)');
classification(170, 550) % female
classification(200, 600) % male
classification(120, 380) % male

% Or we can classify by using posterior probabilities, i.e. we should have
% larger than 90%
P1 = mvnpdf([170; 550], x_m_bar', S_m)
P2 = mvnpdf([170; 550], x_f_bar', S_f)
P1 / (P1 + P2)
P2 / (P1 + P2)

%% QDA
clc; clear;
load('dataset_problem_10_1.mat')

X_m = X(X(:, 1) == 1, :);
X_f = X(X(:, 1) == 2, :);

x_m_bar = mean(X_m(:, 2:end));
S_m = cov(X_m(:, 2:end));

x_f_bar = mean(X_f(:,2:end));
S_f = cov(X_f(:,2:end));
n1 = 37
n2 = 29
p = 2

% Hotelling test statittics

test_statistics = (x_m_bar - x_f_bar) / (S_m/n1 + S_f/n2) * (x_m_bar - x_f_bar)'
p_value = 1 - chi2cdf(test_statistics, p)

%%

% calculate distriminant functions
figure(1)
scatter(X_m(:, 2), X_m(:, 3))
hold on;
scatter(X_f(:, 2), X_f(:, 3))
legend('male', 'female')
grid on;
K = 1/2 * log(det(S_m) / det(S_f)) + 1/2 * (x_m_bar / S_m * x_m_bar' - x_f_bar / S_f * x_f_bar');
fimplicit(@(x,y) (-1/2 * [x, y] * (inv(S_m) - inv(S_f)) * [x; y] + (x_m_bar / S_m ...
                        - x_f_bar / S_f) * [x; y] - K), [100 220 350 700], 'black', 'LineWidth', 2);
xlabel('tail length[mm]')
ylabel('snout to vent length[mm]')

%% Performance evalulation of CDA classification % Apparent
classification = @(x,y)(-1/2 * [x, y] * (inv(S_m) - inv(S_f)) * [x; y] + (x_m_bar / S_m - x_f_bar / S_f) * [x; y] - K);
male = zeros(66, 1);
female = zeros(66, 1);
for i = 1:66
    male(i, classification(X(i, 2), X(i, 3)) > 0) = 1;
    female(i, classification(X(i, 2), X(i, 3)) < 0) = 1;
end
correct_male = sum(male(1:37));
false_male = 37-correct_male;
correct_female = sum(female(38:66));
false_female = 29-correct_female;
confusion_matrix = [[correct_male, false_male];
                    [false_female, correct_female]]

%% AER % Actual error rate
classifications = zeros(66, 1);
for i = 1:66

    X_new = X;
    X_new(i, :) = [];

    X_m = X_new(X_new(:, 1) == 1, :);
    X_f = X_new(X_new(:, 1) == 2, :);
    
    n1 = size(X_m, 1);
    n2 = size(X_f, 1);

    x_m_bar = mean(X_m(:, 2:end));
    S_m = cov(X_m(:, 2:end));
    x_f_bar = mean(X_f(:,2:end));
    S_f = cov(X_f(:,2:end));
    S_pool = ((n1-1)*S_m + (n2-1)*S_f) / (n1+n2-2);

    classification = @(x,y)(x_m_bar - x_f_bar) / S_pool * ([x; y] - 1/2 * (x_m_bar + x_f_bar)');
    
    classifications(i, classification(X(i, 2), X(i, 3)) > 0) = 1;
end
classifications(classifications(:) == 0) = 2;
actual_error_rate = sum(X(:, 1) == classifications)

%% Predictions

classification(170, 550) % female
classification(200, 600) % male
classification(120, 380) % male
