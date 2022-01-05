% 
clc; clear;
format short

% The file %dataset_problem_6_1.mat  contains data for the measured characteristics of 76 young bulls 
% (less than two years old) sold at an auction.  The univariate (scalar) response variable, Y, is the bull%s height in [inches] 
% at shoulder at the time of the sale.  The two predictor variables, collected in the matrix, Z, is the bull s yearling (one year old) 
% height in [inches] at shoulder and the bull’s fat free body weight in [pounds], respectively

load('dataset_problem_6_1.mat')
r = 2
n = 76
Z = [ones(76,1) Z]

% Perform least-squares
B = pinv(Z)*Y

% Sigma 
sigma_2 = ((Y-Z*B)' * (Y-Z*B))/(n-(r+1))

% SS-decomposition
n = 76

% total sum of square
SST = var(Y) * (n-1)

% residual (error) sum of squares
SSE = var(Y-Z*B) * (n-1)

% regression sum of squares
SSR = SST - SSE

% coeeffiicient of determination
R_2 = 1 - SSE / SST
R_2_adjusted = (SSR / (n-(r+1))) / (SST / (n-1))

% 
B_cov = sigma_2 * inv(Z'*Z);

% simulatenous
alpha = 0.05;
B_conf = B + [-sqrt(diag(B_cov)) * sqrt((r+1) * finv(1-alpha, r+1, n-(r+1))) ...
               sqrt(diag(B_cov)) * sqrt((r+1) * finv(1-alpha, r+1, n-(r+1)))]

% bonferroni confidence intervals
p = r+1;
B_conf = B + [-tinv(1-alpha/(2*p), n-(r+1)) * sqrt(diag(B_cov)) ...
               tinv(1-alpha/(2*p), n-(r+1)) * sqrt(diag(B_cov)) ];

%  Perform an Extra Sum of Squares based test for significant dependency of Y on 
%  (z1,z2) (significant regression model), that is test the hypothesis
q = 0;
Z_1 = Z(:, 1);
B_1 = pinv(Z_1)*Y;
SS_EXTRA = (Y-Z_1*B_1)'*(Y-Z_1*B_1) - (Y-Z*B)'*(Y-Z*B);
test_statistics = (SS_EXTRA / (r - q)) / sigma_2; 
critical_value = finv(1-alpha, r-q, n-(r+1));
p_value = 1-fcdf(test_statistics, r-q, n-(r+1));
disp("model check")
disp("p_value: " + string(p_value))

% fischer test of equal variance

% Perform an ANOVA table based test for the same hypothesis and demonstrate that this approach is equivalent to the Extra Sum of Squares approach     
% Perform marginal tests  H଴,୧:  β୧ൌ0,iൌ0,1,2  and compare the results with above    
% Calculate a 95% prediction interval for E[Y(z0)], the mean of new responses at z=z0, for z0 = [z1 z2] = [50.5  970]     
% Calculate a 95% prediction interval for Y(z0), a single new response at z=z0

figure(1)
subplot(2, 3, 1)
stem((Y-Z*B))
subplot(2, 3, 2)
stem(sort(Y-Z*B))
subplot(2, 3, 3)
stem(sort(abs(Y-Z*B)))
subplot(2, 3, 4)
hist((Y-Z*B))
subplot(2, 3, 5)
qqplot((Y-Z*B))
subplot(2, 3, 6)
stem(Y, (Y-Z*B))
xlabel('y')
ylabel('\epsilon')
grid

% The model seems to normal distributed and it seems to be a valid model
% as we can regress about 80% of the model

