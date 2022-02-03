% Load the data
clc; clear; format short

load('dataset_problem_7_1.mat');

% append ones
Z = [ones(62, 1) Z];

% Mean response variables
B = inv(Z'*Z)*Z'*Y;

% Estimated response
Y_hat = Z * B;

% if the covariance is identity, then we can assume ULMR
n = 62;
m = 4;
r = 4;
S = 1/(n-(r+1)) * (Y - Y_hat)' * (Y - Y_hat);

% Perform an uncorrected (not mean-centered) decomposition of (m x m) SS
% matrices
SST = (Y' * Y);
SSR = (Y_hat' * Y_hat);
SSE = SST - SSR;

Y_bar = mean(Y);
SST = (Y - Y_bar)' * (Y - Y_bar);
SSR = (Y_hat - Y_bar)' * (Y_hat - Y_bar);
R_2 = 1 - det(SSE) / det(SST)
R_2_adjusted = 1 - (det(SSE)/(n-(r+1)*m))/(det(SST)/(n-m))

% SST = n-m
% SSR = (r+1)*(m-1)
% SSE = n-(r+1)*m

% Next, for inference for the regression coefficients, indicate the ((r+1) x 1) vector of 
% estimated regression coefficients for Y1
S_Beta1 = S(1, 1) * inv(Z'*Z)

% Simultaneous confidence intervals
alpha = 0.05;
B_conf = B(:, 1) + [-sqrt(diag(S_Beta1)) * sqrt((r+1) * finv(1-alpha, r+1, n-(r+1))) ...
                     sqrt(diag(S_Beta1)) * sqrt((r+1) * finv(1-alpha, r+1, n-(r+1)))]

% Bonferroni confidence intervals
p = r+1;
B_conf = B(:, 1) + [-tinv(1-alpha/(2*p), n-(r+1)) * sqrt(diag(S_Beta1)) ...
                     tinv(1-alpha/(2*p), n-(r+1)) * sqrt(diag(S_Beta1))]


% Test for z0 is neq to zero, but the rest are
% Extra Sum of Squares ;
q = 0;
Z_0 = Z(:, 1);
B_0 = pinv(Z_0) * Y;
SSE_H0 = (Y - Z_0*B_0)' * (Y - Z_0*B_0);

% Wilks Lambda
Lambda = (det(SSE/(n-(r+1))) / (det(SSE_H0/(n-(q+1)))));

% Test-statistics
test_statistics = -(n-(r+1)-(m-r+q+1)/2)*log(Lambda)

% P-value
P = 1-chi2cdf(test_statistics, m*(r-q))

% Perform an Extra Sum of Squares based test for significant predictor
% variable z4
q = 3;
Z_4 = Z(:, 1:4);
B_4 = pinv(Z_4) * Y;
SSE_H4 = (Y - Z_4*B_4)' * (Y - Z_4*B_4);

% Wilks Lambda
Lambda = (det(SSE/(n-(r+1))) / (det(SSE_H4/(n-(q+1)))));

% Test-statistics
test_statistics = -(n-(r+1)-(m-r+q+1)/2)*log(Lambda)

% P-value
P = 1-chi2cdf(test_statistics, m*(r-q))

[beta, Sigma, E, CovB, logL]  = mvregress(Z, Y);


