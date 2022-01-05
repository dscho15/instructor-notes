clc; clear; format short

load("dataset_problem_8_1.mat")

X = X_time

x_bar = mean(X)
S = cov(X)
R = corrcoef(X)
figure(1)
boxplot(X)
figure(2)
plotmatrix(X)

%% Perform a PCA on R using svd or eig
[U, S, V] = svd(R)

% U * S * V'
variation = diag(S)/sum(diag(S))
pareto(variation)

%% Indicate the loadings for PC1 (the PC with highest variance explained/highest 
% eigenvalue) on Z-variables (standardized variables) and X-variables respectively
% The difference is to do PCA on R and S
[U, S, V] = svd(R);
% The maximum loading of R
V = (-1)*V;
V(:,1)

% The loading for X, this is complete garbo unsued
% [U_X, S_X, V_X] = svd(X);
% V_X = (-1)*V_X;

%% Calculate the PC1 scores for the 54 countries 
% The scores for the 54 countries
Z = normalize(X);
V_1 = V(:,1)
scores_1 = Z*V(:,1)

%% Indicate the loadings for PC2 (the PC with highest variance explained/highest 
% eigenvalue) on Z-variables (standardized variables) and X-variables respectively
% The difference is to do PCA on R and S
[U, S, V] = svd(R);
% The maximum loading of R
V = (-1)*V;
V_2 = V(:,2)

%% Calculate the PC2 scores for the 54 countries 
% The scores for the 54 countries
Z = normalize(X);
scores_2 = Z*V(:,2)

%% Make some scatter plots and biplots
str = string(1:54);
figure(3)
textscatter(scores_1, scores_2, str);
figure(4)
vbls = {'100m','200m','400m','800m','1500m', '5000m', '10000m', '42192m'}; % Labels for the variables
biplot([V_1 V_2], 'Scores', [scores_1, scores_2], 'Varlabels', vbls)
grid on

%% Calculate the PC3 scores for the 54 countries 
% The scores for the 54 countries
Z = normalize(X);
scores_3 = Z*V(:,3)
V_3 = V(:, 3)

%% Make some scatter plots
figure(5)
vbls = {'100m','200m','400m','800m','1500m', '5000m', '10000m', '42192m'}; % Labels for the variables
biplot([V_1 V_2 V_3], 'Scores', [scores_1, scores_2, scores_3], 'Varlabels', vbls)
grid on

%% Calculate the sample correlation between PC1(vector) and the original variables 
% and likewise for PC2 and interpret the values
correlation_vector_1 = sqrt(S(1, 1)) * V_1
correlation_vector_2 = sqrt(S(2, 2)) * V_2

% Large Sample 
%
% 



