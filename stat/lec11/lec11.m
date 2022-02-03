% 
load('dataset_problem_12_3.mat')
load('university_names_problem_12_3.mat')

% Since the 6 variables have different units and different dynamic range, standardize
% the raw data matrix, X, to a matrix, Z, who s columns (variables) are all zero-mean and unit-variance

% Single is clostest
% Complete is furtherst away
% Centroid is cluster
% Average is all possible between two clusters

Z = normalize(X)
Z_sgl = linkage(Z, 'single');
Z_avg = linkage(Z, 'average');
Z_com = linkage(Z, 'complete');
Z_cen = linkage(Z, 'centroid');

subplot(2, 2, 1)
dendrogram(Z_sgl)
title('Single')

subplot(2, 2, 2)
cutoff = 4
CL = cluster(Z_avg, 'Cutoff', cutoff, 'Criterion', 'distance');
dendrogram(Z_avg, 'ColorThreshold', cutoff);
title('Average')

subplot(2, 2, 3)
cutoff = 6
CL = cluster(Z_com, 'Cutoff', cutoff, 'Criterion', 'distance');
dendrogram(Z_com, 'ColorThreshold', cutoff);
title('Complete')

subplot(2, 2, 4)
cutoff = 3.5
CL = cluster(Z_cen, 'Cutoff', cutoff, 'Criterion', 'distance');
dendrogram(Z_cen, 'ColorThreshold', cutoff);
title('Centroid')

%%
% K-means clustering (using Z): U
% Use the K-means algorithm to group the universities into K = 3, 4, 5 and 6 clusters, respectively. 
% For each value of K, allow several different initial clusterings ("replicates"), 
% thus allowing the algorithm to output the final clustering with lowest Within-Sum-of-Squares. 
% Indicate which universities are in each of the clusters for each value of K, if possible, comment on the groupings, and compare with above.

Z = normalize(X)
[clusters_3, C_3, ssum_3] = kmeans(Z, 3, 'Distance', 'sqeuclidean', 'Replicates', 10, 'Start', 'sample')
[clusters_4, C_4, ssum_4] = kmeans(Z, 4, 'Distance', 'sqeuclidean', 'Replicates', 10, 'Start', 'sample')
[clusters_5, C_5, ssum_5] = kmeans(Z, 5, 'Distance', 'sqeuclidean', 'Replicates', 10, 'Start', 'sample')
[clusters_6, C_6, ssum_6] = kmeans(Z, 6, 'Distance', 'sqeuclidean', 'Replicates', 10, 'Start', 'sample')
stem([3, 4, 5, 6], [mean(ssum_3), mean(ssum_4), mean(ssum_5), mean(ssum_6)])

%%
% Problem 12.4 (MATLAB)Gaussian Mixture Model clustering of city crime rates dat
clear;
load('dataset_problem_12_4.mat')
load('city_names_problem_12_4.mat')
scatter(X(:, 1), X(:, 2))

%%
gm2 = fitgmdist(X,2)

scatter(X(:,1), X(:,2), 20,'.') % Scatter plot with points of size 10
hold on
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
fcontour(gmPDF,[0 60])

%% 
figure(3)
idx = cluster(gm, X);
c1 = (idx == 1)
hold on;
scatter(X(c1, 1), X(c1, 2), 30)
textscatter(X(c1, 1), X(c1, 2), names(c1))

c2 = (idx == 2)
scatter(X(c2, 1), X(c2, 2), 30)
textscatter(X(c2, 1), X(c2, 2), names(c2))

% mvnpdf(gm.mu(1, :), Sigma(1))
% mvnpdf(gm.mu(2, :), Sigma(2))

%% 
gm3 = fitgmdist(X, 3)
scatter(X(:,1), X(:,2), 100,'x')
hold on;
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm3,[x0 y0]),x,y);
fcontour(gmPDF,[0 60])

[ gm2.AIC gm2.BIC ]
[ gm3.AIC gm3.BIC ]

