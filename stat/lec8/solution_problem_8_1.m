clear;close all;clc;
format compact

load('dataset_problem_8_1.mat');
% NOTE: Country ID numbers is found in the file 'countriesID_problem_8_1.dat'

disp('---------------------------------------------------------------------------------')
disp('Descriptive statistics estimators for male time-data')
disp('---------------------------------------------------------------------------------')
[n p] = size(X_time)
Xmean_time = mean(X_time,1)
S_time = cov(X_time)
R_time = corr(X_time)

figure('Name','Box plot for original male time-data')
boxplot(X_time,'labels',{'100 m [sec]','200 m [sec]','400 m [sec]','800 m [min]','1500 m [min]','5000 m [min]','10000 m [min]','Marathon [min]'}),title('Box plot for original male time-data','Fontsize',14)
figure('Name','Scatter plot for original male time-data')
plotmatrix(X_time),title('Scatter plot for original male time-data','Fontsize',14)
xlabel('100 m                  200 m                   400 m                   800 m                   1500 m                   5000 m                  10000 m                  Marathon','Fontsize',14)
ylabel('Marat.     10000     5000       1500        800        400        200        100','Fontsize',14)

disp('---------------------------------------------------------------------------------')
disp('PCA on R for male time-data')
disp('---------------------------------------------------------------------------------')
[E_R_time LAMBDA_R_time E_R_time] = svd(R_time);
E_R_time = sign(E_R_time(1,1))*E_R_time;
LAMBDA_R_time,E_R_time

disp('---------------------------------------------------------------------------------')
disp('Variance explained by principal components')
disp('---------------------------------------------------------------------------------')
PC_R_time_percent_explained = 100*diag(LAMBDA_R_time)/sum(diag(LAMBDA_R_time))
figure('Name','Pareto plot for male time-data')
bar(1:p,PC_R_time_percent_explained),xlabel('PC'),ylabel('% variance of PC''s'), hold on
title('Pareto plot and scree plot for male time-data','Fontsize',14)
accumulated = zeros(1,p);
accumulated(1) = PC_R_time_percent_explained(1);
for i = 2:p
    accumulated(i) = sum(PC_R_time_percent_explained(1:i));
end
plot(1:p,accumulated),grid
hold off

V = diag(diag(S_time));

disp('---------------------------------------------------------------------------------')
disp('Loadings for PC 1')
disp('---------------------------------------------------------------------------------')
PC_1_R_time_loadings_Z_variables = E_R_time(:,1)'
PC_1_R_time_loadings_X_variables = ((V^(-1/2))*E_R_time(:,1))'

disp('---------------------------------------------------------------------------------')
disp('Ranking of observations by PC 1')
disp('---------------------------------------------------------------------------------')
Z_time = zscore(X_time);
PC_1_R_time_scores_Z_variables = Z_time*PC_1_R_time_loadings_Z_variables';
[time_rank_Z index_Z] = sort(PC_1_R_time_scores_Z_variables);
time_rank_index_Z = (index_Z(1:10))'
PC_1_R_time_scores_X_variables = (X_time - ones(n,1)*mean(X_time))*PC_1_R_time_loadings_X_variables';
[time_rank_X index_X] = sort(PC_1_R_time_scores_X_variables);
time_rank_index_X = (index_X(1:10))'
disp('---------------------------------------------------------------------------------')
disp('Variance explained by PC 1')
disp('---------------------------------------------------------------------------------')
PC_1_R_time_variance = LAMBDA_R_time(1,1)
PC_1_R_time_variance_explained_ratio = PC_1_R_time_variance/trace(LAMBDA_R_time)

disp('---------------------------------------------------------------------------------')
disp('Loadings for PC 2')
disp('---------------------------------------------------------------------------------')
PC_2_R_time_loadings_Z_variables = E_R_time(:,2)'
PC_2_R_time_loadings_X_variables = ((V^(-1/2))*E_R_time(:,2))'

disp('---------------------------------------------------------------------------------')
disp('Ranking of observations by PC 2')
disp('---------------------------------------------------------------------------------')
PC_2_R_time_scores_Z_variables = Z_time*PC_2_R_time_loadings_Z_variables';
[time_rank_Z index_Z] = sort(PC_2_R_time_scores_Z_variables);
time_rank_index_Z = (index_Z(1:10))'
PC_2_R_time_scores_X_variables = (X_time - ones(n,1)*mean(X_time))*PC_2_R_time_loadings_X_variables';
[time_rank_X index_X] = sort(PC_2_R_time_scores_X_variables);
time_rank_index_X = (index_X(1:10))'

disp('---------------------------------------------------------------------------------')
disp('Variance explained by PC 1 and 2')
disp('---------------------------------------------------------------------------------')
PC_2_R_time_variance = LAMBDA_R_time(2,2)
PC_1_2_R_time_variance_explained_ratio = (PC_1_R_time_variance+PC_2_R_time_variance)/trace(LAMBDA_R_time)

figure('Name','Scatterplots for male time-data PC 1 and 2 scores')
plot(PC_1_R_time_scores_Z_variables,PC_2_R_time_scores_Z_variables,'+'),grid,axis equal  %X-scores same result
title('Scatterplots for male time-data PC 1 and 2 scores','Fontsize',14)
xlabel('1st Principal Component scores')
ylabel('2nd Principal Component scores')

figure('Name','Bi-plot for male time-data PC 1 and 2 scores and loadings') 
names = char('100','200','400','800','1500','5000','10000','Marathon');
biplot(E_R_time(:,1:2),'scores',[PC_1_R_time_scores_Z_variables PC_2_R_time_scores_Z_variables],'varlabels',names)
title('Bi-plot for male time-data PC 1 and 2 scores and loadings','Fontsize',14)

figure('Name','Bi-plot for male time-data PC 1, 2 and 3 scores and loadings')
PC_3_R_time_loadings_Z_variables = E_R_time(:,3)';
PC_3_R_time_scores_Z_variables = Z_time*PC_3_R_time_loadings_Z_variables';
biplot(E_R_time(:,1:3),'scores',[PC_1_R_time_scores_Z_variables PC_2_R_time_scores_Z_variables PC_3_R_time_scores_Z_variables],'varlabels',names)
title('Bi-plot for male time-data PC 1, 2 and 3 scores and loadings','Fontsize',14)

disp('---------------------------------------------------------------------------------')
disp('Correlation between PC 1 and 2 and original variables')
disp('---------------------------------------------------------------------------------')
Correlations_PC_1_R_time_versus_Z_or_X_variables = sqrt(LAMBDA_R_time(1,1))*PC_1_R_time_loadings_Z_variables
Correlations_PC_2_R_time_versus_Z_or_X_variables = sqrt(LAMBDA_R_time(2,2))*PC_2_R_time_loadings_Z_variables
figure('Name','Circle diagram for correlation between time-data PC 1 and 2 and original variables') 
t = 0:.01:pi;
x = cos(t);
y = sin(t);
plot(x,y,'k','LineWidth',3),axis square,title('Circle diagram for correlation between time-data PC 1 and 2 and original variables','Fontsize',14),grid,xlim([-1.2 1.2]),ylim([-1.2 1.2])
xlabel('Correlation with PC 1','Fontsize',14),ylabel('Correlation with PC 2','Fontsize',14)
hold on
plot(x,-y,'k','LineWidth',3);
plot(Correlations_PC_1_R_time_versus_Z_or_X_variables, Correlations_PC_2_R_time_versus_Z_or_X_variables,'r+','Linewidth',10)
text(0.93,0.45,'100','Fontsize',12)
text(0.75,0.38,'200','Fontsize',12)
text(0.73,0.28,'400','Fontsize',12)
text(0.77,-0.07,'800','Fontsize',12)
text(0.77,-0.15,'1500','Fontsize',12)
text(0.78,-0.24,'5000','Fontsize',12)
text(0.98,-0.27,'10000','Fontsize',12)
text(0.7,-0.31,'Marath','Fontsize',12)
text(-0.5,0.7,'Since        {\Sigma}_j_=_1_._._p r^2(X_i,PC_j) = 1','Fontsize',14)
text(-0.5,0.55,'it follows   r^2(X_i,PC_1) + r^2(X_i,PC_2) {\leq} 1','Fontsize',14)
hold off

disp('---------------------------------------------------------------------------------')
disp('Approx. simult. Bonferroni CI for variances of time-data PC 1 and 2')
disp('---------------------------------------------------------------------------------')
alpha = 0.1
psim = 2;
alpha_bonf = alpha/psim
CI_approx_lambda_1 = LAMBDA_R_time(1,1)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
CI_approx_lambda_2 = LAMBDA_R_time(2,2)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
disp('---------------------------------------------------------------------------------')

disp('---------------------------------------------------------------------------------')
disp('  ')
%------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------

disp('-------------------------------------------------------------------------------------------------------------------')
disp('-------------------------------------------------------------------------------------------------------------------')
disp('-------------------------------------------------------------------------------------------------------------------')

disp('  ')
disp('---------------------------------------------------------------------------------')
disp('Descriptive statistics estimators for male speed-data')
disp('---------------------------------------------------------------------------------')
X_speed(:,1) = 100./X_time(:,1);
X_speed(:,2) = 200./X_time(:,2);
X_speed(:,3) = 400./X_time(:,3);
X_speed(:,4) = 800./(60*X_time(:,4));
X_speed(:,5) = 1500./(60*X_time(:,5));
X_speed(:,6) = 5000./(60*X_time(:,6));
X_speed(:,7) = 10000./(60*X_time(:,7));
X_speed(:,8) = 42195./(60*X_time(:,8));
Xmean_speed = mean(X_speed,1)
S_speed = cov(X_speed)
R_speed = corr(X_speed)

figure('Name','Box plot for original male speed-data')
boxplot(X_speed,'labels',{'100 m [m/sec]','200 m [m/sec]','400 m [m/sec]','800 m [m/sec]','1500 m [m/sec]','5000 m [m/sec]','10000 m [m/sec]','Marathon [m/sec]'}),title('Box plot for original male speed-data','Fontsize',14)
figure('Name','Scatter plot for original male speed-data')
plotmatrix(X_speed),title('Scatter plot for original male speed-data','Fontsize',14)
xlabel('100 m                  200 m                   400 m                   800 m                   1500 m                   5000 m                  10000 m                  Marathon','Fontsize',14)
ylabel('Marat.     10000     5000       1500        800        400        200        100','Fontsize',14)

disp('---------------------------------------------------------------------------------')
disp('PCA on S for male speed-data')
disp('---------------------------------------------------------------------------------')
[E_S_speed LAMBDA_S_speed E_S_speed] = svd(S_speed);
E_S_speed = sign(E_S_speed(1,1))*E_S_speed;
LAMBDA_S_speed,E_S_speed

disp('---------------------------------------------------------------------------------')
disp('Variance explained by principal components')
disp('---------------------------------------------------------------------------------')
PC_S_speed_percent_explained = 100*diag(LAMBDA_S_speed)/sum(diag(LAMBDA_S_speed))
figure('Name','Pareto plot for male speed-data')
bar(1:p,PC_S_speed_percent_explained),xlabel('PC'),ylabel('% variance of PC''s'), hold on
title('Pareto plot and scree plot for male speed-data','Fontsize',14)
accumulated = zeros(1,p);
accumulated(1) = PC_S_speed_percent_explained(1);
for i = 2:p
    accumulated(i) = sum(PC_S_speed_percent_explained(1:i));
end
plot(1:p,accumulated),grid
hold off

disp('---------------------------------------------------------------------------------')
disp('Loadings for PC 1')
disp('---------------------------------------------------------------------------------')
PC_1_S_speed_loadings_X_variables = (E_S_speed(:,1))'

disp('---------------------------------------------------------------------------------')
disp('Ranking of observations by PC 1')
disp('---------------------------------------------------------------------------------')
PC_1_S_speed_scores_X_variables = (X_speed - ones(n,1)*mean(X_speed))*PC_1_S_speed_loadings_X_variables';
[speed_rank_X index_X] = sort(PC_1_S_speed_scores_X_variables,'descend');
speed_rank_index_X = (index_X(1:10))'

disp('---------------------------------------------------------------------------------')
disp('Variance explained by PC 1')
disp('---------------------------------------------------------------------------------')
PC_1_S_speed_variance = LAMBDA_S_speed(1,1)
PC_1_S_speed_variance_explained_ratio = PC_1_S_speed_variance/trace(LAMBDA_S_speed)

disp('---------------------------------------------------------------------------------')
disp('Loadings for PC 2')
disp('---------------------------------------------------------------------------------')
PC_2_S_speed_loadings_X_variables = (E_S_speed(:,2))'

disp('---------------------------------------------------------------------------------')
disp('Ranking of observations by PC 2')
disp('---------------------------------------------------------------------------------')
PC_2_S_speed_scores_X_variables = (X_speed - ones(n,1)*mean(X_speed))*PC_2_S_speed_loadings_X_variables';
[speed_rank_X index_X] = sort(PC_2_S_speed_scores_X_variables,'descend');
speed_rank_index_X = (index_X(1:10))'

disp('---------------------------------------------------------------------------------')
disp('Variance explained by PC 1 and 2')
disp('---------------------------------------------------------------------------------')
PC_2_S_speed_variance = LAMBDA_S_speed(2,2)
PC_1_2_S_speed_variance_explained_ratio = (PC_1_S_speed_variance+PC_2_S_speed_variance)/trace(LAMBDA_S_speed)

figure('Name','Scatterplots for male speed-data PC 1 and 2 scores')
plot(PC_1_S_speed_scores_X_variables,PC_2_S_speed_scores_X_variables,'+'),grid,axis equal  %X-scores same result
title('Scatterplots for male speed-data PC 1 and 2 scores','Fontsize',14)
xlabel('1st Principal Component scores')
ylabel('2nd Principal Component scores')

figure('Name','Bi-plot for male speed-data PC 1 and 2 scores and loadings') 
names = char('100','200','400','800','1500','5000','10000','Marathon');
biplot(E_S_speed(:,1:2),'scores',[PC_1_S_speed_scores_X_variables PC_2_S_speed_scores_X_variables],'varlabels',names)
title('Bi-plot for male speed-data PC 1 and 2 scores and loadings','Fontsize',14)

figure('Name','Bi-plot for male speed-data PC 1, 2 and 3 scores and loadings')
PC_3_S_speed_loadings_X_variables = (E_S_speed(:,3))';
PC_3_S_speed_scores_X_variables = (X_speed - ones(n,1)*mean(X_speed))*PC_3_S_speed_loadings_X_variables';
biplot(E_S_speed(:,1:3),'scores',[PC_1_S_speed_scores_X_variables PC_2_S_speed_scores_X_variables PC_3_S_speed_scores_X_variables],'varlabels',names)
title('Bi-plot for male speed-data PC 1, 2 and 3 scores and loadings','Fontsize',14)

disp('---------------------------------------------------------------------------------')
disp('Correlation between PC 1 and 2 and original variables')
disp('---------------------------------------------------------------------------------')
Correlations_PC_1_S_speed_versus_X_variables = sqrt(LAMBDA_S_speed(1,1))*PC_1_S_speed_loadings_X_variables./(sqrt(diag(S_speed)))'
Correlations_PC_2_S_speed_versus_X_variables = sqrt(LAMBDA_S_speed(2,2))*PC_2_S_speed_loadings_X_variables./(sqrt(diag(S_speed)))'
figure('Name','Circle diagram for correlation between speed-data PC 1 and 2 and original variables') 
t = 0:.01:pi;
x = cos(t);
y = sin(t);
plot(x,y,'k','LineWidth',3),axis square,title('Circle diagram for correlation between speed-data PC 1 and 2 and original variables','Fontsize',14),grid,xlim([-1.2 1.2]),ylim([-1.2 1.2])
xlabel('Correlation with PC 1','Fontsize',14),ylabel('Correlation with PC 2','Fontsize',14)
hold on
plot(x,-y,'k','LineWidth',3);
plot(Correlations_PC_1_S_speed_versus_X_variables,Correlations_PC_2_S_speed_versus_X_variables,'r+','Linewidth',10)
text(0.92,0.44,'100','Fontsize',12)
text(0.68,0.44,'200','Fontsize',12)
text(0.7,0.38,'400','Fontsize',12)
text(0.75,0.033,'800','Fontsize',12)
text(0.78,-0.05,'1500','Fontsize',12)
text(0.8,-0.18,'5000','Fontsize',12)
text(1,-0.22,'10000','Fontsize',12)
text(0.72,-0.27,'Marath','Fontsize',12)
text(-0.5,0.7,'Since        {\Sigma}_j_=_1_._._p r^2(X_i,PC_j) = 1','Fontsize',14)
text(-0.5,0.55,'it follows   r^2(X_i,PC_1) + r^2(X_i,PC_2) {\leq} 1','Fontsize',14)
hold off

disp('---------------------------------------------------------------------------------')
disp('Approx. simult. Bonferroni CI for variances of speed-data PC 1 and 2')
disp('---------------------------------------------------------------------------------')
alpha = 0.1
psim = 2;
alpha_bonf = alpha/psim
CI_approx_lambda_1 = LAMBDA_S_speed(1,1)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
CI_approx_lambda_2 = LAMBDA_S_speed(2,2)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
disp('---------------------------------------------------------------------------------')
%------------------------------------------------------------------------------------------------------------------------
