clc;clear;close all
format compact
%------------------------------------------------------------------------
% LDA: define three populations, different mu, same SIGMA
%------------------------------------------------------------------------
g = 3;
p = 2;
mu1 = [2 2];
SIGMA1 = [1 0.5; 0.5 1];
n1 = 50;
seed1 = 10;rng(seed1); 
X1 = mvnrnd(mu1,SIGMA1,n1);
mu2 = [5 3];
SIGMA2 = SIGMA1;
n2 = n1;
seed2 = 110;rng(seed2); 
X2 = mvnrnd(mu2,SIGMA2,n2);
mu3 = [3 5];
SIGMA3 = SIGMA1;
n3 = 50;
seed3 = 210;rng(seed3); 
X3 = mvnrnd(mu3,SIGMA3,n3);
%------------------------------------------------------------------------
% descriptive statistics
%------------------------------------------------------------------------
x1_avg = mean(X1);
S1 = cov(X1);
x2_avg = mean(X2);
S2 = cov(X2);
x3_avg = mean(X3);
S3 = cov(X3);
%------------------------------------------------------------------------
% plot samples and estimated prediction ellipses
%------------------------------------------------------------------------
figure('Name','LDA: define three populations, different mu, same SIGMA')
subplot(4,3,1:3)
axis off
text(0,1.1,'Discriminant Analysis (classification) for three (g=3) bivariate (p=2) populations, case I','Fontsize',24)
text(0,0.8,'Training data with known population membership:','Fontsize',18)
text(0,0.56,['Population i=1 sample: n_1 = ' num2str(n1) ' observations'],'Color','r','Fontsize',18)
text(0,0.36,['Population i=2 sample: n_2 = ' num2str(n2) ' observations'],'Color','b','Fontsize',18)
text(0,0.16,['Population i=3 sample: n_3 = ' num2str(n3) ' observations'],'Color','g','Fontsize',18)
alpha = [0.5 0.1];
text(0,-0.1,['Also is shown (1-{\alpha})*100% estimated prediction ellipses for  1-{\alpha}  = ' mat2str(1-alpha) '  (assuming bivariate normal distributions)'],'Fontsize',14)
subplot(4,3,4:12)
y1 = x1_avg(1)-3*sqrt(S1(1,1)):6*sqrt(S1(1,1))/50:x1_avg(1)+3*sqrt(S1(1,1)); 
y2 = x1_avg(2)-3*sqrt(S1(2,2)):6*sqrt(S1(2,2))/50:x1_avg(2)+3*sqrt(S1(2,2));
[Y1,Y2] = meshgrid(y1,y2);
F1 = mvnpdf([Y1(:) Y2(:)],x1_avg,S1);
F1 = reshape(F1,length(y2),length(y1));
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid,axis equal
xlabel('x_i_1_j,  (i=1: j=1,...,n_1)  (i=2: j=1,...,n_2)  (i=3: j=1,...,n_3)','Fontsize',16); ylabel('x_i_2_j','Fontsize',16)
hold on
c = chi2inv(1-alpha,p);
C1 = (1/(2*pi*sqrt(det(S1)))).*exp(-c/2);
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)%,axis square
%---------------------------------
z1 = x2_avg(1)-3*sqrt(S2(1,1)):6*sqrt(S2(1,1))/50:x2_avg(1)+3*sqrt(S2(1,1)); 
z2 = x2_avg(2)-3*sqrt(S2(2,2)):6*sqrt(S2(2,2))/50:x2_avg(2)+3*sqrt(S2(2,2));
[Z1,Z2] = meshgrid(z1,z2);
F2 = mvnpdf([Z1(:) Z2(:)],x2_avg,S2);
F2 = reshape(F2,length(z2),length(z1));
plot(X2(:,1),X2(:,2),'bo','LineWidth',3)
c = chi2inv(1-alpha,p);
C2 = (1/(2*pi*sqrt(det(S2)))).*exp(-c/2);
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)%,axis square
%---------------------------------
u1 = x3_avg(1)-3*sqrt(S3(1,1)):6*sqrt(S3(1,1))/50:x3_avg(1)+3*sqrt(S3(1,1)); 
u2 = x3_avg(2)-3*sqrt(S3(2,2)):6*sqrt(S3(2,2))/50:x3_avg(2)+3*sqrt(S3(2,2));
[U1,U2] = meshgrid(u1,u2);
F3 = mvnpdf([U1(:) U2(:)],x3_avg,S3);
F3 = reshape(F3,length(u2),length(u1));
plot(X3(:,1),X3(:,2),'go','LineWidth',3)
c = chi2inv(1-alpha,p);
C3 = (1/(2*pi*sqrt(det(S3)))).*exp(-c/2);
contour(u1,u2,F3,C3,'Color','g','LineWidth',1)%,axis square
hold off
%------------------------------------------------------------------------
% MVN model check for population 1 sample
%------------------------------------------------------------------------
figure('Name','MVN model check for population 1 sample')
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations from POPULATION 1 possibly being from bivariate normal distribution','Fontsize',20)
subplot(5,2,3)
axis off
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n_1 = ' num2str(n1)],'Fontsize',16)
text(0.5,0.75,['{\mu}_1_,_h_a_t = x_1 = ' mat2str(x1_avg,2)],'Fontsize',16)
text(0.635,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_1_,_h_a_t = S_1 = ' mat2str(S1,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X1,'ro')
xlabel('x_1_1_j                                          x_1_2_j','Fontsize',16)
ylabel('x_1_2_j                           x_1_1_j','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, j=1..n_1 :','Fontsize',18)
text(0,0.7,'d_1_j^2 = (x_1_j - x_1)^T S_1^-^1 (x_1_j - x_1) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.178,0.91,'_','Fontsize',16),warning('off')
text(0.408,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu1_hat = x1_avg;
d1j_sqr = zeros(1,n1);
for j = 1:n1
    d1j_sqr(j) = (X1(j,:)-mu1_hat)*inv(S1)*(X1(j,:)-mu1_hat)';
end
df = 2;
z_j = chi2rnd(df,1,n1);
% z2 = chi2rnd(20,1,n);
qqplot(d1j_sqr,z_j),grid,xlabel('quantiles for d_1_j^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_1_j^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% MVN model check for population 2 sample
%------------------------------------------------------------------------
figure('Name','MVN model check for population 2 sample')
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations from POPULATION 2 possibly being from bivariate normal distribution','Fontsize',20)
subplot(5,2,3)
axis off
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n_2 = ' num2str(n2)],'Fontsize',16)
text(0.5,0.75,['{\mu}_2_,_h_a_t = x_2 = ' mat2str(x2_avg,2)],'Fontsize',16)
text(0.635,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_2_,_h_a_t = S_2 = ' mat2str(S2,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X2,'bo')
xlabel('x_2_1_j                                          x_2_2_j','Fontsize',16)
ylabel('x_2_2_j                           x_2_1_j','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, j=1..n_2 :','Fontsize',18)
text(0,0.7,'d_2_j^2 = (x_2_j - x_2)^T S_2^-^1 (x_2_j - x_2) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.178,0.91,'_','Fontsize',16),warning('off')
text(0.408,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu2_hat = x2_avg;
d2j_sqr = zeros(1,n2);
for j = 1:n2
    d2j_sqr(j) = (X2(j,:)-mu2_hat)*inv(S2)*(X2(j,:)-mu2_hat)';
end
df = 2;
z_j = chi2rnd(df,1,n2);
% z2 = chi2rnd(20,1,n);
qqplot(d2j_sqr,z_j),grid,xlabel('quantiles for d_2_j^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_2_j^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% MVN model check for population 3 sample
%------------------------------------------------------------------------
figure('Name','MVN model check for population 3 sample')
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations from POPULATION 3 possibly being from bivariate normal distribution','Fontsize',20)
subplot(5,2,3)
axis off
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n_3 = ' num2str(n3)],'Fontsize',16)
text(0.5,0.75,['{\mu}_3_,_h_a_t = x_3 = ' mat2str(x3_avg,2)],'Fontsize',16)
text(0.635,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_3_,_h_a_t = S_3 = ' mat2str(S3,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X3,'go')
xlabel('x_3_1_j                                          x_3_2_j','Fontsize',16)
ylabel('x_3_2_j                           x_3_1_j','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, j=1..n_3 :','Fontsize',18)
text(0,0.7,'d_3_j^2 = (x_3_j - x_3)^T S_3^-^1 (x_3_j - x_3) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.178,0.91,'_','Fontsize',16),warning('off')
text(0.408,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu3_hat = x3_avg;
d3j_sqr = zeros(1,n3);
for j = 1:n3
    d3j_sqr(j) = (X3(j,:)-mu3_hat)*inv(S3)*(X3(j,:)-mu3_hat)';
end
df = 2;
z_j = chi2rnd(df,1,n3);
% z2 = chi2rnd(20,1,n);
qqplot(d2j_sqr,z_j),grid,xlabel('quantiles for d_3_j^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_3_j^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% Bartlett test for equal covariance matrices for both populations
%------------------------------------------------------------------------
figure('Name','Bartlett and MANOVA tests')
subplot(5,3,[4 7 10])
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid,axis equal
xlabel('x_i_1_j,  (i=1,2,3)  (j=1,...,n_i)','Fontsize',12); ylabel('x_i_2_j','Fontsize',12)
hold on
plot(X2(:,1),X2(:,2),'bo','LineWidth',3)
plot(X3(:,1),X3(:,2),'go','LineWidth',3)
title({'Thus having "verified" the assumptions';'  ';'X_1 ~ N_2({\mu}_1,{\Sigma}_1)    (RED)';'X_2 ~ N_2({\mu}_2,{\Sigma}_2)   (BLUE)';'X_3 ~ N_2({\mu}_3,{\Sigma}_3)   (GREEN)'},'Fontsize',14)
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)
contour(u1,u2,F3,C3,'Color','g','LineWidth',1)
hold off
subplot(5,3,8:9)
axis off
line([0 1],[0.5 0.5],'Color','k','LineWidth',3)
subplot(5,3,[2 3 5 6])
axis off
Sp = ((n1-1)*S1 + (n2-1)*S2 + (n3-1)*S3)/(n1+n2+n3-3);
T = (n1+n2+n3-3)*log(det(Sp)) - (n1-1)*log(det(S1)) - (n2-1)*log(det(S2)) - (n3-1)*log(det(S3));
correction_factor = 1 - ((2*p^2+3*p-1)/(6*(p+1)*(g-1)))*(1/(n1-1)+1/(n2-1)+1/(n3-1)-1/(n1+n2+n3-3));
test_statistic = correction_factor*T;
df = (g-1)*p*(p+1)/2;
p_value = 1-chi2cdf(test_statistic,df);
text(0,1,'Bartlett test of equal covariance matrices,  H_0: {\Sigma}_1 = {\Sigma}_2 = {\Sigma}_3','Fontsize',20)
text(0,0.85,['S_1 = ' mat2str(S1,2), ',    n_1 = ' num2str(n1)],'Color','r','Fontsize',14)
text(0,0.75,['S_2 = ' mat2str(S2,2), ',    n_2 = ' num2str(n2)],'Color','b','Fontsize',14)
text(0,0.65,['S_3 = ' mat2str(S3,2), ',    n_3 = ' num2str(n3)],'Color','g','Fontsize',14)
text(0,0.55,['S_p = [(n_1-1)S_1 + (n_2-1)S_2 + (n_3-1)S_3]/(n_1+n_2+n_3-3) = ' mat2str(Sp,2)],'Color','k','Fontsize',14)
text(0,0.4,'Test statistic, T = C*[(n_1+n_2+n_3-3)log|S_p| - (n_1-1)log|S_1| - (n_2-1)log|S_2| - (n_3-1)log|S_3|]  ~  {\chi}^2_p_(_p_+_1_)_/_2','Color','k','Fontsize',14)
text(0,0.28,'Finite population correction factor, C = 1 - ((2*p^2+3*p-1)/(6*(p+1)(g-1)))*(1/(n_1-1)+1/(n_2-1)+1/(n_3-1)-1/(n_1+n_+n_3-3))','Color','k','Fontsize',10)
text(0,0.1,['t = ' num2str(test_statistic,3)],'Color','k','Fontsize',14)
text(0,0,['p-value = P(T > t) = ' num2str(p_value,3) '   (where T ~ {\chi}^2_(_g_-_1_)_p_(_p_+_1_)_/_2)'],'Color','k','Fontsize',14)
text(0,-0.2,'CONCLUSION:  H_0 is ACCEPTED  =>  LINEAR Discriminant Analysis case','Fontsize',14)
%------------------------------------------------------------------------
% MANOVA test for significant separation of populations
%------------------------------------------------------------------------
subplot(5,3,[11 12 14 15])
axis off
mu_hat_1 = x1_avg;
mu_hat_2 = x2_avg;
mu_hat_3 = x3_avg;
mu_hat = mean([x1_avg; x2_avg; x3_avg]);
SSB_1 = n1*(mu_hat_1 - mu_hat)'*(mu_hat_1 - mu_hat);
SSB_2 = n2*(mu_hat_2 - mu_hat)'*(mu_hat_2 - mu_hat);
SSB_3 = n3*(mu_hat_3 - mu_hat)'*(mu_hat_3 - mu_hat);
SSB = SSB_1 + SSB_2 + SSB_3;
df_B = g-1;
SB = SSB/df_B;
SSW_1 = (X1 - ones(n1,1)*mu_hat_1)'*(X1 - ones(n1,1)*mu_hat_1);
SSW_2 = (X2 - ones(n2,1)*mu_hat_2)'*(X2 - ones(n2,1)*mu_hat_2);
SSW_3 = (X3 - ones(n3,1)*mu_hat_3)'*(X3 - ones(n3,1)*mu_hat_3);
SSW = SSW_1 + SSW_2 + SSW_3;
n = n1 + n2 + n3;
df_W = n-g;
SW = SSW/df_W;
LAMBDA = det(SSW)/det(SSB + SSW);
test_statistic = ((n-p-2)/p)*(1-sqrt(LAMBDA))/sqrt(LAMBDA);
p_value = 1-fcdf(test_statistic,2*p,2*(n-p-2));
text(0,1.2,'MANOVA test of significant separation of populations,  H_0: {\mu}_1 = {\mu}_2 = {\mu}_3','Fontsize',20)
text(0,1,['x_1_,_a_v_g = ' mat2str(x1_avg,2), ',    n_1 = ' num2str(n1)],'Color','r','Fontsize',14)
text(0,0.9,['x_2_,_a_v_g = ' mat2str(x2_avg,2), ',    n_2 = ' num2str(n2)],'Color','b','Fontsize',14)
text(0,0.8,['x_3_,_a_v_g = ' mat2str(x3_avg,2), ',    n_3 = ' num2str(n3)],'Color','g','Fontsize',14)
text(0,0.65,['{\Lambda*} = |SS_W| / |SS_B + SS_W| = ' num2str(LAMBDA,2)],'Color','k','Fontsize',14)
text(0,0.5,'Test statistic, T = [(n_1 + n_2 + n_3 - p - 2) / p] * [(1-sqrt({\Lambda*}) / sqrt({\Lambda*})]  ~  F_2_p_,_2_(_n_1_+_n_2_+_n_3_-_p_-_2_)','Color','k','Fontsize',14)
text(0,0.33,['t = ' num2str(test_statistic,3)],'Color','k','Fontsize',14)
text(0,0.2,['p-value = P(T > t) = ' num2str(p_value,3) '   (where T ~ F_2_p_,_2_(_n_1_+_n_2_+_n_3_-_p_-_2_))'],'Color','k','Fontsize',14)
text(0,0,'CONCLUSION:  H_0 is REJECTED  =>  Discriminant Analysis is relevant (significant separation of populations)','Fontsize',14)
%------------------------------------------------------------------------
% Linear discriminant function for equal priors and costs
%------------------------------------------------------------------------
figure('Name','Linear discriminant function for equal priors and costs')
subplot(4,7,1:3)
axis off
text(0,1.2,'Classification for g=3 MVN populations with equal misclassification costs and equal priors','Fontsize',24)
text(0,0.8,'Pairwise linear discriminant functions for new observation x_0:','Fontsize',16)
text(0,0.6,'d_i_k(x_0) = (x_i_,_a_v_g- x_k_,_a_v_g)^T S_p^-^1 [x_0 - 0.5(x_i_,_a_v_g+ x_k_,_a_v_g)],  i,k=1,2,3','Fontsize',16)
text(0,0.27,'Decision rule:','Fontsize',16)
text(0.3,0.25,'d_1_2(x_0) > 0 {\wedge} d_1_3(x_0) > 0  {\rightarrow}  d_1','Fontsize',16)
text(0.3,0.03,'d_2_1(x_0) > 0 {\wedge} d_2_3(x_0) > 0  {\rightarrow}  d_2','Fontsize',16)
text(0.3,-0.19,'d_3_1(x_0) > 0 {\wedge} d_3_2(x_0) > 0  {\rightarrow}  d_3','Fontsize',16)
subplot(4,7,4:7:25)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
annotation('arrow',[0.483125 0.559375],...
    [0.572979591836735 0.573979591836735],'HeadLength',20,'HeadWidth',20,...
    'LineWidth',6);
subplot(4,7,[8:10 15:17 22:24])
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid,axis equal,xlim([-1 8]),ylim([-1 8])
xlabel('x_i_1_j,  (i=1,2,3)  (j=1,...,n_i)','Fontsize',12); ylabel('x_i_2_j','Fontsize',12)
hold on
plot(X2(:,1),X2(:,2),'bo','LineWidth',3)
plot(X3(:,1),X3(:,2),'go','LineWidth',3)
text(-4,7.5,'d_1_2(x) = 0','Color','k','Fontsize',20)
text(-4,6.5,'d_1_3(x) = 0','Color','m','Fontsize',20)
text(-4,5.5,'d_2_3(x) = 0','Color','c','Fontsize',20)
annotation('arrow',[0.13125 0.163125],...
    [0.677571428571428 0.677295918367347],'HeadLength',15,'HeadWidth',15,...
    'LineWidth',6);
annotation('arrow',[0.1325 0.164375],...
    [0.612928571428571 0.612653061224489],'HeadLength',15,'HeadWidth',15,...
    'LineWidth',6,...
    'Color',[1 0 1]);
annotation('arrow',[0.13375 0.165625],...
    [0.545326530612243 0.545051020408162],'HeadLength',15,'HeadWidth',15,...
    'LineWidth',6,...
    'Color',[0 1 1]);
text(-4,1.5,'X_1~N_2({\mu}_1,{\Sigma})','Color','r','Fontsize',16)
text(-4,0.8,'X_2~N_2({\mu}_2,{\Sigma})','Color','b','Fontsize',16)
text(-4,0.1,'X_3~N_2({\mu}_3,{\Sigma})','Color','g','Fontsize',16)
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)
contour(u1,u2,F3,C3,'Color','g','LineWidth',1)
v1_min = min([min(y1),min(z1),min(u1)]);
v1_max = max([max(y1),max(z1),max(u1)]);
v1 = v1_min:(v1_max-v1_min)/100:v1_max;
v2_min = min([min(y2),min(z2),min(u2)]);
v2_max = max([max(y2),max(z2),max(u2)]);
v2 = v2_min:(v2_max-v2_min)/100:v2_max;
[V1,V2] = meshgrid(v1,v2);
F12 = (x1_avg-x2_avg)*inv(Sp)*[V1(:)' ; V2(:)'] - 0.5*(x1_avg-x2_avg)*inv(Sp)*(x1_avg+x2_avg)';
F12 = reshape(F12,length(v2),length(v1));
contour(v1,v2,F12,[eps eps],'Color','k','LineWidth',3)
F13 = (x1_avg-x3_avg)*inv(Sp)*[V1(:)' ; V2(:)'] - 0.5*(x1_avg-x3_avg)*inv(Sp)*(x1_avg+x3_avg)';
F13 = reshape(F13,length(v2),length(v1));
contour(v1,v2,F13,[eps eps],'Color','m','LineWidth',3)
F23 = (x2_avg-x3_avg)*inv(Sp)*[V1(:)' ; V2(:)'] - 0.5*(x2_avg-x3_avg)*inv(Sp)*(x2_avg+x3_avg)';
F23 = reshape(F23,length(v2),length(v1));
contour(v1,v2,F23,[eps eps],'Color','c','LineWidth',3)
hold off
%------------------------------------------------------------------------
% Calculate confusion matrix and APER for training data
%------------------------------------------------------------------------
LDA_12_scores_pop1 = (x1_avg-x2_avg)*inv(Sp)*X1' - 0.5*(x1_avg-x2_avg)*inv(Sp)*(x1_avg+x2_avg)';
LDA_13_scores_pop1 = (x1_avg-x3_avg)*inv(Sp)*X1' - 0.5*(x1_avg-x3_avg)*inv(Sp)*(x1_avg+x3_avg)';
OK_classify_pop1 = find(LDA_12_scores_pop1 > 0 & LDA_13_scores_pop1 > 0);
n1C_APER = length(OK_classify_pop1);
not_OK_classify_pop1 = find(LDA_12_scores_pop1 < 0 | LDA_13_scores_pop1 < 0);
LDA_23_scores_pop1 = (x2_avg-x3_avg)*inv(Sp)*X1(not_OK_classify_pop1,:)' - 0.5*(x2_avg-x3_avg)*inv(Sp)*(x2_avg+x3_avg)';
classify_as_2_pop1 = find(LDA_23_scores_pop1 > 0);
n1M2_APER = length(classify_as_2_pop1);
classify_as_3_pop1 = find(LDA_23_scores_pop1 < 0);
n1M3_APER = length(classify_as_3_pop1);
%----------------------------------
LDA_21_scores_pop2 = (x2_avg-x1_avg)*inv(Sp)*X2' - 0.5*(x2_avg-x1_avg)*inv(Sp)*(x2_avg+x1_avg)';
LDA_23_scores_pop2 = (x2_avg-x3_avg)*inv(Sp)*X2' - 0.5*(x2_avg-x3_avg)*inv(Sp)*(x2_avg+x3_avg)';
OK_classify_pop2 = find(LDA_21_scores_pop2 > 0 & LDA_23_scores_pop2 > 0);
n2C_APER = length(OK_classify_pop2);
not_OK_classify_pop2 = find(LDA_21_scores_pop2 < 0 | LDA_23_scores_pop2 < 0);
LDA_13_scores_pop2 = (x1_avg-x3_avg)*inv(Sp)*X2(not_OK_classify_pop2,:)' - 0.5*(x1_avg-x3_avg)*inv(Sp)*(x1_avg+x3_avg)';
classify_as_1_pop2 = find(LDA_13_scores_pop2 > 0);
n2M1_APER = length(classify_as_1_pop2);
classify_as_3_pop2 = find(LDA_13_scores_pop2 < 0);
n2M3_APER = length(classify_as_3_pop2);
%----------------------------------
LDA_31_scores_pop3 = (x3_avg-x1_avg)*inv(Sp)*X3' - 0.5*(x3_avg-x1_avg)*inv(Sp)*(x3_avg+x1_avg)';
LDA_32_scores_pop3 = (x3_avg-x2_avg)*inv(Sp)*X3' - 0.5*(x3_avg-x2_avg)*inv(Sp)*(x3_avg+x2_avg)';
OK_classify_pop3 = find(LDA_31_scores_pop3 > 0 & LDA_32_scores_pop3 > 0);
n3C_APER = length(OK_classify_pop3);
not_OK_classify_pop3 = find(LDA_31_scores_pop3 < 0 | LDA_32_scores_pop3 < 0);
LDA_12_scores_pop3 = (x1_avg-x2_avg)*inv(Sp)*X3(not_OK_classify_pop3,:)' - 0.5*(x1_avg-x2_avg)*inv(Sp)*(x1_avg+x2_avg)';
classify_as_1_pop3 = find(LDA_12_scores_pop3 > 0);
n3M1_APER = length(classify_as_1_pop3);
classify_as_2_pop3 = find(LDA_12_scores_pop3 < 0);
n3M2_APER = length(classify_as_2_pop3);
%----------------------------------
APER = (n1M2_APER+n1M3_APER+n2M1_APER+n2M3_APER+n3M1_APER+n3M2_APER)/(n1+n2+n3);
%------------------------------------------------------------------------
% Calculate confusion matrix and AER_hat for training data using
% cross-validation (hold-out) procedure
%------------------------------------------------------------------------
n1C_AER = 0;
n1M2_AER = 0;
n1M3_AER = 0;
for j1 = 1:n1
    X1_red = X1([1:(j1-1) j1+1:n1],:);
    x1_avg_red = mean(X1_red);
    S1_red = cov(X1_red);
    Sp_red = ((n1-2)*S1_red + (n2-1)*S2 + (n3-1)*S3)/(n1+n2+n3-4);
    LDA_12_score_heldout = (x1_avg_red-x2_avg)*inv(Sp_red)*X1(j1,:)' - 0.5*(x1_avg_red-x2_avg)*inv(Sp_red)*(x1_avg_red+x2_avg)';
    LDA_13_score_heldout = (x1_avg_red-x3_avg)*inv(Sp_red)*X1(j1,:)' - 0.5*(x1_avg_red-x3_avg)*inv(Sp_red)*(x1_avg_red+x3_avg)';
    LDA_23_score_heldout = (x2_avg-x3_avg)*inv(Sp_red)*X1(j1,:)' - 0.5*(x2_avg-x3_avg)*inv(Sp_red)*(x2_avg+x3_avg)';
    if (LDA_12_score_heldout > 0 && LDA_13_score_heldout > 0)
        n1C_AER = n1C_AER + 1;
    else if LDA_23_score_heldout > 0    
        n1M2_AER = n1M2_AER + 1;
        else
            n1M3_AER = n1M3_AER + 1;
        end
    end
end
%----------------------------------
n2C_AER = 0;
n2M1_AER = 0;
n2M3_AER = 0;
for j2 = 1:n2
    X2_red = X2([1:(j2-1) j2+1:n2],:);
    x2_avg_red = mean(X2_red);
    S2_red = cov(X2_red);
    Sp_red = ((n1-1)*S1 + (n2-2)*S2_red+ (n3-1)*S3)/(n1+n2+n3-4);
    LDA_21_score_heldout = (x2_avg_red-x1_avg)*inv(Sp_red)*X2(j2,:)' - 0.5*(x2_avg_red-x1_avg)*inv(Sp_red)*(x2_avg_red+x1_avg)';
    LDA_23_score_heldout = (x2_avg_red-x3_avg)*inv(Sp_red)*X2(j2,:)' - 0.5*(x2_avg_red-x3_avg)*inv(Sp_red)*(x2_avg_red+x3_avg)';
    LDA_13_score_heldout = (x1_avg-x3_avg)*inv(Sp_red)*X2(j2,:)' - 0.5*(x1_avg-x3_avg)*inv(Sp_red)*(x1_avg+x3_avg)';
    if (LDA_21_score_heldout > 0 && LDA_23_score_heldout > 0)
        n2C_AER = n2C_AER + 1;
    else if LDA_13_score_heldout > 0    
        n2M1_AER = n2M1_AER + 1;
        else
            n2M3_AER = n2M3_AER + 1;
        end
    end
end
%----------------------------------
n3C_AER = 0;
n3M1_AER = 0;
n3M2_AER = 0;
for j3 = 1:n3
    X3_red = X3([1:(j3-1) j3+1:n3],:);
    x3_avg_red = mean(X3_red);
    S3_red = cov(X3_red);
    Sp_red = ((n1-1)*S1 + (n2-1)*S2+ (n3-2)*S3_red)/(n1+n2+n3-4);
    LDA_31_score_heldout = (x3_avg_red-x1_avg)*inv(Sp_red)*X3(j3,:)' - 0.5*(x3_avg_red-x1_avg)*inv(Sp_red)*(x3_avg_red+x1_avg)';
    LDA_32_score_heldout = (x3_avg_red-x2_avg)*inv(Sp_red)*X3(j3,:)' - 0.5*(x3_avg_red-x2_avg)*inv(Sp_red)*(x3_avg_red+x2_avg)';
    LDA_12_score_heldout = (x1_avg-x2_avg)*inv(Sp_red)*X3(j3,:)' - 0.5*(x1_avg-x2_avg)*inv(Sp_red)*(x1_avg+x2_avg)';
    if (LDA_31_score_heldout > 0 && LDA_32_score_heldout > 0)
        n3C_AER = n3C_AER + 1;
    else if LDA_12_score_heldout > 0    
        n3M1_AER = n3M1_AER + 1;
        else
            n3M2_AER = n3M2_AER + 1;
        end
    end
end
%----------------------------------
AER_hat = (n1M2_AER+n1M3_AER+n2M1_AER+n2M3_AER+n3M1_AER+n3M2_AER)/(n1+n2+n3);
%------------------------------------------------------------------------
% Printout confusion tables and error rates
%------------------------------------------------------------------------
subplot(4,7,[5:7 12:14 19:21 26:28])
axis off
text(0,0.95,'Performance of classification','Fontsize',22)
text(0,0.91,'(empirical error rates)','Fontsize',22)
text(0,0.85,'------------------------------------------','Fontsize',20)
text(0,0.79,'Confusion matrix for training data:','Fontsize',18)
text(0.37,0.72,'classified as','Fontsize',16)
text(0.4,0.67,'d_1    d_2    d_3','Fontsize',16)
text(0.1,0.59,'true','Fontsize',16)
text(0.05,0.56,'population','Fontsize',16)
text(0.32,0.62,'{\pi}_1','Fontsize',16)
text(0.4,0.62,num2str(n1C_APER),'Fontsize',16)
text(0.5,0.62,num2str(n1M2_APER),'Fontsize',16)
text(0.59,0.62,num2str(n1M3_APER),'Fontsize',16)
text(0.68,0.62,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',14)
text(0.32,0.58,'{\pi}_2','Fontsize',16)
text(0.41,0.58,num2str(n2M1_APER),'Fontsize',16)
text(0.49,0.58,num2str(n2C_APER),'Fontsize',16)
text(0.59,0.58,num2str(n2M3_APER),'Fontsize',16)
text(0.68,0.58,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',14)
text(0.32,0.54,'{\pi}_3','Fontsize',16)
text(0.41,0.54,num2str(n3M1_APER),'Fontsize',16)
text(0.49,0.54,num2str(n3M2_APER),'Fontsize',16)
text(0.58,0.54,num2str(n3C_APER),'Fontsize',16)
text(0.68,0.54,['(sum = n_3 = ' num2str(n3) ')'],'Fontsize',14)
text(0,0.47,['Apparent Error Rate, APER = ' num2str(n1M2_APER+n1M3_APER+n2M1_APER+n2M3_APER+n3M1_APER+n3M2_APER) '/' num2str(n1+n1+n3) ' = ' num2str(APER,3)],'Fontsize',16)
text(0,0.42,'------------------------------------------','Fontsize',20)
text(0,0.36,'Confusion matrix for training data:','Fontsize',18)
text(0,0.32,'using cross-validation (hold-out):','Fontsize',18)
text(0.37,0.25,'classified as','Fontsize',16)
text(0.4,0.2,'d_1    d_2    d_3','Fontsize',16)
text(0.1,0.12,'true','Fontsize',16)
text(0.05,0.09,'population','Fontsize',16)
text(0.32,0.15,'{\pi}_1','Fontsize',16)
text(0.4,0.15,num2str(n1C_AER),'Fontsize',16)
text(0.5,0.15,num2str(n1M2_AER),'Fontsize',16)
text(0.59,0.15,num2str(n1M3_AER),'Fontsize',16)
text(0.68,0.15,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',14)
text(0.32,0.11,'{\pi}_2','Fontsize',16)
text(0.41,0.11,num2str(n2M1_AER),'Fontsize',16)
text(0.49,0.11,num2str(n2C_AER),'Fontsize',16)
text(0.59,0.11,num2str(n2M3_AER),'Fontsize',16)
text(0.68,0.11,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',14)
text(0.32,0.07,'{\pi}_3','Fontsize',16)
text(0.41,0.07,num2str(n3M1_AER),'Fontsize',16)
text(0.49,0.07,num2str(n3M2_AER),'Fontsize',16)
text(0.58,0.07,num2str(n3C_AER),'Fontsize',16)
text(0.68,0.07,['(sum = n_3 = ' num2str(n3) ')'],'Fontsize',14)
text(0,0,['Estimated Actual Error Rate, AER_h_a_t = ' num2str(n1M2_AER+n1M3_AER+n2M1_AER+n2M3_AER+n3M1_AER+n3M2_AER) '/' num2str(n1+n1+n3) ' = ' num2str(AER_hat,3)],'Fontsize',16)
text(0,-0.05,'------------------------------------------','Fontsize',20)
%-------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------
% QDA: define three populations, different mu, different SIGMA
%------------------------------------------------------------------------
g = 3;
p = 2;
mu1 = [2 2];
SIGMA1 = [1 0.5; 0.5 1];
n1 = 50;
seed1 = 310;rng(seed1); 
X1 = mvnrnd(mu1,SIGMA1,n1);
mu2 = [5 3];
SIGMA2 = [1 0; 0 1];
n2 = n1;
seed2 = 410;rng(seed2); 
X2 = mvnrnd(mu2,SIGMA2,n2);
mu3 = [3 5];
SIGMA3 = [1 -0.5; -0.5 1];
n3 = 50;
seed3 = 510;rng(seed3); 
X3 = mvnrnd(mu3,SIGMA3,n3);
%------------------------------------------------------------------------
% descriptive statistics
%------------------------------------------------------------------------
x1_avg = mean(X1);
S1 = cov(X1);
x2_avg = mean(X2);
S2 = cov(X2);
x3_avg = mean(X3);
S3 = cov(X3);
%------------------------------------------------------------------------
% plot samples and estimated prediction ellipses
%------------------------------------------------------------------------
figure('Name','QDA: define three populations, different mu, different SIGMA')
subplot(4,3,1:3)
axis off
text(0,1.1,'Discriminant Analysis (classification) for three (g=3) bivariate (p=2) populations, case II','Fontsize',24)
text(0,0.8,'Training data with known population membership:','Fontsize',18)
text(0,0.56,['Population i=1 sample: n_1 = ' num2str(n1) ' observations'],'Color','r','Fontsize',18)
text(0,0.36,['Population i=2 sample: n_2 = ' num2str(n2) ' observations'],'Color','b','Fontsize',18)
text(0,0.16,['Population i=3 sample: n_3 = ' num2str(n3) ' observations'],'Color','g','Fontsize',18)
alpha = [0.5 0.1];
text(0,-0.1,['Also is shown (1-{\alpha})*100% estimated prediction ellipses for  1-{\alpha}  = ' mat2str(1-alpha) '  (assuming bivariate normal distributions)'],'Fontsize',14)
subplot(4,3,4:12)
y1 = x1_avg(1)-3*sqrt(S1(1,1)):6*sqrt(S1(1,1))/50:x1_avg(1)+3*sqrt(S1(1,1)); 
y2 = x1_avg(2)-3*sqrt(S1(2,2)):6*sqrt(S1(2,2))/50:x1_avg(2)+3*sqrt(S1(2,2));
[Y1,Y2] = meshgrid(y1,y2);
F1 = mvnpdf([Y1(:) Y2(:)],x1_avg,S1);
F1 = reshape(F1,length(y2),length(y1));
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid,axis equal
xlabel('x_i_1_j,  (i=1: j=1,...,n_1)  (i=2: j=1,...,n_2)  (i=3: j=1,...,n_3)','Fontsize',16); ylabel('x_i_2_j','Fontsize',16)
hold on
c = chi2inv(1-alpha,p);
C1 = (1/(2*pi*sqrt(det(S1)))).*exp(-c/2);
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)%,axis square
%---------------------------------
z1 = x2_avg(1)-3*sqrt(S2(1,1)):6*sqrt(S2(1,1))/50:x2_avg(1)+3*sqrt(S2(1,1)); 
z2 = x2_avg(2)-3*sqrt(S2(2,2)):6*sqrt(S2(2,2))/50:x2_avg(2)+3*sqrt(S2(2,2));
[Z1,Z2] = meshgrid(z1,z2);
F2 = mvnpdf([Z1(:) Z2(:)],x2_avg,S2);
F2 = reshape(F2,length(z2),length(z1));
plot(X2(:,1),X2(:,2),'bo','LineWidth',3)
c = chi2inv(1-alpha,p);
C2 = (1/(2*pi*sqrt(det(S2)))).*exp(-c/2);
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)%,axis square
%---------------------------------
u1 = x3_avg(1)-3*sqrt(S3(1,1)):6*sqrt(S3(1,1))/50:x3_avg(1)+3*sqrt(S3(1,1)); 
u2 = x3_avg(2)-3*sqrt(S3(2,2)):6*sqrt(S3(2,2))/50:x3_avg(2)+3*sqrt(S3(2,2));
[U1,U2] = meshgrid(u1,u2);
F3 = mvnpdf([U1(:) U2(:)],x3_avg,S3);
F3 = reshape(F3,length(u2),length(u1));
plot(X3(:,1),X3(:,2),'go','LineWidth',3)
c = chi2inv(1-alpha,p);
C3 = (1/(2*pi*sqrt(det(S3)))).*exp(-c/2);
contour(u1,u2,F3,C3,'Color','g','LineWidth',1)%,axis square
hold off
%------------------------------------------------------------------------
% MVN model check for population 1 sample
%------------------------------------------------------------------------
figure('Name','MVN model check for population 1 sample')
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations from POPULATION 1 possibly being from bivariate normal distribution','Fontsize',20)
subplot(5,2,3)
axis off
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n_1 = ' num2str(n1)],'Fontsize',16)
text(0.5,0.75,['{\mu}_1_,_h_a_t = x_1 = ' mat2str(x1_avg,2)],'Fontsize',16)
text(0.635,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_1_,_h_a_t = S_1 = ' mat2str(S1,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X1,'ro')
xlabel('x_1_1_j                                          x_1_2_j','Fontsize',16)
ylabel('x_1_2_j                           x_1_1_j','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, j=1..n_1 :','Fontsize',18)
text(0,0.7,'d_1_j^2 = (x_1_j - x_1)^T S_1^-^1 (x_1_j - x_1) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.178,0.91,'_','Fontsize',16),warning('off')
text(0.408,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu1_hat = x1_avg;
d1j_sqr = zeros(1,n1);
for j = 1:n1
    d1j_sqr(j) = (X1(j,:)-mu1_hat)*inv(S1)*(X1(j,:)-mu1_hat)';
end
df = 2;
z_j = chi2rnd(df,1,n1);
% z2 = chi2rnd(20,1,n);
qqplot(d1j_sqr,z_j),grid,xlabel('quantiles for d_1_j^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_1_j^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% MVN model check for population 2 sample
%------------------------------------------------------------------------
figure('Name','MVN model check for population 2 sample')
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations from POPULATION 2 possibly being from bivariate normal distribution','Fontsize',20)
subplot(5,2,3)
axis off
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n_2 = ' num2str(n2)],'Fontsize',16)
text(0.5,0.75,['{\mu}_2_,_h_a_t = x_2 = ' mat2str(x2_avg,2)],'Fontsize',16)
text(0.635,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_2_,_h_a_t = S_2 = ' mat2str(S2,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X2,'bo')
xlabel('x_2_1_j                                          x_2_2_j','Fontsize',16)
ylabel('x_2_2_j                           x_2_1_j','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, j=1..n_2 :','Fontsize',18)
text(0,0.7,'d_2_j^2 = (x_2_j - x_2)^T S_2^-^1 (x_2_j - x_2) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.178,0.91,'_','Fontsize',16),warning('off')
text(0.408,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu2_hat = x2_avg;
d2j_sqr = zeros(1,n2);
for j = 1:n2
    d2j_sqr(j) = (X2(j,:)-mu2_hat)*inv(S2)*(X2(j,:)-mu2_hat)';
end
df = 2;
z_j = chi2rnd(df,1,n2);
% z2 = chi2rnd(20,1,n);
qqplot(d2j_sqr,z_j),grid,xlabel('quantiles for d_2_j^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_2_j^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% MVN model check for population 3 sample
%------------------------------------------------------------------------
figure('Name','MVN model check for population 3 sample')
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations from POPULATION 3 possibly being from bivariate normal distribution','Fontsize',20)
subplot(5,2,3)
axis off
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n_3 = ' num2str(n3)],'Fontsize',16)
text(0.5,0.75,['{\mu}_3_,_h_a_t = x_3 = ' mat2str(x3_avg,2)],'Fontsize',16)
text(0.635,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_3_,_h_a_t = S_3 = ' mat2str(S3,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X3,'go')
xlabel('x_3_1_j                                          x_3_2_j','Fontsize',16)
ylabel('x_3_2_j                           x_3_1_j','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, j=1..n_3 :','Fontsize',18)
text(0,0.7,'d_3_j^2 = (x_3_j - x_3)^T S_3^-^1 (x_3_j - x_3) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.178,0.91,'_','Fontsize',16),warning('off')
text(0.408,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu3_hat = x3_avg;
d3j_sqr = zeros(1,n3);
for j = 1:n3
    d3j_sqr(j) = (X3(j,:)-mu3_hat)*inv(S3)*(X3(j,:)-mu3_hat)';
end
df = 2;
z_j = chi2rnd(df,1,n3);
% z2 = chi2rnd(20,1,n);
qqplot(d2j_sqr,z_j),grid,xlabel('quantiles for d_3_j^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_3_j^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% Bartlett test for equal covariance matrices for both populations
%------------------------------------------------------------------------
figure('Name','Bartlett and MANOVA tests')
subplot(5,3,[4 7 10])
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid,axis equal
xlabel('x_i_1_j,  (i=1,2,3)  (j=1,...,n_i)','Fontsize',12); ylabel('x_i_2_j','Fontsize',12)
hold on
plot(X2(:,1),X2(:,2),'bo','LineWidth',3)
plot(X3(:,1),X3(:,2),'go','LineWidth',3)
title({'Thus having "verified" the assumptions';'  ';'X_1 ~ N_2({\mu}_1,{\Sigma}_1)    (RED)';'X_2 ~ N_2({\mu}_2,{\Sigma}_2)   (BLUE)';'X_3 ~ N_2({\mu}_3,{\Sigma}_3)   (GREEN)'},'Fontsize',14)
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)
contour(u1,u2,F3,C3,'Color','g','LineWidth',1)
hold off
subplot(5,3,8:9)
axis off
line([0 1],[0.5 0.5],'Color','k','LineWidth',3)
subplot(5,3,[2 3 5 6])
axis off
Sp = ((n1-1)*S1 + (n2-1)*S2 + (n3-1)*S3)/(n1+n2+n3-3);
T = (n1+n2+n3-3)*log(det(Sp)) - (n1-1)*log(det(S1)) - (n2-1)*log(det(S2)) - (n3-1)*log(det(S3));
correction_factor = 1 - ((2*p^2+3*p-1)/(6*(p+1)*(g-1)))*(1/(n1-1)+1/(n2-1)+1/(n3-1)-1/(n1+n2+n3-3));
test_statistic = correction_factor*T;
df = (g-1)*p*(p+1)/2;
p_value = 1-chi2cdf(test_statistic,df);
text(0,1,'Bartlett test of equal covariance matrices,  H_0: {\Sigma}_1 = {\Sigma}_2 = {\Sigma}_3','Fontsize',20)
text(0,0.85,['S_1 = ' mat2str(S1,2), ',    n_1 = ' num2str(n1)],'Color','r','Fontsize',14)
text(0,0.75,['S_2 = ' mat2str(S2,2), ',    n_2 = ' num2str(n2)],'Color','b','Fontsize',14)
text(0,0.65,['S_3 = ' mat2str(S3,2), ',    n_3 = ' num2str(n3)],'Color','g','Fontsize',14)
text(0,0.55,['S_p = [(n_1-1)S_1 + (n_2-1)S_2 + (n_3-1)S_3]/(n_1+n_2+n_3-3) = ' mat2str(Sp,2)],'Color','k','Fontsize',14)
text(0,0.4,'Test statistic, T = C*[(n_1+n_2+n_3-3)log|S_p| - (n_1-1)log|S_1| - (n_2-1)log|S_2| - (n_3-1)log|S_3|]  ~  {\chi}^2_p_(_p_+_1_)_/_2','Color','k','Fontsize',14)
text(0,0.28,'Finite population correction factor, C = 1 - ((2*p^2+3*p-1)/(6*(p+1)(g-1)))*(1/(n_1-1)+1/(n_2-1)+1/(n_3-1)-1/(n_1+n_+n_3-3))','Color','k','Fontsize',10)
text(0,0.1,['t = ' num2str(test_statistic,3)],'Color','k','Fontsize',14)
text(0,0,['p-value = P(T > t) = ' num2str(p_value,3) '   (where T ~ {\chi}^2_(_g_-_1_)_p_(_p_+_1_)_/_2)'],'Color','k','Fontsize',14)
text(0,-0.2,'CONCLUSION:  H_0 is REJECTED  =>  QUADRATIC Discriminant Analysis case','Fontsize',14)
%------------------------------------------------------------------------
% MANOVA test for significant separation of populations
%------------------------------------------------------------------------
subplot(5,3,[11 12 14 15])
axis off
mu_hat_1 = x1_avg;
mu_hat_2 = x2_avg;
mu_hat_3 = x3_avg;
mu_hat = mean([x1_avg; x2_avg; x3_avg]);
SSB_1 = n1*(mu_hat_1 - mu_hat)'*(mu_hat_1 - mu_hat);
SSB_2 = n2*(mu_hat_2 - mu_hat)'*(mu_hat_2 - mu_hat);
SSB_3 = n3*(mu_hat_3 - mu_hat)'*(mu_hat_3 - mu_hat);
SSB = SSB_1 + SSB_2 + SSB_3;
df_B = g-1;
SB = SSB/df_B;
SSW_1 = (X1 - ones(n1,1)*mu_hat_1)'*(X1 - ones(n1,1)*mu_hat_1);
SSW_2 = (X2 - ones(n2,1)*mu_hat_2)'*(X2 - ones(n2,1)*mu_hat_2);
SSW_3 = (X3 - ones(n3,1)*mu_hat_3)'*(X3 - ones(n3,1)*mu_hat_3);
SSW = SSW_1 + SSW_2 + SSW_3;
n = n1 + n2 + n3;
df_W = n-g;
SW = SSW/df_W;
LAMBDA = det(SSW)/det(SSB + SSW);
test_statistic = ((n-p-2)/p)*(1-sqrt(LAMBDA))/sqrt(LAMBDA);
p_value = 1-fcdf(test_statistic,2*p,2*(n-p-2));
text(0,1.2,'MANOVA test of significant separation of populations,  H_0: {\mu}_1 = {\mu}_2 = {\mu}_3','Fontsize',20)
text(0,1,['x_1_,_a_v_g = ' mat2str(x1_avg,2), ',    n_1 = ' num2str(n1)],'Color','r','Fontsize',14)
text(0,0.9,['x_2_,_a_v_g = ' mat2str(x2_avg,2), ',    n_2 = ' num2str(n2)],'Color','b','Fontsize',14)
text(0,0.8,['x_3_,_a_v_g = ' mat2str(x3_avg,2), ',    n_3 = ' num2str(n3)],'Color','g','Fontsize',14)
text(0,0.65,['{\Lambda*} = |SS_W| / |SS_B + SS_W| = ' num2str(LAMBDA,2)],'Color','k','Fontsize',14)
text(0,0.5,'Test statistic, T = [(n_1 + n_2 + n_3 - p - 2) / p] * [(1-sqrt({\Lambda*}) / sqrt({\Lambda*})]  ~  F_2_p_,_2_(_n_1_+_n_2_+_n_3_-_p_-_2_)','Color','k','Fontsize',14)
text(0,0.33,['t = ' num2str(test_statistic,3)],'Color','k','Fontsize',14)
text(0,0.2,['p-value = P(T > t) = ' num2str(p_value,3) '   (where T ~ F_2_p_,_2_(_n_1_+_n_2_+_n_3_-_p_-_2_))'],'Color','k','Fontsize',14)
text(0,0,'CONCLUSION:  H_0 is REJECTED  =>  Discriminant Analysis is relevant (significant separation of populations)','Fontsize',14)
%------------------------------------------------------------------------
% Quadratic discriminant function for equal priors and costs
%------------------------------------------------------------------------
figure('Name','Quadratic discriminant function for equal priors and costs')
subplot(4,7,1:3)
axis off
text(-0.1,1.2,'Classification for g=3 MVN populations with equal misclassification costs and equal priors','Fontsize',24)
text(-0.1,0.93,'Pairwise quadratic discriminant functions for new observation x_0:','Fontsize',16)
text(-0.1,0.7,'d_i_k(x_0) = -0.5*x_0^T(S_i^-^1- S_k^-^1)x_0 + (x_i_,_a_v_g^TS_i^-^1- x_k_,_a_v_g^TS_k^-^1)x_0 - K_i_k,  i,k=1,2,3','Fontsize',16)
text(-0.1,0.45,' where  K_i_k = 0.5*log(|S_i|/|S_k|) + 0.5*(x_i_,_a_v_g^TS_i^-^1x_i_,_a_v_g- x_k_,_a_v_g^TS_k^-^1x_k_,_a_v_g)','Fontsize',10)
text(-0.1,0.21,'Decision rule:','Fontsize',16)
text(0.2,0.19,'d_1_2(x_0) > 0 {\wedge} d_1_3(x_0) > 0  {\rightarrow}  d_1','Fontsize',16)
text(0.2,-0.03,'d_2_1(x_0) > 0 {\wedge} d_2_3(x_0) > 0  {\rightarrow}  d_2','Fontsize',16)
text(0.2,-0.25,'d_3_1(x_0) > 0 {\wedge} d_3_2(x_0) > 0  {\rightarrow}  d_3','Fontsize',16)
subplot(4,7,4:7:25)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
annotation('arrow',[0.483125 0.559375],...
    [0.572979591836735 0.573979591836735],'HeadLength',20,'HeadWidth',20,...
    'LineWidth',6);
subplot(4,7,[8:10 15:17 22:24])
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid,axis equal,xlim([-1 8]),ylim([-1 8])
xlabel('x_i_1_j,  (i=1,2,3)  (j=1,...,n_i)','Fontsize',12); ylabel('x_i_2_j','Fontsize',12)
hold on
plot(X2(:,1),X2(:,2),'bo','LineWidth',3)
plot(X3(:,1),X3(:,2),'go','LineWidth',3)
text(-4,7.5,'d_1_2(x) = 0','Color','k','Fontsize',20)
text(-4,6.5,'d_1_3(x) = 0','Color','m','Fontsize',20)
text(-4,5.5,'d_2_3(x) = 0','Color','c','Fontsize',20)
annotation('arrow',[0.13125 0.163125],...
    [0.677571428571428 0.677295918367347],'HeadLength',15,'HeadWidth',15,...
    'LineWidth',6);
annotation('arrow',[0.1325 0.164375],...
    [0.612928571428571 0.612653061224489],'HeadLength',15,'HeadWidth',15,...
    'LineWidth',6,...
    'Color',[1 0 1]);
annotation('arrow',[0.13375 0.165625],...
    [0.545326530612243 0.545051020408162],'HeadLength',15,'HeadWidth',15,...
    'LineWidth',6,...
    'Color',[0 1 1]);
text(-4,1.5,'X_1~N_2({\mu}_1,{\Sigma})','Color','r','Fontsize',16)
text(-4,0.8,'X_2~N_2({\mu}_2,{\Sigma})','Color','b','Fontsize',16)
text(-4,0.1,'X_3~N_2({\mu}_3,{\Sigma})','Color','g','Fontsize',16)
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)
contour(u1,u2,F3,C3,'Color','g','LineWidth',1)
v1_min = min([min(y1),min(z1),min(u1)]);
v1_max = max([max(y1),max(z1),max(u1)]);
v1 = v1_min:(v1_max-v1_min)/100:v1_max;
v2_min = min([min(y2),min(z2),min(u2)]);
v2_max = max([max(y2),max(z2),max(u2)]);
v2 = v2_min:(v2_max-v2_min)/100:v2_max;
[V1,V2] = meshgrid(v1,v2);
V1_vector = V1(:);
V2_vector = V2(:);
K12 = 0.5*log(det(S1)/det(S2)) + 0.5*(x1_avg*inv(S1)*x1_avg' - x2_avg*inv(S2)*x2_avg');
F12_QDA = zeros(1,length(V1_vector));
for q = 1:length(V1_vector)
    F12_QDA(q) = -0.5*[V1(q) V2(q)]*(inv(S1) - inv(S2))*[V1(q); V2(q)] + (x1_avg*inv(S1) - x2_avg*inv(S2))*[V1(q); V2(q)] - K12;
end
F12_QDA = reshape(F12_QDA,length(v2),length(v1));
contour(v1,v2,F12_QDA,[eps eps],'Color','k','LineWidth',3)
K13 = 0.5*log(det(S1)/det(S3)) + 0.5*(x1_avg*inv(S1)*x1_avg' - x3_avg*inv(S3)*x3_avg');
F13_QDA = zeros(1,length(V1_vector));
for q = 1:length(V1_vector)
    F13_QDA(q) = -0.5*[V1(q) V2(q)]*(inv(S1) - inv(S3))*[V1(q); V2(q)] + (x1_avg*inv(S1) - x3_avg*inv(S3))*[V1(q); V2(q)] - K13;
end
F13_QDA = reshape(F13_QDA,length(v2),length(v1));
contour(v1,v2,F13_QDA,[eps eps],'Color','m','LineWidth',3)
K23 = 0.5*log(det(S2)/det(S3)) + 0.5*(x2_avg*inv(S2)*x2_avg' - x3_avg*inv(S3)*x3_avg');
F23_QDA = zeros(1,length(V1_vector));
for q = 1:length(V1_vector)
    F23_QDA(q) = -0.5*[V1(q) V2(q)]*(inv(S2) - inv(S3))*[V1(q); V2(q)] + (x2_avg*inv(S2) - x3_avg*inv(S3))*[V1(q); V2(q)] - K23;
end
F23_QDA = reshape(F23_QDA,length(v2),length(v1));
contour(v1,v2,F23_QDA,[eps eps],'Color','c','LineWidth',3)
hold off
%------------------------------------------------------------------------
% Calculate confusion matrix and APER for training data
%------------------------------------------------------------------------
QDA_12_scores_pop1 = zeros(1,n1);
QDA_13_scores_pop1 = zeros(1,n1);
for q = 1:n1
    QDA_12_scores_pop1(q) = -0.5*X1(q,:)*(inv(S1) - inv(S2))*X1(q,:)' + (x1_avg*inv(S1) - x2_avg*inv(S2))*X1(q,:)' - K12;
    QDA_13_scores_pop1(q) = -0.5*X1(q,:)*(inv(S1) - inv(S3))*X1(q,:)' + (x1_avg*inv(S1) - x3_avg*inv(S3))*X1(q,:)' - K13;
end 
OK_classify_pop1 = find(QDA_12_scores_pop1 > 0 & QDA_13_scores_pop1 > 0);
n1C_APER = length(OK_classify_pop1);
not_OK_classify_pop1 = find(QDA_12_scores_pop1 < 0 | QDA_13_scores_pop1 < 0);
QDA_23_scores_pop1 = zeros(1,length(not_OK_classify_pop1));
for q = 1:length(not_OK_classify_pop1)
    QDA_23_scores_pop1(q) = -0.5*X1(not_OK_classify_pop1(q),:)*(inv(S2) - inv(S3))*X1(not_OK_classify_pop1(q),:)' + (x2_avg*inv(S2) - x3_avg*inv(S3))*X1(not_OK_classify_pop1(q),:)' - K23;
end 
classify_as_2_pop1 = find(QDA_23_scores_pop1 > 0);
n1M2_APER = length(classify_as_2_pop1);
classify_as_3_pop1 = find(QDA_23_scores_pop1 < 0);
n1M3_APER = length(classify_as_3_pop1);
%----------------------------------
QDA_21_scores_pop2 = zeros(1,n2);
QDA_23_scores_pop2 = zeros(1,n2);
for q = 1:n2
    QDA_21_scores_pop2(q) = -0.5*X2(q,:)*(inv(S2) - inv(S1))*X2(q,:)' + (x2_avg*inv(S2) - x1_avg*inv(S1))*X2(q,:)' + K12;
    QDA_23_scores_pop2(q) = -0.5*X2(q,:)*(inv(S2) - inv(S3))*X2(q,:)' + (x2_avg*inv(S2) - x3_avg*inv(S3))*X2(q,:)' - K23;
end 
OK_classify_pop2 = find(QDA_21_scores_pop2 > 0 & QDA_23_scores_pop2 > 0);
n2C_APER = length(OK_classify_pop2);
not_OK_classify_pop2 = find(QDA_21_scores_pop2 < 0 | QDA_23_scores_pop2 < 0);
QDA_13_scores_pop2 = zeros(1,length(not_OK_classify_pop2));
for q = 1:length(not_OK_classify_pop2)
    QDA_13_scores_pop2(q) = -0.5*X2(not_OK_classify_pop2(q),:)*(inv(S1) - inv(S3))*X2(not_OK_classify_pop2(q),:)' + (x1_avg*inv(S1) - x3_avg*inv(S3))*X2(not_OK_classify_pop2(q),:)' - K13;
end 
classify_as_1_pop2 = find(QDA_13_scores_pop2 > 0);
n2M1_APER = length(classify_as_1_pop2);
classify_as_3_pop2 = find(QDA_13_scores_pop2 < 0);
n2M3_APER = length(classify_as_3_pop2);
%----------------------------------
QDA_31_scores_pop3 = zeros(1,n3);
QDA_32_scores_pop3 = zeros(1,n3);
for q = 1:n3
    QDA_31_scores_pop3(q) = -0.5*X3(q,:)*(inv(S3) - inv(S1))*X3(q,:)' + (x3_avg*inv(S3) - x1_avg*inv(S1))*X3(q,:)' + K13;
    QDA_32_scores_pop3(q) = -0.5*X3(q,:)*(inv(S3) - inv(S2))*X3(q,:)' + (x3_avg*inv(S3) - x2_avg*inv(S2))*X3(q,:)' + K23;
end 
OK_classify_pop3 = find(QDA_31_scores_pop3 > 0 & QDA_32_scores_pop3 > 0);
n3C_APER = length(OK_classify_pop3);
not_OK_classify_pop3 = find(QDA_31_scores_pop3 < 0 | QDA_32_scores_pop3 < 0);
QDA_12_scores_pop3 = zeros(1,length(not_OK_classify_pop3));
for q = 1:length(not_OK_classify_pop3)
    QDA_12_scores_pop3(q) = -0.5*X3(not_OK_classify_pop3(q),:)*(inv(S1) - inv(S2))*X3(not_OK_classify_pop3(q),:)' + (x1_avg*inv(S1) - x2_avg*inv(S2))*X3(not_OK_classify_pop3(q),:)' - K12;
end 
classify_as_1_pop3 = find(QDA_12_scores_pop3 > 0);
n3M1_APER = length(classify_as_1_pop3);
classify_as_2_pop3 = find(QDA_12_scores_pop3 < 0);
n3M2_APER = length(classify_as_2_pop3);
%----------------------------------
APER = (n1M2_APER+n1M3_APER+n2M1_APER+n2M3_APER+n3M1_APER+n3M2_APER)/(n1+n2+n3);
%------------------------------------------------------------------------
% Calculate confusion matrix and AER_hat for training data using
% cross-validation (hold-out) procedure
%------------------------------------------------------------------------
n1C_AER = 0;
n1M2_AER = 0;
n1M3_AER = 0;
for q = 1:n1
    X1_red = X1([1:(q-1) q+1:n1],:);
    x1_avg_red = mean(X1_red);
    S1_red = cov(X1_red);
    K12_red = 0.5*log(det(S1_red)/det(S2)) + 0.5*(x1_avg_red*inv(S1_red)*x1_avg_red' - x2_avg*inv(S2)*x2_avg');
    K13_red = 0.5*log(det(S1_red)/det(S3)) + 0.5*(x1_avg_red*inv(S1_red)*x1_avg_red' - x3_avg*inv(S3)*x3_avg');
    QDA_12_score_heldout = -0.5*X1(q,:)*(inv(S1_red) - inv(S2))*X1(q,:)' + (x1_avg_red*inv(S1_red) - x2_avg*inv(S2))*X1(q,:)' - K12_red;
    QDA_13_score_heldout = -0.5*X1(q,:)*(inv(S1_red) - inv(S3))*X1(q,:)' + (x1_avg_red*inv(S1_red) - x3_avg*inv(S3))*X1(q,:)' - K13_red;
    QDA_23_score_heldout = -0.5*X1(q,:)*(inv(S2) - inv(S3))*X1(q,:)' + (x2_avg*inv(S2) - x3_avg*inv(S3))*X1(q,:)' - K23;
    if (QDA_12_score_heldout > 0 && QDA_13_score_heldout > 0)
        n1C_AER = n1C_AER + 1;
    else if QDA_23_score_heldout > 0
        n1M2_AER = n1M2_AER + 1;
        else
            n1M3_AER = n1M3_AER + 1;
        end
    end
end
%----------------------------------
n2C_AER = 0;
n2M1_AER = 0;
n2M3_AER = 0;
for q = 1:n2
    X2_red = X2([1:(q-1) q+1:n1],:);
    x2_avg_red = mean(X2_red);
    S2_red = cov(X2_red);
    K12_red = 0.5*log(det(S1)/det(S2_red)) + 0.5*(x1_avg*inv(S1)*x1_avg' - x2_avg_red*inv(S2_red)*x2_avg_red');
    K23_red = 0.5*log(det(S2_red)/det(S3)) + 0.5*(x2_avg_red*inv(S2_red)*x2_avg_red' - x3_avg*inv(S3)*x3_avg');
    QDA_21_score_heldout = -0.5*X2(q,:)*(inv(S2_red) - inv(S1))*X2(q,:)' + (x2_avg_red*inv(S2_red) - x1_avg*inv(S1))*X2(q,:)' + K12_red;
    QDA_23_score_heldout = -0.5*X2(q,:)*(inv(S2_red) - inv(S3))*X2(q,:)' + (x2_avg_red*inv(S2_red) - x3_avg*inv(S3))*X2(q,:)' - K23_red;
    QDA_13_score_heldout = -0.5*X2(q,:)*(inv(S1) - inv(S3))*X2(q,:)' + (x1_avg*inv(S1) - x3_avg*inv(S3))*X2(q,:)' - K13;
    if (QDA_21_score_heldout > 0 && QDA_23_score_heldout > 0)
        n2C_AER = n2C_AER + 1;
    else if QDA_13_score_heldout > 0
        n2M1_AER = n2M1_AER + 1;
        else
            n2M3_AER = n2M3_AER + 1;
        end
    end
end
%----------------------------------
n3C_AER = 0;
n3M1_AER = 0;
n3M2_AER = 0;
for q = 1:n3
    X3_red = X3([1:(q-1) q+1:n1],:);
    x3_avg_red = mean(X3_red);
    S3_red = cov(X3_red);
    K13_red = 0.5*log(det(S1)/det(S3_red)) + 0.5*(x1_avg*inv(S1)*x1_avg' - x3_avg_red*inv(S3_red)*x3_avg_red');
    K23_red = 0.5*log(det(S2)/det(S3_red)) + 0.5*(x2_avg*inv(S2)*x2_avg' - x3_avg_red*inv(S3_red)*x3_avg_red');
    QDA_31_score_heldout = -0.5*X3(q,:)*(inv(S3_red) - inv(S1))*X3(q,:)' + (x3_avg_red*inv(S3_red) - x1_avg*inv(S1))*X3(q,:)' + K13_red;
    QDA_32_score_heldout = -0.5*X3(q,:)*(inv(S3_red) - inv(S2))*X3(q,:)' + (x3_avg_red*inv(S3_red) - x2_avg*inv(S2))*X3(q,:)' + K23_red;
    QDA_12_score_heldout = -0.5*X3(q,:)*(inv(S1) - inv(S2))*X3(q,:)' + (x1_avg*inv(S1) - x2_avg*inv(S2))*X3(q,:)' - K12;
    if (QDA_31_score_heldout > 0 && QDA_32_score_heldout > 0)
        n3C_AER = n3C_AER + 1;
    else if QDA_12_score_heldout > 0
        n3M1_AER = n3M1_AER + 1;
        else
            n3M2_AER = n3M2_AER + 1;
        end
    end
end
%----------------------------------
AER_hat = (n1M2_AER+n1M3_AER+n2M1_AER+n2M3_AER+n3M1_AER+n3M2_AER)/(n1+n2+n3);
%------------------------------------------------------------------------
% Printout confusion tables and error rates
%------------------------------------------------------------------------
subplot(4,7,[5:7 12:14 19:21 26:28])
axis off
text(0,0.95,'Performance of classification','Fontsize',22)
text(0,0.91,'(empirical error rates)','Fontsize',22)
text(0,0.85,'------------------------------------------','Fontsize',20)
text(0,0.79,'Confusion matrix for training data:','Fontsize',18)
text(0.37,0.72,'classified as','Fontsize',16)
text(0.4,0.67,'d_1    d_2    d_3','Fontsize',16)
text(0.1,0.59,'true','Fontsize',16)
text(0.05,0.56,'population','Fontsize',16)
text(0.32,0.62,'{\pi}_1','Fontsize',16)
text(0.4,0.62,num2str(n1C_APER),'Fontsize',16)
text(0.5,0.62,num2str(n1M2_APER),'Fontsize',16)
text(0.59,0.62,num2str(n1M3_APER),'Fontsize',16)
text(0.68,0.62,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',14)
text(0.32,0.58,'{\pi}_2','Fontsize',16)
text(0.41,0.58,num2str(n2M1_APER),'Fontsize',16)
text(0.49,0.58,num2str(n2C_APER),'Fontsize',16)
text(0.59,0.58,num2str(n2M3_APER),'Fontsize',16)
text(0.68,0.58,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',14)
text(0.32,0.54,'{\pi}_3','Fontsize',16)
text(0.41,0.54,num2str(n3M1_APER),'Fontsize',16)
text(0.49,0.54,num2str(n3M2_APER),'Fontsize',16)
text(0.58,0.54,num2str(n3C_APER),'Fontsize',16)
text(0.68,0.54,['(sum = n_3 = ' num2str(n3) ')'],'Fontsize',14)
text(0,0.47,['Apparent Error Rate, APER = ' num2str(n1M2_APER+n1M3_APER+n2M1_APER+n2M3_APER+n3M1_APER+n3M2_APER) '/' num2str(n1+n1+n3) ' = ' num2str(APER,3)],'Fontsize',16)
text(0,0.42,'------------------------------------------','Fontsize',20)
text(0,0.36,'Confusion matrix for training data:','Fontsize',18)
text(0,0.32,'using cross-validation (hold-out):','Fontsize',18)
text(0.37,0.25,'classified as','Fontsize',16)
text(0.4,0.2,'d_1    d_2    d_3','Fontsize',16)
text(0.1,0.12,'true','Fontsize',16)
text(0.05,0.09,'population','Fontsize',16)
text(0.32,0.15,'{\pi}_1','Fontsize',16)
text(0.4,0.15,num2str(n1C_AER),'Fontsize',16)
text(0.5,0.15,num2str(n1M2_AER),'Fontsize',16)
text(0.59,0.15,num2str(n1M3_AER),'Fontsize',16)
text(0.68,0.15,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',14)
text(0.32,0.11,'{\pi}_2','Fontsize',16)
text(0.41,0.11,num2str(n2M1_AER),'Fontsize',16)
text(0.49,0.11,num2str(n2C_AER),'Fontsize',16)
text(0.59,0.11,num2str(n2M3_AER),'Fontsize',16)
text(0.68,0.11,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',14)
text(0.32,0.07,'{\pi}_3','Fontsize',16)
text(0.41,0.07,num2str(n3M1_AER),'Fontsize',16)
text(0.49,0.07,num2str(n3M2_AER),'Fontsize',16)
text(0.58,0.07,num2str(n3C_AER),'Fontsize',16)
text(0.68,0.07,['(sum = n_3 = ' num2str(n3) ')'],'Fontsize',14)
text(0,0,['Estimated Actual Error Rate, AER_h_a_t = ' num2str(n1M2_AER+n1M3_AER+n2M1_AER+n2M3_AER+n3M1_AER+n3M2_AER) '/' num2str(n1+n1+n3) ' = ' num2str(AER_hat,3)],'Fontsize',16)
text(0,-0.05,'------------------------------------------','Fontsize',20)
%------------------------------------------------------------------------
% Comparison with (un-appropriate) linear discriminant function
%------------------------------------------------------------------------
figure('Name','Linear discriminant function for equal priors and costs')
subplot(4,7,1:3)
axis off
text(0,1.2,'Classification for g=3 MVN populations with equal misclassification costs and equal priors','Fontsize',24)
text(0,0.8,'Pairwise linear discriminant functions for new observation x_0:','Fontsize',16)
text(0,0.6,'d_i_k(x_0) = (x_i_,_a_v_g- x_k_,_a_v_g)^T S_p^-^1 [x_0 - 0.5(x_i_,_a_v_g+ x_k_,_a_v_g)],  i,k=1,2,3','Fontsize',16)
text(0,0.27,'Decision rule:','Fontsize',16)
text(0.3,0.25,'d_1_2(x_0) > 0 {\wedge} d_1_3(x_0) > 0  {\rightarrow}  d_1','Fontsize',16)
text(0.3,0.03,'d_2_1(x_0) > 0 {\wedge} d_2_3(x_0) > 0  {\rightarrow}  d_2','Fontsize',16)
text(0.3,-0.19,'d_3_1(x_0) > 0 {\wedge} d_3_2(x_0) > 0  {\rightarrow}  d_3','Fontsize',16)
subplot(4,7,4:7:25)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
annotation('arrow',[0.483125 0.559375],...
    [0.572979591836735 0.573979591836735],'HeadLength',20,'HeadWidth',20,...
    'LineWidth',6);
subplot(4,7,[8:10 15:17 22:24])
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid,axis equal,xlim([-1 8]),ylim([-1 8])
xlabel('x_i_1_j,  (i=1,2,3)  (j=1,...,n_i)','Fontsize',12); ylabel('x_i_2_j','Fontsize',12)
hold on
plot(X2(:,1),X2(:,2),'bo','LineWidth',3)
plot(X3(:,1),X3(:,2),'go','LineWidth',3)
text(-4,7.5,'d_1_2(x) = 0','Color','k','Fontsize',20)
text(-4,6.5,'d_1_3(x) = 0','Color','m','Fontsize',20)
text(-4,5.5,'d_2_3(x) = 0','Color','c','Fontsize',20)
annotation('arrow',[0.13125 0.163125],...
    [0.677571428571428 0.677295918367347],'HeadLength',15,'HeadWidth',15,...
    'LineWidth',6);
annotation('arrow',[0.1325 0.164375],...
    [0.612928571428571 0.612653061224489],'HeadLength',15,'HeadWidth',15,...
    'LineWidth',6,...
    'Color',[1 0 1]);
annotation('arrow',[0.13375 0.165625],...
    [0.545326530612243 0.545051020408162],'HeadLength',15,'HeadWidth',15,...
    'LineWidth',6,...
    'Color',[0 1 1]);
text(-4,1.5,'X_1~N_2({\mu}_1,{\Sigma})','Color','r','Fontsize',16)
text(-4,0.8,'X_2~N_2({\mu}_2,{\Sigma})','Color','b','Fontsize',16)
text(-4,0.1,'X_3~N_2({\mu}_3,{\Sigma})','Color','g','Fontsize',16)
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)
contour(u1,u2,F3,C3,'Color','g','LineWidth',1)
v1_min = min([min(y1),min(z1),min(u1)]);
v1_max = max([max(y1),max(z1),max(u1)]);
v1 = v1_min:(v1_max-v1_min)/100:v1_max;
v2_min = min([min(y2),min(z2),min(u2)]);
v2_max = max([max(y2),max(z2),max(u2)]);
v2 = v2_min:(v2_max-v2_min)/100:v2_max;
[V1,V2] = meshgrid(v1,v2);
F12 = (x1_avg-x2_avg)*inv(Sp)*[V1(:)' ; V2(:)'] - 0.5*(x1_avg-x2_avg)*inv(Sp)*(x1_avg+x2_avg)';
F12 = reshape(F12,length(v2),length(v1));
contour(v1,v2,F12,[eps eps],'Color','k','LineWidth',3)
F13 = (x1_avg-x3_avg)*inv(Sp)*[V1(:)' ; V2(:)'] - 0.5*(x1_avg-x3_avg)*inv(Sp)*(x1_avg+x3_avg)';
F13 = reshape(F13,length(v2),length(v1));
contour(v1,v2,F13,[eps eps],'Color','m','LineWidth',3)
F23 = (x2_avg-x3_avg)*inv(Sp)*[V1(:)' ; V2(:)'] - 0.5*(x2_avg-x3_avg)*inv(Sp)*(x2_avg+x3_avg)';
F23 = reshape(F23,length(v2),length(v1));
contour(v1,v2,F23,[eps eps],'Color','c','LineWidth',3)
hold off
%------------------------------------------------------------------------
% Calculate confusion matrix and APER for training data
%------------------------------------------------------------------------
LDA_12_scores_pop1 = (x1_avg-x2_avg)*inv(Sp)*X1' - 0.5*(x1_avg-x2_avg)*inv(Sp)*(x1_avg+x2_avg)';
LDA_13_scores_pop1 = (x1_avg-x3_avg)*inv(Sp)*X1' - 0.5*(x1_avg-x3_avg)*inv(Sp)*(x1_avg+x3_avg)';
OK_classify_pop1 = find(LDA_12_scores_pop1 > 0 & LDA_13_scores_pop1 > 0);
n1C_APER = length(OK_classify_pop1);
not_OK_classify_pop1 = find(LDA_12_scores_pop1 < 0 | LDA_13_scores_pop1 < 0);
LDA_23_scores_pop1 = (x2_avg-x3_avg)*inv(Sp)*X1(not_OK_classify_pop1,:)' - 0.5*(x2_avg-x3_avg)*inv(Sp)*(x2_avg+x3_avg)';
classify_as_2_pop1 = find(LDA_23_scores_pop1 > 0);
n1M2_APER = length(classify_as_2_pop1);
classify_as_3_pop1 = find(LDA_23_scores_pop1 < 0);
n1M3_APER = length(classify_as_3_pop1);
%----------------------------------
LDA_21_scores_pop2 = (x2_avg-x1_avg)*inv(Sp)*X2' - 0.5*(x2_avg-x1_avg)*inv(Sp)*(x2_avg+x1_avg)';
LDA_23_scores_pop2 = (x2_avg-x3_avg)*inv(Sp)*X2' - 0.5*(x2_avg-x3_avg)*inv(Sp)*(x2_avg+x3_avg)';
OK_classify_pop2 = find(LDA_21_scores_pop2 > 0 & LDA_23_scores_pop2 > 0);
n2C_APER = length(OK_classify_pop2);
not_OK_classify_pop2 = find(LDA_21_scores_pop2 < 0 | LDA_23_scores_pop2 < 0);
LDA_13_scores_pop2 = (x1_avg-x3_avg)*inv(Sp)*X2(not_OK_classify_pop2,:)' - 0.5*(x1_avg-x3_avg)*inv(Sp)*(x1_avg+x3_avg)';
classify_as_1_pop2 = find(LDA_13_scores_pop2 > 0);
n2M1_APER = length(classify_as_1_pop2);
classify_as_3_pop2 = find(LDA_13_scores_pop2 < 0);
n2M3_APER = length(classify_as_3_pop2);
%----------------------------------
LDA_31_scores_pop3 = (x3_avg-x1_avg)*inv(Sp)*X3' - 0.5*(x3_avg-x1_avg)*inv(Sp)*(x3_avg+x1_avg)';
LDA_32_scores_pop3 = (x3_avg-x2_avg)*inv(Sp)*X3' - 0.5*(x3_avg-x2_avg)*inv(Sp)*(x3_avg+x2_avg)';
OK_classify_pop3 = find(LDA_31_scores_pop3 > 0 & LDA_32_scores_pop3 > 0);
n3C_APER = length(OK_classify_pop3);
not_OK_classify_pop3 = find(LDA_31_scores_pop3 < 0 | LDA_32_scores_pop3 < 0);
LDA_12_scores_pop3 = (x1_avg-x2_avg)*inv(Sp)*X3(not_OK_classify_pop3,:)' - 0.5*(x1_avg-x2_avg)*inv(Sp)*(x1_avg+x2_avg)';
classify_as_1_pop3 = find(LDA_12_scores_pop3 > 0);
n3M1_APER = length(classify_as_1_pop3);
classify_as_2_pop3 = find(LDA_12_scores_pop3 < 0);
n3M2_APER = length(classify_as_2_pop3);
%----------------------------------
APER = (n1M2_APER+n1M3_APER+n2M1_APER+n2M3_APER+n3M1_APER+n3M2_APER)/(n1+n2+n3);
%------------------------------------------------------------------------
% Calculate confusion matrix and AER_hat for training data using
% cross-validation (hold-out) procedure
%------------------------------------------------------------------------
n1C_AER = 0;
n1M2_AER = 0;
n1M3_AER = 0;
for j1 = 1:n1
    X1_red = X1([1:(j1-1) j1+1:n1],:);
    x1_avg_red = mean(X1_red);
    S1_red = cov(X1_red);
    Sp_red = ((n1-2)*S1_red + (n2-1)*S2 + (n3-1)*S3)/(n1+n2+n3-4);
    LDA_12_score_heldout = (x1_avg_red-x2_avg)*inv(Sp_red)*X1(j1,:)' - 0.5*(x1_avg_red-x2_avg)*inv(Sp_red)*(x1_avg_red+x2_avg)';
    LDA_13_score_heldout = (x1_avg_red-x3_avg)*inv(Sp_red)*X1(j1,:)' - 0.5*(x1_avg_red-x3_avg)*inv(Sp_red)*(x1_avg_red+x3_avg)';
    LDA_23_score_heldout = (x2_avg-x3_avg)*inv(Sp_red)*X1(j1,:)' - 0.5*(x2_avg-x3_avg)*inv(Sp_red)*(x2_avg+x3_avg)';
    if (LDA_12_score_heldout > 0 && LDA_13_score_heldout > 0)
        n1C_AER = n1C_AER + 1;
    else if LDA_23_score_heldout > 0    
        n1M2_AER = n1M2_AER + 1;
        else
            n1M3_AER = n1M3_AER + 1;
        end
    end
end
%----------------------------------
n2C_AER = 0;
n2M1_AER = 0;
n2M3_AER = 0;
for j2 = 1:n2
    X2_red = X2([1:(j2-1) j2+1:n2],:);
    x2_avg_red = mean(X2_red);
    S2_red = cov(X2_red);
    Sp_red = ((n1-1)*S1 + (n2-2)*S2_red+ (n3-1)*S3)/(n1+n2+n3-4);
    LDA_21_score_heldout = (x2_avg_red-x1_avg)*inv(Sp_red)*X2(j2,:)' - 0.5*(x2_avg_red-x1_avg)*inv(Sp_red)*(x2_avg_red+x1_avg)';
    LDA_23_score_heldout = (x2_avg_red-x3_avg)*inv(Sp_red)*X2(j2,:)' - 0.5*(x2_avg_red-x3_avg)*inv(Sp_red)*(x2_avg_red+x3_avg)';
    LDA_13_score_heldout = (x1_avg-x3_avg)*inv(Sp_red)*X2(j2,:)' - 0.5*(x1_avg-x3_avg)*inv(Sp_red)*(x1_avg+x3_avg)';
    if (LDA_21_score_heldout > 0 && LDA_23_score_heldout > 0)
        n2C_AER = n2C_AER + 1;
    else if LDA_13_score_heldout > 0    
        n2M1_AER = n2M1_AER + 1;
        else
            n2M3_AER = n2M3_AER + 1;
        end
    end
end
%----------------------------------
n3C_AER = 0;
n3M1_AER = 0;
n3M2_AER = 0;
for j3 = 1:n3
    X3_red = X3([1:(j3-1) j3+1:n3],:);
    x3_avg_red = mean(X3_red);
    S3_red = cov(X3_red);
    Sp_red = ((n1-1)*S1 + (n2-1)*S2+ (n3-2)*S3_red)/(n1+n2+n3-4);
    LDA_31_score_heldout = (x3_avg_red-x1_avg)*inv(Sp_red)*X3(j3,:)' - 0.5*(x3_avg_red-x1_avg)*inv(Sp_red)*(x3_avg_red+x1_avg)';
    LDA_32_score_heldout = (x3_avg_red-x2_avg)*inv(Sp_red)*X3(j3,:)' - 0.5*(x3_avg_red-x2_avg)*inv(Sp_red)*(x3_avg_red+x2_avg)';
    LDA_12_score_heldout = (x1_avg-x2_avg)*inv(Sp_red)*X3(j3,:)' - 0.5*(x1_avg-x2_avg)*inv(Sp_red)*(x1_avg+x2_avg)';
    if (LDA_31_score_heldout > 0 && LDA_32_score_heldout > 0)
        n3C_AER = n3C_AER + 1;
    else if LDA_12_score_heldout > 0    
        n3M1_AER = n3M1_AER + 1;
        else
            n3M2_AER = n3M2_AER + 1;
        end
    end
end
%----------------------------------
AER_hat = (n1M2_AER+n1M3_AER+n2M1_AER+n2M3_AER+n3M1_AER+n3M2_AER)/(n1+n2+n3);
%------------------------------------------------------------------------
% Printout confusion tables and error rates
%------------------------------------------------------------------------
subplot(4,7,[5:7 12:14 19:21 26:28])
axis off
text(0,0.95,'Performance of classification','Fontsize',22)
text(0,0.91,'(empirical error rates)','Fontsize',22)
text(0,0.85,'------------------------------------------','Fontsize',20)
text(0,0.79,'Confusion matrix for training data:','Fontsize',18)
text(0.37,0.72,'classified as','Fontsize',16)
text(0.4,0.67,'d_1    d_2    d_3','Fontsize',16)
text(0.1,0.59,'true','Fontsize',16)
text(0.05,0.56,'population','Fontsize',16)
text(0.32,0.62,'{\pi}_1','Fontsize',16)
text(0.4,0.62,num2str(n1C_APER),'Fontsize',16)
text(0.5,0.62,num2str(n1M2_APER),'Fontsize',16)
text(0.59,0.62,num2str(n1M3_APER),'Fontsize',16)
text(0.68,0.62,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',14)
text(0.32,0.58,'{\pi}_2','Fontsize',16)
text(0.41,0.58,num2str(n2M1_APER),'Fontsize',16)
text(0.49,0.58,num2str(n2C_APER),'Fontsize',16)
text(0.59,0.58,num2str(n2M3_APER),'Fontsize',16)
text(0.68,0.58,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',14)
text(0.32,0.54,'{\pi}_3','Fontsize',16)
text(0.41,0.54,num2str(n3M1_APER),'Fontsize',16)
text(0.49,0.54,num2str(n3M2_APER),'Fontsize',16)
text(0.58,0.54,num2str(n3C_APER),'Fontsize',16)
text(0.68,0.54,['(sum = n_3 = ' num2str(n3) ')'],'Fontsize',14)
text(0,0.47,['Apparent Error Rate, APER = ' num2str(n1M2_APER+n1M3_APER+n2M1_APER+n2M3_APER+n3M1_APER+n3M2_APER) '/' num2str(n1+n1+n3) ' = ' num2str(APER,3)],'Fontsize',16)
text(0,0.42,'------------------------------------------','Fontsize',20)
text(0,0.36,'Confusion matrix for training data:','Fontsize',18)
text(0,0.32,'using cross-validation (hold-out):','Fontsize',18)
text(0.37,0.25,'classified as','Fontsize',16)
text(0.4,0.2,'d_1    d_2    d_3','Fontsize',16)
text(0.1,0.12,'true','Fontsize',16)
text(0.05,0.09,'population','Fontsize',16)
text(0.32,0.15,'{\pi}_1','Fontsize',16)
text(0.4,0.15,num2str(n1C_AER),'Fontsize',16)
text(0.5,0.15,num2str(n1M2_AER),'Fontsize',16)
text(0.59,0.15,num2str(n1M3_AER),'Fontsize',16)
text(0.68,0.15,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',14)
text(0.32,0.11,'{\pi}_2','Fontsize',16)
text(0.41,0.11,num2str(n2M1_AER),'Fontsize',16)
text(0.49,0.11,num2str(n2C_AER),'Fontsize',16)
text(0.59,0.11,num2str(n2M3_AER),'Fontsize',16)
text(0.68,0.11,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',14)
text(0.32,0.07,'{\pi}_3','Fontsize',16)
text(0.41,0.07,num2str(n3M1_AER),'Fontsize',16)
text(0.49,0.07,num2str(n3M2_AER),'Fontsize',16)
text(0.58,0.07,num2str(n3C_AER),'Fontsize',16)
text(0.68,0.07,['(sum = n_3 = ' num2str(n3) ')'],'Fontsize',14)
text(0,0,['Estimated Actual Error Rate, AER_h_a_t = ' num2str(n1M2_AER+n1M3_AER+n2M1_AER+n2M3_AER+n3M1_AER+n3M2_AER) '/' num2str(n1+n1+n3) ' = ' num2str(AER_hat,3)],'Fontsize',16)
text(0,-0.05,'------------------------------------------','Fontsize',20)
% %-------------------------------------------------------------------------------------------------------