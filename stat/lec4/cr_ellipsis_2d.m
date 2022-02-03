function cr_ellipsis_2d(mu_diff_hat, SIGMA_pooled_hat, alpha, n1, n2)
t = 0:0.01:2*pi;
N = length(t);
p = 2;
[V,D] = eig(SIGMA_pooled_hat);
a = sqrt(D(2,2))*sqrt((p*(n1+n2-2)*(1/n1+1/n2)/(n1+n2-p-1))*finv(1-alpha,p,n1+n2-p-1));
b = sqrt(D(1,1))*sqrt((p*(n1+n2-2)*(1/n1+1/n2)/(n1+n2-p-1))*finv(1-alpha,p,n1+n2-p-1));
P  = [a*cos(t); b*sin(t)];
theta = atan2(V(2,2),V(1,2));
T = [cos(theta) -sin(theta); sin(theta) cos(theta)];
P_rot = T*P + mu_diff_hat*ones(1,N);
plot(P_rot(1,:),P_rot(2,:),'LineWidth',3,'Color','k', 'DisplayName','Confidence region'),grid