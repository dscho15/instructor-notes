function cr_ellipsis (mu_hat, SIGMA_hat, alpha, n);
t = 0:0.01:2*pi;
N = length(t);
p = 2;
[V,D] = eig(SIGMA_hat);
a = sqrt(D(2,2))*sqrt((p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p));
b = sqrt(D(1,1))*sqrt((p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p));
P  = [a*cos(t); b*sin(t)];
theta = atan2(V(2,2),V(1,2));
T = [cos(theta) -sin(theta); 
     sin(theta) cos(theta)];
P_rot = T*P + mu_hat*ones(1,N);
plot(P_rot(1,:),P_rot(2,:),'LineWidth',3,'Color','k', 'DisplayName','Confidence region'),grid