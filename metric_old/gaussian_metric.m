function x=gaussian_metric(y1,y2,Sigma,H1,H2,sigma)
L=inv(inv(Sigma)+H1'*H1/sigma^2);
M=inv(inv(Sigma)+H2'*H2/sigma^2);
alpha=-L*(H1'*y1)/sigma^2;
beta=-M*(H2'*y2)/sigma^2;
x=(alpha-beta)'*inv(L+M)*(alpha-beta)+log(det(L+M));