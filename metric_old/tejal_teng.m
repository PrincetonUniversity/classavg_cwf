clear all
for ii=1:100
n=11;
d=5; % Number of clean images
x=zeros(n,d);

%Draw signal from sinusoids


% Draw signal from a standard gaussian\
U=orth(randn(n));
%V=[10,10,10,10,10,10,1,1,1,1,1];
%V=[10,10,1,1,1,1,1,1,1,1,1];
V=ones(1,n);
Sigma=U*diag(V)*U'; % SVD

for i=1:d
    x(:,i)=(U*diag(V.^0.5))*randn(n,1);
end
Xc=zeros(n,d);
Y=zeros(n,d);
H=zeros(n,d);


for i=1:d
    for k= i*10-9:i*10
        Xc(:,k)=x(:,i); % 10 copies of each clean signal
    end
end
it=1;
sigma =0.3;

H1d=zeros(n,n,d*10);
for i=1:d*10
    [h, h1d]=CTF(n,3,0.02, 1 + 3.*rand(1,1),2,100, 0.07);  %CTF function from Jane
    H1d(:,:,i)=diag(h1d);
	H1d(:,:,i)=diag(rand(1,n)*2-1);
    N= (sigma^2).*randn(n,1); % mean 0, standard deviation sigma
    Y(:,i)=H1d(:,:,i)*Xc(:,i)+N;
	signal(:,i)=H1d(:,:,i)*Xc(:,i);
end
snr=norm(signal,'fro')/sigma/(d*10*n)^0.5;
meansnr(it)=snr;

affinity=zeros(d*10,d*10);affinity1=zeros(d*10,d*10);
affinity2=zeros(d*10,d*10);
for i=1:d*10
    for k=1:d*10
affinity(i,k)=gaussian_metric(Y(:,i),Y(:,k),Sigma,H1d(:,:,i),H1d(:,:,k),sigma);
affinity1(i,k)=gaussian_metric(Y(:,i),Y(:,k),eye(n),H1d(:,:,i),H1d(:,:,k),sigma);
affinity2(i,k)=gaussian_metric(Y(:,i),Y(:,k),eye(n)*10,H1d(:,:,i),H1d(:,:,k),sigma);
	end
end

Yc=Y;
% Phase flipping
for i=1:d*10
     Yc(:,i)=sign(H1d(:,:,i))*Y(:,i);
end

crosscorr=zeros(d*10,d*10);
for i=1:d*10
    for k=1:d*10
       crosscorr(i,k)=dot(Yc(:,i),Yc(:,k))/(norm(Yc(:,i))*norm(Yc(:,k)));
	   crosscorr1(i,k)=dot(Yc(:,i)-mean(Yc,2),Yc(:,k)-mean(Yc,2))/(norm(Yc(:,i)-mean(Yc,2))*norm(Yc(:,k)-mean(Yc,2)));
    end
end

temp=repmat(1:d,1,10);
true_neighbor=zeros(10*d,10*d);
for i=1:10*d
for j=1:10*d
true_neighbor(i,j)=(floor((i-1)/10)==floor((j-1)/10));
end
end
error(1)=find_neighbors(-affinity,true_neighbor,10);
error(2)=find_neighbors(-affinity1,true_neighbor,10);
error(3)=find_neighbors(-affinity2,true_neighbor,10);
error(4)=find_neighbors(crosscorr,true_neighbor,10);
error(5)=find_neighbors(crosscorr1,true_neighbor,10);
error
allerror(ii,:)=error;
end