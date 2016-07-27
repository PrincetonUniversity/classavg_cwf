clear all

noise_vals=300; % Number of values of noise variance sigma for which the code calculates metric
false_cross=zeros(noise_vals,1);
false_aff=zeros(noise_vals,1);


n=11;
d=5; % Number of clean images
j=sqrt(-1);

x=zeros(n,d);

%Draw signal from sinusoids


% Draw signal from a standard gaussian
for i=1:d
    x(:,i)=randn(n,1);
end

%F=dftmtx(n); %DFT matrix add this after for modified code with complex
%normal distribution

F=eye(n);

X=F*x;  %Columns give the signal in Fourier space

Xc=zeros(n,d);
for i=1:d
    for k= i*10-9:i*10
        Xc(:,k)=X(:,i); % 10 copies of each clean signal
    end
end

Y=zeros(n,d);
H=zeros(n,d);
n_level=0.1:0.005:3;


mu=zeros(n,1);
U=RandOrthMat(n);
V=[10,10,10,10,10,10,1,1,1,1,1];
Sigma=U*diag(V)*U'; % SVD
H1d=zeros(n,n,d*10);
K=zeros(n,n,d*10);
Filter=zeros(n,n);

wrong_Sigma1=eye(n);
wrong_Sigma2=10*eye(n);

for it=1:noise_vals
sigma =n_level(it);
for i=1:d*10
    [h, h1d]=CTF(n,3,0.02, 1 + 3.*rand(1,1),2,100, 0.07);  %CTF function from Jane
    H1d(:,:,i)=diag(h1d);
    N= (sigma^2).*randn(n,1); % mean 0, standard deviation sigma
    Y(:,i)=H1d(:,:,i)*Xc(:,i)+N;
    K(:,:,i)=H1d(:,:,i)*Sigma*H1d(:,:,i)' + sigma^2*eye(n,n);
    signal=H1d(:,:,i)*Xc(:,i);
    Filter=Filter+(H1d(:,:,i).*H1d(:,:,i));
end

half=0.5*eye(n,n);
Filter=power(Filter,half);
snr=mean(signal.*signal)/sigma; % signal: excluding noise
meansnr(it)=sum(snr)/(d*10);
%% Compare 100C2 images generated wrt the metric

Nbs=zeros(d*10,d*10);
affinity=zeros(d*10,d*10);

t1=tic;
for i=1:d*10
    parfor k=1:d*10
        L=wrong_Sigma1 - wrong_Sigma1*H1d(:,:,i)'*(H1d(:,:,i)*wrong_Sigma1*H1d(:,:,i)' + sigma^2*eye(n))^(-1)*H1d(:,:,i)*wrong_Sigma1;
        M=wrong_Sigma1 - wrong_Sigma1*H1d(:,:,k)'*(H1d(:,:,k)*wrong_Sigma1*H1d(:,:,k)' + sigma^2*eye(n))^(-1)*H1d(:,:,k)*wrong_Sigma1;
        alpha= mu + wrong_Sigma1*H1d(:,:,i)'*(H1d(:,:,i)*wrong_Sigma1*H1d(:,:,i)'+ sigma^2*eye(n))^(-1)*(Y(:,i)-H1d(:,:,i)*mu);
        beta=  mu + wrong_Sigma1*H1d(:,:,k)'*(H1d(:,:,k)*wrong_Sigma1*H1d(:,:,k)'+ sigma^2*eye(n))^(-1)*(Y(:,k)-H1d(:,:,k)*mu);
       
        
        affinity(i,k)= -0.5*log(det(Filter*(L+M)*Filter')) -0.5*((Filter*(alpha-beta))'*((Filter*(L+M)*Filter')\(Filter*(alpha-beta))));
        
        
    end
end
t1=toc(t1);    

%% Compare with normalized cross correlation
Yc=Y;
% Phase flipping
for i=1:d*10
     Yc(:,i)=sign(H1d(:,:,i))*Y(:,i);
end

crosscorr=zeros(d*10,d*10);
t2=tic;
for i=1:d*10
    for k=1:d*10
       crosscorr(i,k)=dot(Yc(:,i)-mean(Yc,2),Yc(:,k)-mean(Yc,2))/(norm(Yc(:,i)-mean(Yc,2))*norm(Yc(:,k)-mean(Yc,2)));
    end
end
t2=toc(t2);

%% Statistics of identification of neighbors

%%%%%Cross correlation
Sort_aff=zeros(d*10,d*10);
dim_aff=zeros(d*10,d*10);
for i=1:d*10
[Sort_aff(i,:),dim_aff(i,:)]=sort(affinity(i,:));
end
Sort_aff=flipdim(Sort_aff, 2);
dim_aff=flipdim(dim_aff,2); % Sorted rowwise in descending order

%%%%%Cross correlation
Sort_cross=zeros(d*10,d*10);
dim_cross=zeros(d*10,d*10);
for i=1:d*10
[Sort_cross(i,:),dim_cross(i,:)]=sort(crosscorr(i,:));
end
Sort_cross=flipdim(Sort_cross, 2);
dim_cross=flipdim(dim_cross,2); % Sorted rowwise in descending order


t3=tic;
for p=1:10:d*10
for i=p:p+9
    for k=1:10
        if dim_aff(i,k) > p+9
            false_aff(it)=false_aff(it)+1;
        end
        if dim_cross(i,k) > p+9
            false_cross(it)=false_cross(it)+1;
        end
    end
end
end
t3=toc(t3);
%imagesc(log(affinity))
end



    