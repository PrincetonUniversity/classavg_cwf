clear all

noise_vals=1; % Number of values of noise variance sigma for which the code calculates metric
true_cross=zeros(noise_vals,1);
true_aff=zeros(noise_vals,1);
true_w=zeros(noise_vals,1);

for ii=1:20
    n=3;
    d=8; % Number of clean images
    %j=sqrt(-1);
    
    x=zeros(n,d);
    
    
    U=orth(randn(n));
    V=[10,1,1];
    Sigma=U*diag(V)*U'; % SVD
    
    % Draw signal from a standard gaussian
    for i=1:d
        X(:,i)=(U*diag(V.^0.5))*randn(n,1);
    end
    
    
    Xc=zeros(n,d);
    for i=1:d
        for k= i*10-9:i*10
            Xc(:,k)=X(:,i); % 10 copies of each clean signal
        end
    end
    
    Y=zeros(n,d);
    H=zeros(n,d);
    n_level=0.1;%:0.005:3;
    
    
    mu=zeros(n,1);
    
    H1d=zeros(n,n,d*10);
    K=zeros(n,n,d*10);
    
    
    %wrong_Sigma1=eye(n);
    %wrong_Sigma1=10*eye(n);
    wrong_Sigma1=Sigma;
    
    for it=1:length(n_level)
        sigma =n_level(it);
        for i=1:d*10
            [h, h1d]=CTF(n,3,0.02, 1 + 3.*rand(1,1),2,100, 0.07);  %CTF function from Jane
            H1d(:,:,i)=diag(h1d);
            N= (sigma^2).*randn(n,1); % mean 0, standard deviation sigma
            Y(:,i)=H1d(:,:,i)*Xc(:,i)+N;
            K(:,:,i)=H1d(:,:,i)*Sigma*H1d(:,:,i)' + sigma^2*eye(n,n);
            signal=H1d(:,:,i)*Xc(:,i);
            
        end
        
        snr=mean(signal.*signal)/sigma; % signal: excluding noise
        meansnr(it)=sum(snr)/(d*10);
        %% Compare 100C2 images generated wrt the metric
        
        Nbs=zeros(d*10,d*10);
        affinity=zeros(d*10,d*10);
        affinity2=zeros(d*10,d*10);
        t1=tic;
        for i=1:d*10
            for k=1:d*10
                L=inv(inv(Sigma)+H1d(:,:,i)'*H1d(:,:,i)/sigma^2);
                M=inv(inv(Sigma)+H1d(:,:,k)'*H1d(:,:,k)/sigma^2);
                alpha=-L*(H1d(:,:,i)'*Y(:,i))/sigma^2;
                beta=-M*(H1d(:,:,k)'*Y(:,k))/sigma^2;
                
                affinity2(i,k)= -0.5*log(det(L+M)) -0.5*((alpha-beta)'*inv(L+M)*(alpha-beta));
                affinity3(i,k)= -0.5*log(det(L+M+20*eye(size(L,1)))) -0.5*((alpha-beta)'*inv(L+M)*(alpha-beta));
                affinity(i,k)= -0.5*((alpha-beta)'*inv(L+M)*(alpha-beta));%-0.5*log(det(L+M));
                wiener(i,k)=norm(alpha-beta);
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
        K=10;
        true_neighbor=zeros(K*d,K*d);
        for i=1:K*d
            for j=1:K*d
                true_neighbor(i,j)=(floor((i-1)/K)==floor((j-1)/K));
            end
        end
        error(1)=find_neighbors(affinity,true_neighbor,K);
        error(2)=find_neighbors(affinity2,true_neighbor,K);
        error(5)=find_neighbors(affinity3,true_neighbor,K);
        error(3)=find_neighbors(-wiener,true_neighbor,K);
        % error(4)=find_neighbors(crosscorr,true_neighbor,K);
        % error(5)=find_neighbors(crosscorr1,true_neighbor,K);
        error(4)=find_neighbors(crosscorr,true_neighbor,K);
        
        %
        %     %%%%%Metric
        %     Sort_aff=zeros(d*10,d*10);
        %     dim_aff=zeros(d*10,d*10);
        %     for i=1:d*10
        %         [Sort_aff(i,:),dim_aff(i,:)]=sort(affinity(i,:));
        %     end
        %     Sort_aff=flipdim(Sort_aff, 2);
        %     dim_aff=flipdim(dim_aff,2); % Sorted rowwise in descending order
        %
        %     %%%%%Wiener
        %     Sort_w=zeros(d*10,d*10);
        %     dim_w=zeros(d*10,d*10);
        %     for i=1:d*10
        %         [Sort_w(i,:),dim_w(i,:)]=sort(wiener(i,:));
        %     end
        %     %Sort_w=flipdim(Sort_w, 2);
        %     %dim_w=flipdim(dim_w,2); % Sorted rowwise in descending order
        %
        %     %%%%%Cross correlation
        %     Sort_cross=zeros(d*10,d*10);
        %     dim_cross=zeros(d*10,d*10);
        %     for i=1:d*10
        %         [Sort_cross(i,:),dim_cross(i,:)]=sort(crosscorr(i,:));
        %     end
        %     Sort_cross=flipdim(Sort_cross, 2);
        %     dim_cross=flipdim(dim_cross,2); % Sorted rowwise in descending order
        %
        %
        %     t3=tic;
        %     for p=1:10:d*10
        %         for i=p:p+9
        %             for k=1:10
        %                 if dim_aff(i,k) <= p+9 && dim_aff(i,k) >= p
        %                     true_aff(it)=true_aff(it)+1;
        %                 end
        %                 if dim_cross(i,k) <= p+9 && dim_cross(i,k) >= p
        %                     true_cross(it)=true_cross(it)+1;
        %                 end
        %                 if dim_w(i,k) <= p+9 && dim_w(i,k) >= p
        %                     true_w(it)=true_w(it)+1;
        %                 end
        %             end
        %         end
        %     end
        %
        %
        %
        %
        %     t3=toc(t3);
        %imagesc(log(affinity))
    end
    
    error
    allerror(ii,:)=error;
end

