function x=find_neighbors(X,true_neighbor,d)
[d1,d2]=size(X);
X=X+10*max(max(abs(X)))*eye(size(X,1));
X1=zeros(size(X));
for i=1:d1
temp=sort(X(i,:),'descend');
cutoff=temp(d);
X1(i,:)=(X(i,:)>=cutoff);
end
x=sum(sum(X1.*true_neighbor))/sum(sum(true_neighbor));