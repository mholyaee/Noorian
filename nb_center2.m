function[Rating]=nb_center2(P,nb_RatingMatrix,k)
[~,C]=size(nb_RatingMatrix);
Rating=zeros(k,C);
%Trust=zeros(k,C);
nc=zeros(k,C);
for i=1:length(P)
    Rating(P(i)+1,:)=nb_RatingMatrix(i,:)+Rating(P(i)+1,:);
    %Trust(P(i)+1,:)=nb_TrustNetwork(i,:)+Trust(P(i)+1,:);
    ind=find(nb_RatingMatrix(i,:)>0);
    nc(P(i)+1,ind)=nc(P(i)+1,ind)+1;
end
for i=1:k
    Rating(i,:)=Rating(i,:)./nc(i,:);
end
Rating(isnan(Rating))=0;

