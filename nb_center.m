function[Rating]=nb_center(P,nb_RatingMatrix,k)
[~,C]=size(nb_RatingMatrix);
Rating=zeros(k,C);
%Trust=zeros(k,C);
nc=zeros(1,k);
for i=1:length(P)
    Rating(P(i)+1,:)=nb_RatingMatrix(i,:)+Rating(P(i)+1,:);
    %Trust(P(i)+1,:)=nb_TrustNetwork(i,:)+Trust(P(i)+1,:);
    nc(P(i)+1)=nc(P(i)+1)+1;
end
for i=1:k
    Rating(i,:)=Rating(i,:)/nc(i);
end

