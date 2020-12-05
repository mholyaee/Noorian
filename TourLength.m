% Based on the Eqs 11 and 12.
function fit=TourLength(tour,model,TR,Ratings,tau,k)

%n=numel(tour);
[~,n]=size(Ratings);
R_hat=zeros(1,n);
%tour=[tour tour(1)];
ind=find((TR)>0);
nR=length(ind);
for i=1:nR
    %L=L+model.D(tour(i),tour(i+1));
    try
    %xx=Ratings(tour,ind(i))';
    %yy=tau(k,:);
        %L=sum(Ratings(tour,i)'.*tau(k,:));
        L=0;
        M=0;
        for t=1:length(tour)-1
        L=(Ratings(tour(t),ind(i))*tau(tour(t),tour(t+1)))+L;
        M=tau(tour(t),tour(t+1))+M;
        end
    catch
        fprintf('\n**error**\n');
    end
    R_hat(ind(i))=L/M;
end
%R_hat=round(R_hat);

fit=sum(abs(TR(ind)-R_hat(ind)))/length(ind);

end