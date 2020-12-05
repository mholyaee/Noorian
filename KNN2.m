function RNNname=KNN2(neighbors,RatingMatrix,TrustNetwork,K,tempdis)  
nF=length(neighbors);
nb_RatingMatrix=RatingMatrix(neighbors,:);
nb_TrustNetwork=TrustNetwork(neighbors,:);
sim_all=zeros(nF,nF);
%T_i_j=zeros(nF,nF);% trust statement among the user i and user j
Degree2=sum(nb_TrustNetwork()'~=0);
AvgDegree2=mean(Degree2);
%[AllU,~]=size(TrustNetwork);
Dmax2=log(nF)/log(AvgDegree2);
W_i_j=zeros(nF,nF);
for i=1:nF
    nRatingi=sum(nb_RatingMatrix(i,:)~=0);
    sumRatingi=sum(nb_RatingMatrix(i,:));
    avgi=sumRatingi/nRatingi;
    Disi=tempdis(:,neighbors(i));    % tempdis obtained in line 44
    T_i_j=(Dmax2-Disi+1)/Dmax2;
    T_i_j(isinf(T_i_j))=0;% Is it true?!

    for j=i+1:nF
        temp_i=nb_RatingMatrix(i,:);
        nRatingj=sum(nb_RatingMatrix(j,:)~=0);
        sumRatingj=sum(nb_RatingMatrix(j,:));
        avgj=sumRatingj/nRatingj;
        %        
        temp_j=nb_RatingMatrix(j,:);
        cmmn_items=temp_i.*temp_j;
        cmmn_itm_ind=find(cmmn_items~=0);
        temp_i=temp_i(cmmn_itm_ind)-avgi;
        temp_j=temp_j(cmmn_itm_ind)-avgj;
        num=temp_i*temp_j';
        den=sqrt(temp_i*temp_i')*sqrt(temp_j*temp_j');
        sim_all(i,j)=num/den;
        if isnan(sim_all(i,j))==1% It should be checked. Two users which there isnt any common elements between them.
            sim_all(i,j)=0;
        end
        if sim_all(i,j)+T_i_j(j)~=0 && sim_all(i,j)*T_i_j(j)~=0
            W_i_j(i,j)=2*sim_all(i,j)*T_i_j(j)/(sim_all(i,j)+T_i_j(j));
        elseif sim_all(i,j)==0 && T_i_j(j)~=0
            W_i_j(i,j)=T_i_j(j);
        elseif sim_all(i,j)~=0 && T_i_j(j)==0
            W_i_j(i,j)=sim_all(i,j);
        else
            'ggg';
        end
        if W_i_j(i,j)<0 % based on the comment of author
            W_i_j(i,j)=0;
        end
    end
end
W_i_j=W_i_j+W_i_j';
%*****************************************
[~,I]=sort(W_i_j,2,'descend'); 
knn=I(:,1:K);

NN=zeros(nF,nF);
for i=1:nF
    for j=1:K
        NN(i,knn(i,j))=1;% he length of row i is M which is 0/1. if NN[i,j]=1 then j is Knn of i 
    end
end
RNN=NN';% Reverse NN

RNNname='neighbors.tab';
fid_RNN=fopen(RNNname,'w+');
for i=1:nF     
    for j=1:nF
        if RNN(i,j)==1
            fprintf(fid_RNN,'%d ',j);
        end
    end
    fprintf(fid_RNN,'\n');
end
fclose(fid_RNN);
