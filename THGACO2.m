function [MAE]=THGACO2(varargin)
% Parameters
teta=0.63;% Based on the paper (Page 13, First paragraph)
k_EndUser=10;
F_Dmax=1;
if nargin==0
    [datafile,path]=uigetfile('D:\My Research\My projects\MSNoorian\1_SourceCode\*.txt','Select Dataset');
    [trustfile,~]=uigetfile('D:\My Research\My projects\MSNoorian\1_SourceCode\*.txt','Select Dataset');
    datafile=[path,datafile];
    trustfile=[path,trustfile];
    
    fid_data=fopen(datafile);
    file_data=textscan(fid_data,'%s','delimiter','\n');
    content=char(file_data{1,1});
    [R,~]=size(content);
    for r=1:R
        line=textscan(content(r,:), '%s','delimiter',' ');
        line=line{1};
        Rating(r,:)=[str2num(line{1}),str2num(line{2}),str2num(line{3})];
    end
    fclose(fid_data);
    Sp_Rating=spconvert(Rating);   % Convert 2 Sparse format
    RatingMatrix=full(Sp_Rating);
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fid_trust=fopen(trustfile);
    file_trust=textscan(fid_trust,'%s','delimiter','\n');
    content=char(file_trust{1,1});
    [R,~]=size(content);
    
    for r=1:R
        line=textscan(content(r,:), '%s','delimiter',' ');
        line=line{1};
        Trust(r,:)=[str2num(line{1}),str2num(line{2}),str2num(line{3})];
    end
    fclose(fid_trust);
    Sp_Trust=spconvert(Trust);   % Convert 2 Sparse format
    TrustNetwork=full(Sp_Trust);
    [nU,~]=size(RatingMatrix);
    TargetUser=round((nU-1).*(rand(1))) + 1;
    TargetUser=1444;
else
    TargetUser=varargin{1};
    RatingMatrix=varargin{2};
    TrustNetwork=varargin{3};
    Sp_Trust=varargin{4};
end
% Step1: Ranking
% Compute trust based similarity value
Degree=sum(TrustNetwork()'~=0);
AvgDegree=mean(Degree);
[AllU,~]=size(TrustNetwork);
tempdis=graphallshortestpaths(Sp_Trust);
tempdis2=tempdis;
tempdis2(isinf(tempdis2))=0;
if F_Dmax==0
    Dmax=log(AllU)/log(AvgDegree);% equal to the average path length of the graph as
else
    Dmax=max(max(tempdis2))*2;
end
DisTargetUser=tempdis(:,TargetUser);% Trust distance to target user
T_a_u=(Dmax-DisTargetUser+1)/Dmax;
T_a_u(isinf(T_a_u))=0;% Is it true. based on the commands in lines 75-82
% Computing PCC
% Computing r_hat_a
nRatingTarget=sum(RatingMatrix(TargetUser,:)~=0);
sumRatingTarget=sum(RatingMatrix(TargetUser,:));
avgTarget=sumRatingTarget/nRatingTarget;
[nU,~]=size(RatingMatrix);
sim_i=zeros(1,nU);
W_a_u=zeros(1,nU);
for i=1:nU
    if i~=TargetUser
        % Computing r_hat_u
        nRating=sum(RatingMatrix(i,:)~=0);
        sumRating=sum(RatingMatrix(i,:));
        avgi=sumRating/nRating;
        % Finding common rated items between u and a
        temp_a=RatingMatrix(TargetUser,:);
        temp_u=RatingMatrix(i,:);
        cmmn_items=temp_a.*temp_u;
        cmmn_itm_ind=find(cmmn_items~=0);
        temp_a=temp_a(cmmn_itm_ind)-avgTarget;
        temp_u=temp_u(cmmn_itm_ind)-avgi;
        % Computing Eq.6
        num=temp_a*temp_u';
        den=sqrt(temp_a*temp_a')*sqrt(temp_u*temp_u');
        sim_i(i)=num/den;
        if isnan(sim_i(i))==1
            sim_i(i)=0;
        end
        % Computing trust-based similarity values
        if sim_i(i)+T_a_u(i)~=0 && sim_i(i)*T_a_u(i)~=0
            W_a_u(i)=2*sim_i(i)*T_a_u(i)/(sim_i(i)+T_a_u(i));
        elseif sim_i(i)==0 && T_a_u(i)~=0
            W_a_u(i)=T_a_u(i);
        elseif sim_i(i)~=0 && T_a_u(i)==0
            W_a_u(i)=sim_i(i);
        end
    end
end
% neighbors of the target user
neighbors=find(W_a_u>=teta);
while length(neighbors)<k_EndUser
    teta=teta-0.1;
    neighbors=find(W_a_u>=teta);
end
% Creating neighbour distance matrix
nF=length(neighbors);
nb_RatingMatrix=RatingMatrix(neighbors,:);
nb_TrustNetwork=TrustNetwork(neighbors,:);
%% Hypergraph section
K=5;  %  parameter k in the SkNN, one can change it.
C=2;   %  parameter cs the support count threshold, one can change it.
RNNname=KNN2(neighbors,RatingMatrix,TrustNetwork,K,tempdis);
system(['fpgrowth.exe -tm -s-',num2str(C),' ',RNNname,' ',RNNname,'.out']); %Call an external executable file fpgrowth.exe
fid_RNNout=fopen([RNNname,'.out'],'r');
fid_graph=fopen(['neighbors','.Graph'],'w+');
fprintf(fid_graph,'                             \n');
NUMedge=0;
while ~feof(fid_RNNout)
    tline=fgetl(fid_RNNout);
    if tline==-1
        return;
    end
    index=strfind(tline,'(');
    newtline=[tline(index+1:end-1),' ',tline(1:index-2),'\n'];
    fprintf(fid_graph,newtline);
    NUMedge=NUMedge+1;
end
frewind(fid_graph);% set pointer to the start of the file
fprintf(fid_graph,'%d %d 1',NUMedge,nF);%%% check
fclose(fid_RNNout);
fclose(fid_graph);
Nparts=k_EndUser;% #############Number of partitions###############
UBfactor=1;
Nruns=10;  %default:10 range:1...10
CType=5;   %default:1  range:1...5
RType=3;   %default:1  range:1...3
Vcycle=1;  %default:1  range:1...3
Reconst=0; %default:0  range:0...1
dbglvl=0; %range:0 or 24
system(['hmetis ','neighbors','.Graph ',num2str(Nparts),' ',num2str(UBfactor),' ',num2str(Nruns),' ',num2str(CType),' ',num2str(RType),' ',num2str(Vcycle),' ',num2str(Reconst),' ',num2str(dbglvl)]);
%Call an external executable file hmetis. e.g. hmetis 1.Graph 2 1 10 1 1 1 0 24
P=load(['neighbors','.Graph.part.',num2str(Nparts)]);% results of paritioning
[nb_RatingMatrix]=nb_center2(P,nb_RatingMatrix,k_EndUser);
% The constructed hypergraph is partitioned into "k_EndUser" clusters.
% Now we should find the center of each cluster and get to aco
% Create User graph
nF=k_EndUser;
%% Step II – Weighing
sim_all=zeros(nF,nF);
%T_i_j=zeros(nF,nF);% trust statement among the user i and user j
Degree2=sum(nb_TrustNetwork()'~=0);
AvgDegree2=mean(Degree2);
%[AllU,~]=size(TrustNetwork);
%Dmax2=log(nF)/log(AvgDegree2);
%W_i_j=zeros(nF,nF);
for i=1:nF
    nRatingi=sum(nb_RatingMatrix(i,:)~=0);
    sumRatingi=sum(nb_RatingMatrix(i,:));
    avgi=sumRatingi/nRatingi;
    %     Disi=tempdis(:,neighbors(i));    % tempdis obtained in line 44
    %     T_i_j=(Dmax2-Disi+1)/Dmax2;
    %     T_i_j(isinf(T_i_j))=0;% Is it true?!
    
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
        %         if sim_all(i,j)+T_i_j(j)~=0 && sim_all(i,j)*T_i_j(j)~=0
        %             W_i_j(i,j)=2*sim_all(i,j)*T_i_j(j)/(sim_all(i,j)+T_i_j(j));
        %         elseif sim_all(i,j)==0 && T_i_j(j)~=0
        %             W_i_j(i,j)=T_i_j(j);
        %         elseif sim_all(i,j)~=0 && T_i_j(j)==0
        %             W_i_j(i,j)=sim_all(i,j);
        %         else
        %             'ggg';
        %         end
        %         if W_i_j(i,j)<0 % based on the comment of author
        %             W_i_j(i,j)=0;
        %         end
    end
end
sim_all=sim_all+sim_all';
%% Hypergraph section
%  K=5;  %  parameter k in the SkNN, one can change it.
%  C=2;   %  parameter cs the support count threshold, one can change it.
%  RNNname=KNN2(neighbors,RatingMatrix,TrustNetwork,K,tempdis);
%  system(['fpgrowth.exe -tm -s-',num2str(C),' ',RNNname,' ',RNNname,'.out']); %Call an external executable file fpgrowth.exe
%  fid_RNNout=fopen([RNNname,'.out'],'r');
%  fid_graph=fopen(['neighbors','.Graph'],'w+');
%  fprintf(fid_graph,'                             \n');
%  NUMedge=0;
%  while ~feof(fid_RNNout)
%      tline=fgetl(fid_RNNout);
%      if tline==-1
%          return;
%      end
%      index=strfind(tline,'(');
%      newtline=[tline(index+1:end-1),' ',tline(1:index-2),'\n'];
%      fprintf(fid_graph,newtline);
%      NUMedge=NUMedge+1;
%  end
%  frewind(fid_graph);% set pointer to the start of the file
%  fprintf(fid_graph,'%d %d 1',NUMedge,nF);%%% check
%  fclose(fid_RNNout);
%  fclose(fid_graph);
%  Nparts=k_EndUser;% #############Number of partitions###############
%  UBfactor=1;
%  Nruns=10;  %default:10 range:1...10
%  CType=5;   %default:1  range:1...5
%  RType=3;   %default:1  range:1...3
%  Vcycle=1;  %default:1  range:1...3
%  Reconst=0; %default:0  range:0...1
%  dbglvl=0; %range:0 or 24
%  system(['hmetis ','neighbors','.Graph ',num2str(Nparts),' ',num2str(UBfactor),' ',num2str(Nruns),' ',num2str(CType),' ',num2str(RType),' ',num2str(Vcycle),' ',num2str(Reconst),' ',num2str(dbglvl)]);
%  %Call an external executable file hmetis. e.g. hmetis 1.Graph 2 1 10 1 1 1 0 24
%  P=load(['neighbors','.Graph.part.',num2str(Nparts)]);% results of paritioning
%  nb_center();
% The constructed hypergraph is partitioned into "k_EndUser" clusters.
% Now we should find the center of each cluster and get to aco
% Create User graph
%% Applying ACO
W=sim_all;
W(W<=0)=0.0001;  %Will be checked
%W=1./W;
%[BestSol,BestAnt,tau]=aco(W,RatingMatrix(TargetUser,:),nb_RatingMatrix);
%% Prediction
%BestWeight=tau(BestAnt,:);
% BestL=length(BestSol.Tour);
% BestWeight=zeros(1,BestL-1);
% for b=1:BestL-1
%     BestWeight(b)=tau(BestSol.Tour(b),BestSol.Tour(b+1));
% end
% [s,ind]=sort(BestWeight,'descend');
% for j=1:length(s)-1
%     sorted(j)=BestSol.Tour(ind(j));
% end
% temp=BestSol.Tour;
% Selecting k users
% K_users=sorted(1:k_EndUser-1);
% K_weights=s(1:k_EndUser-1);
[K_users,nR]=size(nb_RatingMatrix);
R_hat=zeros(1,nR);
%M=K_users;
for i=1:nR
    %xx=nb_RatingMatrix(K_users,i)';
    %yy=K_weights;
    L=sum(nb_RatingMatrix(:,i));
    if L>0
        M=sum(nb_RatingMatrix(:,i)>0);
        R_hat(i)=L/M;
    end
end
%R_hat=round(R_hat);
%% Evaluation
% Mean Average Error(MAE)
TR=RatingMatrix(TargetUser,:);
ind=find((TR)>0);
MAE=sum(abs(TR(ind)-R_hat(ind)))/length(ind);

end