function []=THGACO()
% Parameters
teta=0;
[datafile,path]=uigetfile('D:\My Research\My projects\MSNoorian\0_Dataset\*.txt','Select Dataset');
[trustfile,~]=uigetfile('D:\My Research\My projects\MSNoorian\0_Dataset\*.txt','Select Dataset');
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
TargetUser=(nU-1).*round(rand(1)) + 1;
% Step1: Ranking
% Compute trust based similarity value
Degree=sum(TrustNetwork()'~=0);
AvgDegree=mean(Degree);
[AllU,~]=size(TrustNetwork);
Dmax=log(AllU)/log(AvgDegree);% equal to the average path length of the graph as
tempdis=graphallshortestpaths(Sp_Trust);
DisTargetUser=tempdis(:,TargetUser);
T_a_u=(Dmax-DisTargetUser+1)/Dmax;
% Computing PCC
% Computing r_hat_a
nRatingTarget=sum(RatingMatrix(TargetUser,:)~=0);
sumRatingTarget=sum(RatingMatrix(TargetUser,:));
avgTarget=sumRatingTarget/nRatingTarget;
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
       num=temp_a*temp_u';
       den=sqrt(temp_a*temp_a')*sqrt(temp_u*temp_u');
       sim_i(i)=num/den;
       % Computing trust-based similarity values
       if sim_i(i)+T_a_u(i)~=0 || sim_i(i)*T_a_u(i)~=0
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
%% Step II – Weighing
% Hypergraph section

% Create User graph


end