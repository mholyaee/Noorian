%Cold start
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
n=1;
for i=1:nU
    if sum(RatingMatrix(i,:)>0)<5
        ColdUsers(n)=i;
        n=n+1;
    end
end
n=length(ColdUsers);
MAE=zeros(2,n);
for i=1:n
    MAE(1,i)=THGACO2(ColdUsers(i),RatingMatrix,TrustNetwork,Sp_Trust);
    MAE(2,i)=TACO2(ColdUsers(i),RatingMatrix,TrustNetwork,Sp_Trust);
end