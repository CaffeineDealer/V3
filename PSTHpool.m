
plaid=0;
if plaid
    path='C:\research\data\PlaidSpkTrains\';
    names={'slu048d','slu047e','slu046e','slu045e','slu044d','slu022d','slu017c',...
    'ytu326a','ytu331b','ytu332d','ytu334b','ytu336a'};%grating tasks
else
    path='C:\research\data\SuperTuneSpkTrains\';
    names={'slu048a','slu047a','slu046b','slu045b','slu044a','slu022a','slu017b',...
    'ytu326c','ytu331a','ytu332b','ytu334c','ytu336c'};%dots tasks
end
bin=100;
%motiontype=1;
speed=1;
latV3=[];
latMT=[];
for motiontype=1:3
for j=1:length(names)
load(['C:\research\V3 things\V3 categorized\',names{1,j}(1:end-1),'_V3categ.mat']);
V3units=v3categ((v3categ(:,3)<=4),1:2);%|v3categ(:,3)==4
MTunits=v3categ(v3categ(:,3)==5,1:2);

for ci=1:size(V3units,1)+size(MTunits,1)%:length(chunum(:,1))%[1,3:8,11];
    if ci<=size(V3units,1)
        ch=V3units(ci,1);
        unit=V3units(ci,2);
    else
        ch=MTunits(ci-size(V3units,1),1);
        unit=MTunits(ci-size(V3units,1),2);
    end
    %if ch~=stimCh && ch~=stimCh+1
% baseline=load([path,name(1:end-1),condition,num2str(ch),num2str(unit),'spktrain_bl.mat']);
% spktrain_bl=baseline.spktrain_bl;
% psthbl=reshape(spktrain_bl,[4279 864]);
% Fs=baseline.Fs;
load([path,names{1,j}(1:end),num2str(ch),num2str(unit),'spktrain.mat']);
if motiontype<=size(spktrain,3)
if plaid
    spktrain=spktrain(:,:,speed,:,:,:,:);
    soso=reshape(spktrain,[size(spktrain,1) size(spktrain,2)*size(spktrain,3)*...
    size(spktrain,4)*size(spktrain,5)*size(spktrain,6)*size(spktrain,7)]); 
else
    spktrain=spktrain(:,:,motiontype,:,:,:,:);
    soso=reshape(spktrain,[size(spktrain,1) size(spktrain,2)*size(spktrain,3)*...
    size(spktrain,4)*size(spktrain,5)]);
end

%%
% fsfs=load([path2,name(1:end-1),'S',num2str(ch),num2str(unit),'firingMat.mat']);
% firingS=fsfs.firing;
% fnfn=load([path2,name(1:end-1),'N',num2str(ch),num2str(unit),'firingMat.mat']);
% firingN=fnfn.firing;
%%
count=1;
for k=1:bin:size(soso,1)-bin
    lala=sum(soso((k:k+bin),:),1);
    psthst(count)=mean(lala);
    count=count+1;
%      tata=sum(psthbl,1);
%     psthbl(mm,k)=mean(tata);   
end
peak=max(psthst);
index=find(psthst==peak(1));
if peak(1)>2*mean(psthst) && index(1)-5<10
[~,latency]=max(psthst);
    if ci<=size(V3units,1)
        latV3=[latV3,latency];
    else
        latMT=[latMT,latency];
    end
end

end
%(psthst/max(psthst))
%title([name,', ch unit',num2str(ch),num2str(unit)])
%hold on
end

end
latV3M{motiontype,1}=latV3;
latMTM{motiontype,1}=latMT;
latV3=[];
latMT=[];

end
%%
mean1=mean(latV3M{1,1})*10;
mean2=mean(latV3M{2,1})*10;
mean3=mean(latV3M{3,1})*10;
figure
histogram(latV3M{1,1},30)
hold on
histogram(latV3M{2,1},30)
hold on
histogram(latV3M{3,1},20)
title('V3 peak response latency')
legend(['Translation, mean=',num2str(mean1),'ms'],['Spirals, mean=',num2str(mean2),'ms'],['Deformation, mean=',num2str(mean3),'ms'])

mean11=mean(latMTM{1,1})*10;
mean22=mean(latMTM{2,1})*10;
mean33=mean(latMTM{3,1})*10;
figure
histogram(latMTM{1,1},30)
hold on
histogram(latMTM{2,1},30)
hold on
histogram(latMTM{3,1},20)
title('MT peak response latency')
legend(['Translation, mean=',num2str(mean11),'ms'],['Spirals, mean=',num2str(mean22),'ms'],['Deformation, mean=',num2str(mean33),'ms'])
