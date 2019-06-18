
names={'slu062a','slu058c','slu055a','slu050a','slu048a','slu047a','slu046b','slu045b','slu044a','slu023a','slu022a','slu017b'};
%{'slu062a','slu058c','slu055a','slu050a','slu048a','slu047a','slu046b','slu045b','slu044a','slu023a','slu022a','slu017b',...
%'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c'};
path='C:\research\data\SuperTuneSpkTrains\';
Fs=10000;
DprimeV3=[];
DprimeMT=[];
DprimeMST=[];
TunWMT=[];
TunWV3=[];
blV3=[];
blMT=[];
tuningWidthAll=[];
count=1;
count2=1;
motiontype1=2;%1 or 2
motiontype2=2:3;%1 or 2:3
for j=1:length(names)
params=load(['C:\research\data\RFiles\',names{1,j},'_TrialStructure.mat']);
if exist(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'],'file')==2
    load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'])
else
    load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end),'Chunum.mat'])
end
load(['C:\research\V3 things\V3 categorized\',names{1,j}(1:end-1),'_V3categ2.mat']);
v3categ=sortrows(v3categ2);
V3units=v3categ((v3categ(:,3)<=4),1:2);%|v3categ(:,3)==4
MTunits=v3categ(v3categ(:,3)==5,1:2);
tuningW=zeros(length(v3categ(:,1)),2);
for ci=1:size(V3units,1)+size(MTunits,1)
    if ci<=size(V3units,1)
        ch=V3units(ci,1);
        u=V3units(ci,2);
    else
        ch=MTunits(ci-size(V3units,1),1);
        u=MTunits(ci-size(V3units,1),2);
    end
firing=load(['C:\research\data\SuperTuneFiringMatrix\',names{1,j},num2str(ch),num2str(u),'firingMat']);
if size(firing.firing,2)<3 %******* make it more general*************************************
    firing1=firing.firing(:,motiontype1,:,:,:);
else
   firing1=firing.firing(:,motiontype2,:,:,:); 
end
spktrain=load([path,names{1,j},num2str(ch),num2str(u),'spktrain.mat']);
% time,directions,numMotion,rows*columns,trialsPerFeature,sizes,coherences
spktrainbl=load([path,names{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
baseline=squeeze(sum(spktrainbl.spktrain_bl,1))*Fs/size(spktrainbl.spktrain_bl,1);
allstimfir=squeeze(sum(spktrain.spktrain,1))*Fs/size(spktrain.spktrain,1);

%
% TrialsSpks=spktrain.spktrain(:,ddir,typ,dpos,:,dsiz,dcoh);
% Trialsfiring= squeeze(sum(TrialsSpks,1))*Fs/size(TrialsSpks,1);
% TrialsSpksbl=spktrainbl.spktrain_bl(:,ddir,typ,dpos,:,dsiz,dcoh);
% Trialsfiringbl= squeeze(sum(TrialsSpksbl,1))*Fs/size(TrialsSpksbl,1);
 [h,p] = ttest(baseline(:),allstimfir(:));
 keepCriteria=(p<=0.05);
%
bl=mean(baseline(:));
[maxfir,I]=max(firing1(:));
[ddir,typ,dpos,dsiz,dcoh]= ind2sub(size(firing1),I);
dirfir=squeeze(firing1(:,typ,dpos,dsiz,dcoh));
 %keepCriteria= maxfir>3*bl;
 
 count2=count2+1;
 if keepCriteria
    goodch(count)=ch;
    goodunit(count)=u;
    muPrefdirD=dirfir(ddir);
    %VarPrefdirD=var(dirdotspktrain(ddir,:))?
    if ddir<=params.file.taskDialogValues.numberOfDirections/2
        muNulldirD=dirfir(ddir+(params.file.taskDialogValues.numberOfDirections/2));
       % VarNulldirD=var(dirdotspktrain(ddir+(dotsparams.file.taskDialogValues.numberOfDirections/2),:));
    else
        muNulldirD=dirfir(ddir-(params.file.taskDialogValues.numberOfDirections/2));
        %VarNulldirD=var(dirdotspktrain(ddir-(dotsparams.file.taskDialogValues.numberOfDirections/2),:));
    end
%     dprimeDots(count)=(muPrefdirD- muNulldirD)/sqrt((VarPrefdirD+VarNulldirD)/2);
      dprimeDots=(muPrefdirD- muNulldirD)/((muPrefdirD-bl)+ (muNulldirD-bl));
      
    a=mtfit(dirfir');
    xaxis=0:360/length(dirfir):(length(dirfir)-1)*(360/length(dirfir));
    fitted=vonMises(a,xaxis*pi/180);
    %TunWidth=1/a(1);
    peak=find(fitted==max(fitted));
    k=floor(length(dirfir)/2)-peak;
    if k>0
        fittedshift=circshift(fitted,k);
    else
        fittedshift=circshift(fitted,8+k);
    end
    [pks,locs,TunWidth,~] = findpeaks(fittedshift,xaxis);
    tuningW(ci,:)=[TunWidth,ch];
    if ci<=size(V3units,1)
        DprimeV3=[DprimeV3,dprimeDots];
        TunWV3=[TunWV3,TunWidth];
        blV3=[blV3,bl];
    else
        DprimeMT=[DprimeMT,dprimeDots];
        TunWMT=[TunWMT,TunWidth];
          blMT=[blMT,bl];
    end
    count=count+1;
%       figure
%       plot(xaxis,dirfir)
%       hold on
%       plot(xaxis,fitted)

 end


end
% figure
%  scatter(tuningW(:,2),tuningW(:,1))
%  
 tuningWidthAll=[tuningWidthAll;tuningW]; 
%plot clustering for this recording
%contour
end
%%
figure
 scatter(tuningWidthAll(:,2)*0.15,tuningWidthAll(:,1))
 ylabel('Tuning Width')
 xlabel('depth [mm]')
 
figure
scatter(TunWV3,blV3)
hold on
scatter(TunWMT,blMT)
legend('V3','MT')
ylabel('Tuning Width')
xlabel('mean baseline')
 
figure
histogram(blV3,50)
hold on
histogram(blMT,50)
legend('V3','MT')
title('mean baseline')
%%
meanMT=mean(TunWMT);
meanV3=mean(TunWV3);
figure
histogram(TunWMT,60)
hold on
histogram(TunWV3,40)
legend(['MT tuning width, mean =',num2str(meanMT)],['V3 tuning width, mean =',num2str(meanV3)])
title('Direction Tuning width for complex motion (VonMises)')
% DprimeV3=[DprimeV3,dprimeDots(1:16)];
% DprimeMT=[DprimeMT,dprimeDots([17:23])];%,dprimeDots([34:35])];
% DprimeMST=[DprimeMST,dprimeDots(24:29)];%, dprimeDots(36:48)];
% figure
% scatter(goodch,dprimeDots)
% figure
% hist(dprimeDots(1:20),20)
% hold on
% histogram(dprimeDots(21:24),8)
% hold on
% histogram(dprimeDots(25:end),8)
% title(['optic flow direction selectivity']);
% legend('V3','MT','MST')
% figure
% histogram(DprimeMT,15)
% hold on
% histogram(DprimeV3,15)
% legend('MT','V3')
% hold on
% histogram(DprimeMST,10)
% title(['optic flow direction selectivity']);
% legend('V3','MT','MST')
datasettt=[TunWMT;TunWV3(1:149)]';
[p,tbl,stats] = anova2(datasettt,149)
[h,p] = ttest(TunWMT,TunWV3(1:149))
%%
datasettt=[TunWMT1;TunWMT]';
[p,tbl,stats] = anova2(datasettt,149)
[h,p] = ttest(TunWMT1,TunWMT)
%%
datasettt=[TunWV31;TunWV3]';
[p,tbl,stats] = anova2(datasettt,190)
[h,p] = ttest(TunWV31,TunWV3)
%%
datasettt=[TunWMT1;TunWMT]';
[p,tbl,stats] = anova2(datasettt,149)
[h,p] = ttest(TunWMT1,TunWMT)
%%
datasettt=[DprimeV31;DprimeV3]';
[p,tbl,stats] = anova2(datasettt,190)
[h,p] = ttest(DprimeV31,DprimeV3)
%%
datasettt=[DprimeMT1;DprimeMT]';
[p,tbl,stats] = anova2(datasettt,149)
[h,p] = ttest(DprimeMT1,DprimeMT)
%%
datasettt=[DprimeMT;DprimeV3(1:149)]';
[p,tbl,stats] = anova2(datasettt,149)
[h,p] = ttest(DprimeMT,DprimeV3(1:149))
%% other tests
d = computeCohen_d(TunWMT, TunWV3)

[p,h,stats] = ranksum(TunWMT,TunWV3)
figure
histogram(TunWMT,40)
hold on
histogram(TunWV3,40)
legend(['MT tuning width'],['V3 tuning width'])
title(['Direction Tuning width for complex motion (VonMises), p value ',num2str(p)])
%
[p,h,stats] = ranksum(TunWMT1,TunWMT)
figure
histogram(TunWMT1,40)
hold on
histogram(TunWMT,40)
legend(['translational motion'],['complex motion'])
title(['Direction Tuning width for MT (VonMises), p value ',num2str(p)])
%
[p,h,stats] = ranksum(TunWV31,TunWV3)
figure
histogram(TunWV31,40)
hold on
histogram(TunWV3,40)
legend(['translational motion'],['complex motion'])
title(['Direction Tuning width for V3 (VonMises), p value ',num2str(p)])
%
[p,h,stats] = ranksum(TunWMT1,TunWV31)
figure
histogram(TunWMT1,40)
hold on
histogram(TunWV31,40)
legend(['MT tuning width'],['V3 tuning width'])
title(['Direction Tuning width for simple motion (VonMises), p value ',num2str(p)])
%
%[p,h,stats] = ranksum(TunWMT1-TunWMT,TunWV31-TunWV3)
% [h,p] = ttest(TunWMT1-TunWMT,TunWV31(1:149)-TunWV3(1:149))
% figure
% histogram(TunWMT1-TunWMT,40)
% hold on
% histogram(TunWV31-TunWV3,40)
% legend(['MT tuning width'],['V3 tuning width'])
% title(['Direction Tuning width, difference between simple and complex motion (VonMises), p value ',num2str(p)])
% median(TunWMT1-TunWMT)
% median(TunWV31-TunWV3)
%
[p,h,stats] = ranksum([TunWMT1,TunWMT],[TunWV31,TunWV3])
figure
histogram([TunWMT1,TunWMT],40)
hold on
histogram([TunWV31,TunWV3],40)
legend(['MT tuning width'],['V3 tuning width'])
title(['Direction Tuning width, all motion types (VonMises), p value ',num2str(p)])
median([TunWMT1,TunWMT])
median([TunWV31,TunWV3])