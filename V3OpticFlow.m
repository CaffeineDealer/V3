
names={'slu048a','slu047a','slu046b','slu045b','slu044a','slu022a','slu017b',...
    'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c'};
path='C:\research\data\SuperTuneSpkTrains\';
Fs=10000;
DprimeV3=[];
DprimeMT=[];
DprimeMST=[];
TunWMT=[];
TunWV3=[];
count=1;
for j=1:length(names)
params=load(['C:\research\data\RFiles\',names{1,j},'_TrialStructure.mat']);
if exist(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'],'file')==2
    load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'])
else
    load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end),'Chunum.mat'])
end
load(['C:\research\V3 things\V3 categorized\',names{1,j}(1:end-1),'_V3categ.mat']);
V3units=v3categ((v3categ(:,3)<=4),1:2);%|v3categ(:,3)==4
MTunits=v3categ(v3categ(:,3)==5,1:2);
for ci=1:size(V3units,1)+size(MTunits,1)
    if ci<=size(V3units,1)
        ch=V3units(ci,1);
        u=V3units(ci,2);
    else
        ch=MTunits(ci-size(V3units,1),1);
        u=MTunits(ci-size(V3units,1),2);
    end
firing=load(['C:\research\data\SuperTuneFiringMatrix\',names{1,j},num2str(ch),num2str(u),'firingMat']);
if size(firing.firing,2)<3
    firing1=squeeze(firing.firing(:,2,:));
else
   firing1=squeeze(firing.firing(:,2:3,:)); 
end
%spktrain=load([path,name,num2str(ch),num2str(u),'spktrain.mat']);
spktrainbl=load([path,names{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
baseline=squeeze(sum(spktrainbl.spktrain_bl,1))*Fs/size(spktrainbl.spktrain_bl,1);
bl=mean(baseline(:));
[maxfir,I]=max(firing1(:));
[ddir,typ,dpos]= ind2sub(size(firing1),I);
dirfir=squeeze(firing1(:,typ,dpos));
 keepCriteria= maxfir>2*bl;
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
    TunWidth=a(1);
      count=count+1;
%       figure
%       plot(xaxis,dirfir)
%       hold on
%       plot(xaxis,fitted)
 end
if ci<=size(V3units,1)
    DprimeV3=[DprimeV3,dprimeDots];
    TunWV3=[TunWV3,TunWidth];
else
    DprimeMT=[DprimeMT,dprimeDots];
    TunWMT=[TunWMT,TunWidth];
end
end
end
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