% get firing rate matrices
% compute d prime
% compute d prime difference per unit? then scatter plot for both V3 and
% MT?
dotsnames={'slu048b','slu047b','slu046c','slu045c','slu044b','slu022b','slu017b',...
    'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c'};
gratnames={'slu048d','slu047e','slu046e','slu045e','slu044d','slu022d','slu017c',...
    'ytu326a','ytu331b','ytu332d','ytu334b','ytu336a'};
% gratnames={'slu017c'};
% dotsnames={'slu017b'};
count=1;
dots=1;
grat=1;
comparison=0;
contrast=100;%40;%100
coherence=100;%55;%100
Fs=10000;
timeWin=(0.05*Fs:0.4*Fs);
DIThresh=0:0.05:1.5;

pathD='C:\research\data\SuperTuneSpkTrains\';
pathG='C:\research\data\PlaidSpkTrains\';
V3GratdP=[]; MTGratdP=[]; V3DotsdP=[]; MTDotsdP=[];

for j=1:length(gratnames)
load(['C:\research\V3 things\V3 categorized\',dotsnames{1,j}(1:end-1),'_V3categ.mat']);
V3units=v3categ((v3categ(:,3)<=4),1:2);%|v3categ(:,3)==4
MTunits=v3categ(v3categ(:,3)==5,1:2);
dotsparams=load(['C:\research\data\RFiles\',dotsnames{1,j},'_TrialStructure.mat']);
gratparams=load(['C:\research\data\RFiles\',gratnames{1,j},'_TrialStructure.mat']);
%%
mincoh=dotsparams.file.taskDialogValues.maxCoherence-(dotsparams.file.taskDialogValues.numberOfCoherences-1)*...
    (dotsparams.file.taskDialogValues.coherenceStep);
coherences=mincoh:dotsparams.file.taskDialogValues.coherenceStep:dotsparams.file.taskDialogValues.maxCoherence;
coherenceidx=find(coherences==coherence);

contrasts=gratparams.file.taskDialogValues.contrastArray;
contrastidx=find(contrasts>=contrast-15 & contrasts<=contrast+15);
%%
for ci=1:size(V3units,1)+size(MTunits,1)
    if ci<=size(V3units,1)
        ch=V3units(ci,1);
        u=V3units(ci,2);
    else
        ch=MTunits(ci-size(V3units,1),1);
        u=MTunits(ci-size(V3units,1),2);
    end
if dots
% dotfiring=load(['C:\research\data\SuperTuneFiringMatrix\',dotsnames{1,j},num2str(ch),num2str(u),'firingMat']);
% dotfiring1=squeeze(dotfiring.firing(:,1,:));
dotspktrain=load([pathD,dotsnames{1,j},num2str(ch),num2str(u),'spktrain.mat']);
dotfiring=sum(dotspktrain.spktrain(timeWin,:,:,:,:,:,:),1)*Fs/length(timeWin);%size(gratspktrain.spktrain,1);
    dotfiring=mean(dotfiring,5);
sizedotstim=size(dotfiring);%direction,motiontype,pos
speed=dotsparams.file.taskDialogValues.speedStepDegPerSec^dotsparams.file.taskDialogValues.minSpeedDegPerSec;

numsizdot=dotsparams.file.taskDialogValues.numberOfRadii;
numgridposdot=dotsparams.file.taskDialogValues.superTuneColumns*dotsparams.file.taskDialogValues.superTuneRows;

dotspktrainbl=load([pathD,dotsnames{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
baseline=squeeze(sum(dotspktrainbl.spktrain_bl,1))*Fs/size(dotspktrainbl.spktrain_bl,1);
blDots=mean(baseline(:));
 [maxfirD,ID]=max(dotfiring(:));
 [dtim,ddir,motTyp,dpos,dtrial,dsiz,dcoh]= ind2sub(size(dotfiring),ID);
 dirfirdot=squeeze(dotfiring(:,:,motTyp,dpos,dtrial,dsiz,coherenceidx));
%  if numgridposdot>1
%      [ddir,dpos]= ind2sub(size(dotfiring),ID);
%      dirfirdot=squeeze(dotfiring(:,dpos));
%      dirdotspktrain=squeeze(sum(dotspktrain.spktrain(:,:,1,dpos,:),1))*Fs/size(dotspktrain.spktrain,1);
%  elseif numgridposdot==1
%      [ddir]= ID;
%      dirfirdot=squeeze(dotfiring);
%      dirdotspktrain=squeeze(sum(dotspktrain.spktrain(:,:,1,:),1))*Fs/size(dotspktrain.spktrain,1);
%  end
 keepCriteriaD= maxfirD>2*blDots;%mean(dotfiring(:));
end


if grat
%     gratingfiring=load(['C:\research\data\PlaidFiringMatrix\',gratnames{1,j},num2str(ch),num2str(u),'firingMat.mat']);
%     gratingfiring=gratingfiring.firing;
    gratspktrain=load([pathG,gratnames{1,j},num2str(ch),num2str(u),'spktrain.mat']);
    gratingfiring=sum(gratspktrain.spktrain(timeWin,:,:,:,:,:,:),1)*Fs/length(timeWin);%size(gratspktrain.spktrain,1);
    gratingfiring=mean(gratingfiring,7);
    sizegratstim=size(gratingfiring); %time,directions,speed,pos,size,contrast, tria#
    gratspktrainbl=load([pathG,gratnames{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
    baseline2=squeeze(sum(gratspktrainbl.spktrain_bl,1))*Fs/size(gratspktrainbl.spktrain_bl,1);
    blGrat=mean(baseline2(:));
    speedsgrat=zeros(sizegratstim(3),1);
    for i=1:sizegratstim(3)
        speedsgrat(i)=(gratparams.file.taskDialogValues.TFbase^...
            (gratparams.file.taskDialogValues.cyclesPerSecond+i-1))/...
            (2*gratparams.file.taskDialogValues.spatialFrequency1);
        %temporal freq/spatial freq
        %spatial frequesncy multiplied by 2 because it's multiplied by 2 in
        %monkeylab code (some crazy person programmed it that way)
    end
    
    [maxfirG,IG]=max(gratingfiring(:));
     [gtime,gdir,gspd,gpos,gsiz,gcont]= ind2sub(size(gratingfiring),IG); 
    
    if comparison==1
        spdIdx=find(speedsgrat==speed);
        if isempty(spdIdx)
            display('stimulus speeds do not match')
            sorry
        end
    else
        spdIdx=gspd;%speedsgrat(bestSpdGrat);
    end

     dirfirgrat=squeeze(gratingfiring(:,:,gspd,gpos,gsiz,contrastidx));
     keepCriteriaG= maxfirG>2*blGrat;%mean(gratingfiring(:));
end    

if comparison
    keepCriteria = ddir==gdir;
elseif dots && ~grat
    keepCriteria=keepCriteriaD;    
elseif grat && ~dots
    keepCriteria=keepCriteriaG;    
elseif grat && dots
    keepCriteria=keepCriteriaG||keepCriteriaD;
end

 if keepCriteria
    goodch(count)=ch;
    goodunit(count)=u;
    if dots && ~isempty(dirfirdot)
    muPrefdirD=dirfirdot(ddir);
    %VarPrefdirD=var(dirdotspktrain(ddir,:))?
    if ddir<=dotsparams.file.taskDialogValues.numberOfDirections/2
        muNulldirD=dirfirdot(ddir+(dotsparams.file.taskDialogValues.numberOfDirections/2));
       % VarNulldirD=var(dirdotspktrain(ddir+(dotsparams.file.taskDialogValues.numberOfDirections/2),:));
    else
        muNulldirD=dirfirdot(ddir-(dotsparams.file.taskDialogValues.numberOfDirections/2));
        %VarNulldirD=var(dirdotspktrain(ddir-(dotsparams.file.taskDialogValues.numberOfDirections/2),:));
    end
%     dprimeDots(count)=(muPrefdirD- muNulldirD)/sqrt((VarPrefdirD+VarNulldirD)/2);
    if ci<=size(V3units,1)
        dprimeDotsV3=(muPrefdirD- muNulldirD)/((muPrefdirD-blDots)+ (muNulldirD-blDots));
    else
        dprimeDotsMT=(muPrefdirD- muNulldirD)/((muPrefdirD-blDots)+ (muNulldirD-blDots));
    end
    else
        if ci<=size(V3units,1)
        dprimeDotsV3=[];
        else
        dprimeDotsMT=[];
         end
    end
    
    if grat && ~isempty(dirfirgrat)
    muPrefdirG=dirfirgrat(gdir);
    %VarPrefdirG=var(dirgratspktrain(gdir,:))?
    if gdir<=gratparams.file.taskDialogValues.numberOfDirections/2
        muNulldirG=dirfirgrat(gdir+(gratparams.file.taskDialogValues.numberOfDirections/2));
       %VarNulldirG=var(dirgratspktrain(gdir+(gratparams.file.taskDialogValues.numberOfDirections/2),:));
    else
        muNulldirG=dirfirgrat(gdir-(gratparams.file.taskDialogValues.numberOfDirections/2));
       %VarNulldirG=var(dirgratspktrain(gdir-(gratparams.file.taskDialogValues.numberOfDirections/2),:));
    end
%     dprimeGrat(count)=(muPrefdirG- muNulldirG)/sqrt((VarPrefdirG+VarNulldirG)/2);
    if ci<=size(V3units,1)
        dprimeGratV3=(muPrefdirG- muNulldirG)/((muPrefdirG-blGrat)+ (muNulldirG-blGrat));
    else
        dprimeGratMT=(muPrefdirG- muNulldirG)/((muPrefdirG-blGrat)+ (muNulldirG-blGrat));
    end
    else
        if ci<=size(V3units,1)
        dprimeGratV3=[];
        else
        dprimeGratMT=[];
         end
    end

    count=count+1;
    if ci<=size(V3units,1)
        if grat
            V3GratdP=[V3GratdP,dprimeGratV3];
        end
        if dots
            V3DotsdP=[V3DotsdP,dprimeDotsV3];
        end
    else
        if grat
            MTGratdP=[MTGratdP,dprimeGratMT];
        end
        if dots
            MTDotsdP=[MTDotsdP,dprimeDotsMT];
        end
    end

        
 end

clear dprimeGratV3 dprimeGratMT dprimeDotsV3 dprimeDotsMT
 end
end
%%
% figure
% scatter(goodch(V3range),dprimeGrat(V3range),'r')
% hold on
% scatter(goodch(V3range),dprimeDots(V3range),'b')
% legend(['grating'],['dots'])
% title('V3')
% figure
% scatter(goodch(V3range),dprimeGrat(V3range),'r')
% hold on
% scatter(goodch(V3range),dprimeDots(V3range),'b')
% legend(['grating'],['dots'])
% title('V3')
% figure
% scatter(dprimeDots(MTrange),dprimeGrat(MTrange),'r')
% hold on
% scatter(dprimeDots(V3range),dprimeGrat(V3range),'b')
% legend(['MT'],['V3'])
% xlabel('directional selectivity index - dots')
% ylabel('directional selectivity index - grating')
% %scatter dp grating versus dots for each area
% %subtract spontaneous activity when computing dprime
% figure
% scatter(goodch,dprimeGrat,'r')
% hold on
% scatter(goodch,dprimeDots,'b')
% legend(['grating'],['dots'])
% 
% figure
% ROC_data = roc_curve(dprimeDots(V3range),dprimeDots(MTrange));
% figure
% ROC_data2 = roc_curve(dprimeGrat(MTrange),dprimeGrat(V3range));
% clearvars -except V3GratdP MTGratdP V3DotsdP MTDotsdP



% figure
% histogram(V3GratdP(V3GratdP>=DIThresh),50)
% hold on
% histogram(V3DotsdP(V3DotsdP>=DIThresh),50)
% legend('grat','dots')
% title('V3')
% xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')
% 
% figure
% histogram(MTGratdP(MTGratdP>=DIThresh),20)
% hold on
% histogram(MTDotsdP(MTDotsdP>=DIThresh),20)
% legend('grat','dots')
% title('MT')
% xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')
for di=1:length(DIThresh)
    thresh=DIThresh(di);
    nDmt(di)=length(MTDotsdP(MTDotsdP<=thresh))/length(MTDotsdP);
    nGmt(di)=length(MTGratdP(MTGratdP<=thresh))/length(MTGratdP);
    nDv3(di)=length(V3DotsdP(V3DotsdP<=thresh))/length(V3DotsdP);
    nGv3(di)=length(V3GratdP(V3GratdP<=thresh))/length(V3GratdP);
end

figure
plot(DIThresh,nDmt)
hold on
plot(DIThresh,nDv3)
legend('MT','v3')
title(['Dots - coherence ',num2str(coherence)])
xlabel('threshold')
ylabel('ratio of neurons with selectivity <= threshold')
figure
plot(DIThresh,nGmt)
hold on
plot(DIThresh,nGv3)
legend('MT','v3')
title(['Grating - contrast ',num2str(contrast)])
xlabel('threshold')
ylabel('ratio of neurons with selectivity <= threshold')

figure
histogram(MTDotsdP(MTDotsdP<=10),50)
hold on
histogram(V3DotsdP(V3DotsdP<=10),50)
legend('MT','V3')
title(['Dots, coherence = ',num2str(coherence)])
xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')
figure
histogram(MTGratdP(MTGratdP<=10),50)
hold on
histogram(V3GratdP(V3GratdP<=10),50)
legend('MT','V3')
title(['grat, contrast = ',num2str(contrast)])
xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')

% figure
% histogram(MTDotsdP(MTDotsdP>DIThresh & MTDotsdP<=10),50)
% hold on
% histogram(V3DotsdP(V3DotsdP>DIThresh & V3DotsdP<=10),50)
% legend(['MT , ',num2str(length(MTDotsdP(MTDotsdP>DIThresh & MTDotsdP<=10))*100/length(MTDotsdP(MTDotsdP<=10))),'% > thresh']...
%     ,['V3 ,',num2str(length(V3DotsdP(V3DotsdP>DIThresh & V3DotsdP<=10))*100/length(V3DotsdP(V3DotsdP<=10))),'% > thresh'])
% title(['Dots, coherence = ',num2str(coherence), ', thresh = ', num2str(DIThresh)])
% xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')
% figure
% histogram(MTGratdP(MTGratdP>DIThresh & MTGratdP<=10),50)
% hold on
% histogram(V3GratdP(V3GratdP>DIThresh & V3GratdP<=10),50)
% legend(['MT , ',num2str(length(MTGratdP(MTGratdP>DIThresh & MTGratdP<=10))*100/length(MTGratdP(MTGratdP<=10))),'% > thresh']...
%     ,['V3 ,',num2str(length(V3GratdP(V3GratdP>DIThresh & V3GratdP<=10))*100/length(V3GratdP(V3GratdP<=10))),'% > thresh'])
% title(['grat, contrast = ',num2str(contrast), ', thresh = ', num2str(DIThresh)])
% xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')


figure
ROC_data = roc_curve(V3DotsdP,MTDotsdP);
annotation('textbox',[.83 .5 .1 .2],'String','Can responses to dots distinguish V3 from MT?','EdgeColor','none')
figure
ROC_data2 = roc_curve(MTGratdP,V3GratdP);
annotation('textbox',[0.83 .5 .1 .2],'String','Can responses to grating distinguish V3 from MT?','EdgeColor','none')
figure
ROC_data3 = roc_curve(V3DotsdP,V3GratdP);
annotation('textbox',[.83 .5 .1 .2],'String','Can V3 responses distinguish dots stim from grating stim?','EdgeColor','none')
figure
ROC_data4 = roc_curve(MTGratdP,MTDotsdP);
annotation('textbox',[.83 .5 .1 .2],'String','Can MT responses distinguish dots stim from grating stim?','EdgeColor','none')
