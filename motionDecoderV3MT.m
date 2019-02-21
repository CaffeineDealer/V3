%direction decoder
% firing rates from V3 and firing rates from MT which pool can decode better
% inpu firing rate of the neuron
% output decision of the direction of the stimulus
%training data firing rates for rightward motion and firing rates for
%leftward motion (labels known) (you can use trials for training and trials for testing)
dotsnames={'slu048b','slu047b','slu046c','slu045c','slu044b','slu022b','slu017b',...
    'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c'};
gratnames={'slu048d','slu047e','slu046e','slu045e','slu044d','slu022d','slu017c',...
    'ytu326a','ytu331b','ytu332d','ytu334b','ytu336a'};
pathD='C:\research\data\SuperTuneSpkTrains\';
pathG='C:\research\data\PlaidSpkTrains\';

Fs=10000;
timeWin=(0.05*Fs:0.4*Fs);
contrast=100;%40;%100
coherence=100;%55;%100

V3gratlabel1=[]; V3dotslabel1=[]; V3gratlabel2=[]; V3dotslabel2=[];
MTgratlabel1=[]; MTdotslabel1=[]; MTgratlabel2=[]; MTdotslabel2=[];

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
        dotspktrain=load([pathD,dotsnames{1,j},num2str(ch),num2str(u),'spktrain.mat']);
        dotfiring=sum(dotspktrain.spktrain(timeWin,:,:,:,:,:,:),1)*Fs/length(timeWin);%size(gratspktrain.spktrain,1);
        dotfiring=mean(dotfiring,5); %[dtim,ddir,motTyp,dpos,dtrial,dsiz,dcoh]
        [maxfirD,ID]=max(dotfiring(:));
        [dtim,ddir,motTyp,dpos,dtrial,dsiz,dcoh]= ind2sub(size(dotfiring),ID);
        dirfirdot=squeeze(dotfiring(:,:,motTyp,dpos,dtrial,dsiz,coherenceidx));
        dotslabel1=dirfirdot(1);
        dotslabel2=dirfirdot(1+size(dotfiring,2)/2);
        
        gratspktrain=load([pathG,gratnames{1,j},num2str(ch),num2str(u),'spktrain.mat']);
        gratingfiring=sum(gratspktrain.spktrain(timeWin,:,:,:,:,:,:),1)*Fs/length(timeWin);%size(gratspktrain.spktrain,1);
        gratingfiring=mean(gratingfiring,7); %[gtime,gdir,gspd,gpos,gsiz,gcont]
        [maxfirG,IG]=max(gratingfiring(:));
        [gtime,gdir,gspd,gpos,gsiz,gcont]= ind2sub(size(gratingfiring),IG); 
        dirfirgrat=squeeze(gratingfiring(:,:,gspd,gpos,gsiz,contrastidx));
        gratlabel1=dirfirgrat(1);
        gratlabel2=dirfirgrat(1+size(gratingfiring,2)/2);

        
        if ci<=size(V3units,1)
            V3gratlabel1=[V3gratlabel1,gratlabel1];
            V3dotslabel1=[V3dotslabel1,dotslabel1];
            V3gratlabel2=[V3gratlabel2,gratlabel2];
            V3dotslabel2=[V3dotslabel2,dotslabel2];
        else
            MTgratlabel1=[MTgratlabel1,gratlabel1];
            MTdotslabel1=[MTdotslabel1,dotslabel1];
            MTgratlabel2=[MTgratlabel2,gratlabel2];
            MTdotslabel2=[MTdotslabel2,dotslabel2];
        end

    end
end
        % get the firingrates for specific stimuli from Nfiles
        %V3training firingrates and labels
v3Gtraining=[V3gratlabel1' ones(length(V3gratlabel1),1); V3gratlabel2' 2*ones(length(V3gratlabel2),1)];
v3Gtraining = v3Gtraining(randperm(size(v3Gtraining,1)),:);

v3Dtraining=[V3dotslabel1' ones(length(V3dotslabel1),1); V3dotslabel2' 2*ones(length(V3dotslabel2),1)];
v3Dtraining = v3Dtraining(randperm(size(v3Dtraining,1)),:);
        
%         load(['C:\research\data\SuperTuneFiringMatrix\',name,num2str(ch),num2str(u),'firingMat']);
        %MTtraining firingrates and labels
mtGtraining=[MTgratlabel1' ones(length(MTgratlabel1),1); MTgratlabel2' 2*ones(length(MTgratlabel2),1)];
mtGtraining = mtGtraining(randperm(size(mtGtraining,1)),:);

mtDtraining=[MTdotslabel1' ones(length(MTdotslabel1),1); MTdotslabel2' 2*ones(length(MTdotslabel2),1)];
mtDtraining = mtDtraining(randperm(size(mtDtraining,1)),:);       
                % fitting a decoder [Mdl,FitInfo] = fitclinear(X,Ystats)
Xv3g= v3Gtraining(:,1);  Yv3g= v3Gtraining(:,2); 
Xv3d= v3Dtraining(:,1);  Yv3d= v3Dtraining(:,2); 
Xmtg= mtGtraining(:,1);  Ymtg= mtGtraining(:,2); 
Xmtd= mtDtraining(:,1);  Ymtd= mtDtraining(:,2); 
%[Mdl,FitInfo] = fitclinear(X,Y);
Lambda = logspace(-6,-0.5,11);

CVMdlv3g = fitclinear(Xv3g,Yv3g,'KFold',5,'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8);
cev3g = kfoldLoss(CVMdlv3g);

CVMdlv3d = fitclinear(Xv3d,Yv3d,'KFold',5,'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8);
cev3d = kfoldLoss(CVMdlv3d);

CVMdlmtg = fitclinear(Xmtg,Ymtg,'KFold',5,'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8);
cemtg = kfoldLoss(CVMdlmtg);

CVMdlmtd = fitclinear(Xmtd,Ymtd,'KFold',5,'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8);
cemtd = kfoldLoss(CVMdlmtd);


figure;plot(log10(Lambda),log10(cev3g))
hold on 
plot(log10(Lambda),log10(cev3d))
hold on 
plot(log10(Lambda),log10(cemtg))
hold on 
plot(log10(Lambda),log10(cemtd))
ylabel('log_{10} classification error')
xlabel('log_{10} Lambda')
legend('v3g','v3d','mtg','mtd')