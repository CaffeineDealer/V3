% v3 RF mapping

names={'slu023b'};%'slu022e', 'slu023b','slu047d', 'slu050b', 'slu055c', 'slu062b', 'slu060c'};%'slu039?','slu046a',
path='C:\research\data\SuperTuneSpkTrains\';
Fs=10000;
count=1;
count2=1;
probesize=10;
centroidV3=[];
sizeV3=[];
centroidMT=[];
sizeMT=[];
for j=1:length(names)
    load(['C:\research\V3 things\V3 categorized\',names{1,j}(1:end-1),'_V3categ.mat']);
    v3categ=sortrows(v3categ);
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
        bestdir=squeeze(max(firing.firing,[],1));
        bestdir=reshape(bestdir,[4 6]);
        RFsmooth = imgaussfilt(bestdir,0.7);

        
        matrix=RFsmooth/sum(RFsmooth(:));%/sum(RFsmooth(:));
        [m,n]=size(matrix);
        [I,J]=ndgrid(-m/2:m/(m-1):m/2,-n/2:n/(n-1):n/2);
        I=I*-1;
%         centroid=[dot(I(:),matrix(:)),  dot(J(:),matrix(:))];
%         centerx=ceil(centroid(1)+(m/2)*m/(m-1));
%         centery=ceil(centroid(2)+(n/2)*n/(n-1));
%         RFcenter=sqrt(centroid(1)^2+centroid(2)^2)*probesize;
%         mu=matrix(centerx,centery);
%         sigma=sqrt(sum((matrix(:)-mu).^2)/length(matrix(:)));
%         size1=2*sigma*probesize;
%         figure
%         imagesc(bestdir)
%         figure
%         imagesc(matrix)
        %title(['centroid ',num2str(centroid),',size ',num2str(size1)])
        [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(I,J,matrix);
        centroid2=sqrt(fitresult(5)^2+fitresult(6)^2);%*probesize;
        sigma2=sqrt(fitresult(3)^2+fitresult(4)^2);
        size2=2*sigma2;%*probesize; fix thissssssssssssssssssssssssssssssssssssssssssss
        figure
        imagesc(zfit)
%         
        if ci<=size(V3units,1)
            centroidV3=[centroidV3,centroid2];
            sizeV3=[sizeV3,size2];
        else
            centroidMT=[centroidMT,centroid2];
            sizeMT=[sizeMT,size2];
        end
    end
end
figure; scatter(centroidV3,sizeV3)
hold on
scatter(centroidMT,sizeMT)
legend('V3','MT')
