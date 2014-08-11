function [w,BETA,scoresPlot]=findWormMasses(sqPeaks,CMZ,knownMZ)

%Outputs
% w     Index in CMZ for best masses
% BETA  Betas for those masses
% scoresPlot   Image for you to use looking pretty
if nargin==2
    knownMZ = 542.5;
end
allPeaks=cat(2,sqPeaks{:}); %CMZ x scan
[~,knownMZidx]=min(abs(CMZ-knownMZ));
y=allPeaks(knownMZidx,:);

sz=size(sqPeaks);
figure,imshow(reshape(y,sz),[])
isWorm=reshape(mat2gray(imfilter(mat2gray(reshape(y,size(sqPeaks))),fspecial('gaussian'))),[],1);
[XL,~,XS,~,~,PCTVAR] = plsregress(zscore(allPeaks'),isWorm,10);

numPC=find(cumsum(PCTVAR(1,:))>.90,1,'first');

[XL,~,XS,~,BETA,PCTVAR] = plsregress(zscore(allPeaks'),isWorm,numPC);

%Figure out if worm is positive or negative scores
if sum(XS(isWorm>0.2))>0
    wormIsPos = 1;
else
    wormIsPos = 0;
end

%Generate composite scores plot
% by summing all components (looks nicer)
scoresPlot=reshape(sum(XS(:,1:numPC),2),sz);
% or by weighting based on the percent variance explained in each component
%scoresPlot=reshape(XS(:,1:numPC)*PCTVAR(1,:)',sz);

figure,imagesc(scoresPlot), axis equal
colormap('redblue')
title(['Scores Plot, PCs: 1-',num2str(numPC)])

%find best loadings (m/zs). cutoff = zscore > 1.5
BETA=BETA(2:end);
[BETA_sorted,idx]=sort(BETA);
if wormIsPos
    w=idx(zscore(BETA_sorted)>1.5);
else
    w=idx(zscore(BETA_sorted)<1.5);
end
CMZ(w); %show m/zs
BETA=BETA(w);