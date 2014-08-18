%% Install proteowizard
% http://proteowizard.sourceforge.net/downloads.shtml
% Make sure you include vendor reader support (windows only).
% Add the folder containing msconvert to your windows path. E.g. C:\Program Files (x86)\ProteoWizard\ProteoWizard 3.0.6452

%% Select all the raw files in the current directory
d=dir('*.raw');
%d=dir('*.mzXML');
d={d.name};

%% or manually pick them
d=uigetfile('*.raw','multiselect','on');
if ~iscell(d)
    d={d};
end
%% Convert to mzxml and get list of spot IDs
spotIDs=[];
for i=1:length(d)
    system(['msconvert ' d{i} ' --mzXML --32 --mz32']);
    spotIDs{i} = readSpotIDs(d{i});
end
%% Read in mzXML raw data, generate common m/z
scanIdx=[];
AllPeaks=[];
for i=1:length(d)
    [~,fname]=fileparts(d{i});
    msStruct=mzxmlread([fname '.mzXML']);
    [peaks,time] = mzxml2peaks(msStruct);
    %remove all m/z > 1000 (not the intensity!, the m/z range)
    peaks=cellfun(@(x) x(x(:,1)<1000,:), peaks, 'uniformoutput',false);
    AllPeaks = [AllPeaks; peaks];
    scanIdx = [scanIdx; repmat(i,length(peaks),1)];
end
[CMZ, AllPeaks] = mspalign(AllPeaks);
%for each scan, fill in zeros where the mass is in the CMZ but not in the scan
for i=1:length(AllPeaks)
    AllPeaks{i}=accumarray([find(ismember(CMZ,AllPeaks{i}(:,1)))';(1:length(CMZ))'],[AllPeaks{i}(:,2);zeros(length(CMZ),1)],[],@sum);
end
%split back into groups based on the sample number
peaks=[];
for i=1:length(d)
    peaks{i}=AllPeaks(scanIdx==i);
end

%% Do work
worm=[];
for sample=1:length(d)
    %% Sum overlapping scans. Sum with a window size of 2.
    % Using middle of laser, for each pixel, sum the scans that overlap that pixel.
    % If spotDiameter = 100, overlap = 25. Window size = spotDiameter / (overlap * 2)

    % sqPeaks = 2D cell array (height x width of image). Each cell contains a vector [1xN],
    %   length(N) = length(CMZ)

    windowSize = 2; 
    spotID=spotIDs{sample}(:,[3,4]);
    spotID(:,1)=grp2idx(spotID(:,1));
    spotID(:,2)=grp2idx(spotID(:,2));
    sqPeaks=cell(max(spotID));
    sqPeaks=cellfun(@(x) single(zeros(length(CMZ),1)),sqPeaks,'uniformoutput',0); %
    peak=peaks{sample};
    for i=1:length(spotID)
        xInd=spotID(i,1)-windowSize:spotID(i,1)+windowSize;
        xInd=xInd(xInd>0 & xInd<=max(spotID(:,1)));
        yInd=spotID(i,2)-windowSize:spotID(i,2)+windowSize;
        yInd=yInd(yInd>0 & yInd<=max(spotID(:,2)));
        sqPeaks(xInd,yInd)=cellfun(@(x) sum([x,peak{i}],2),sqPeaks(xInd,yInd),'uni',0);
    end
    allPeaks=cat(2,sqPeaks{:}); %CMZ x scan
    sz=size(sqPeaks);
    %% plot mean of each scan
    % **** Decide if you want to normalize the total intensity of each pixel across all masses*****
    % figure,imshow(reshape(mean(allPeaks,1),sz),[])
    
    %% Pull out image of a specific mass
    [~,knownMZidx]=min(abs(CMZ-542.5));
    y=allPeaks(knownMZidx,:);
    figure,imshow(reshape(y,sz),[])
    %% Find worm masses
    [w,beta,scoresPlot]=findWormMasses(sqPeaks,CMZ,542.5);
    s=allPeaks(w,:);
    s=s'*beta; %weight based on betas
    figure,imshow(reshape(s,sz),[])
    CMZ(w)
    %% Divide worm into three segments
    [segmentScan, wormRGB]=segmentWorms(scoresPlot,1);
    %% Pick the head
    segment=0;
    while ~(segment==1 || segment==3)
        segment=pickTheHead(segmentScan, wormRGB, [num2str(sample,'%05d'),'.jpg']);
    end
    % If the head is segment 3, make the head segment 1. And make the head 'red' in the RGB image
    if segment == 3
        segmentScan([1,3])=segmentScan([3,1]);
        wormRGB(:,:,[1,3])=wormRGB(:,:,[3,1]);
    end
    worm(sample).segment=segment; % don't really need this anymore
    worm(sample).wormRGB=wormRGB;
    worm(sample).segmentScan=segmentScan;
    worm(sample).size=sz;
    worm(sample).allPeaks=allPeaks;
    worm(sample).sqPeaks=sqPeaks;
    close all
end
%% Make Y-vect for PLS specific to head/mid/tail
y=[];
for sample = 1:length(d)
    ySample=reshape(zeros(worm(sample).size),[],1); 
    ySample=ySample-5;
    ySample(worm(sample).segmentScan{1})=-10;
    ySample(worm(sample).segmentScan{2})=-10;
    ySample(worm(sample).segmentScan{3})=10;
    y=[y;ySample];
    if sample==1
    worm(sample).segs=[1 size(y,1)];  
    else
    worm(sample).segs=[worm(sample-1).segs(2)+1 size(y,1)];  
    end
end

figure,imshow(reshape(ySample,worm(sample).size),[]), title('Y-vect for only the last sample')

% *** Or specific to just worm (no segments) ***
% y(y==-10)=10;
%% Find masses related to sections
X = cat(2,worm.allPeaks)'; % rows are scans, columns are m/z from CMZ. Same number of rows as len(y)
[~,~,~,~,~,PCTVAR] = plsregress(zscore(X),y,10);
numPC=find(cumsum(PCTVAR(1,:))>.90,1,'first');
[XL,~,XS,~,BETA,PCTVAR] = plsregress(zscore(X),y,numPC);

%% plot scores plot for sample 1
Sample=2;

XS1=sum(XS(worm(Sample).segs(1):worm(Sample).segs(2),:),2);

scoresPlot=reshape(XS1,worm(Sample).size);
figure,imagesc(scoresPlot), axis equal
colormap('redblue')

%% find best m/zs using betas
BETA=BETA(2:end);
[BETA_sorted,idx]=sort(BETA);
w=idx(zscore(BETA_sorted)>2);
BETAW=BETA(w);
CMZ(w)
%make image of those masses. Pull out masses by largest Betas, weight image by betas
s=cellfun(@(x) (x(w,:))'*BETAW, {worm.allPeaks}, 'uni', 0);
%% Make image of those masses in sample S
sample=3;
imshow(reshape(s{sample},worm(sample).size),[])
colormap(jet)
