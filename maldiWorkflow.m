%% Install proteowizard
% http://proteowizard.sourceforge.net/downloads.shtml
% Make sure you include vendor reader support (windows only).
% Add the folder containing msconvert to your windows path. E.g. C:\Program Files (x86)\ProteoWizard\ProteoWizard 3.0.6452

%Required files: 

%% Select all the raw files in the current directory
d=dir('*.raw');
d={d.name};

%% or manually pick them
d=uigetfile('*.raw','multiselect','on')
if ~iscell(d)
    d={d}
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
    %remove peaks with m/z greater than 1000
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

%% Sum overlapping scans. Sum with a window size of 2.
% Using middle of laser, for each pixel, sum the scans that overlap that pixel.
% If spotDiameter = 100, overlap = 25. Window size = spotDiameter / (overlap * 2)

% sqPeaks = 2D cell array (height x width of image). Each cell contains a vector [1xN],
%   length(N) = length(CMZ)


%% *** Change to change sample number
sample=1;
%% *** Build square array by averaging scans based on spot size
windowSize = 2;
spotID=spotIDs{sample}(:,[3,4]);
spotID(:,1)=grp2idx(spotID(:,1));
spotID(:,2)=grp2idx(spotID(:,2));
sqPeaks=cell(max(spotID));
sqPeaks=cellfun(@(x) single(zeros(length(CMZ),1)),sqPeaks,'uni',0); %
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
% Decide if you want to normalize the total intensity of each pixel across all masses
imshow(reshape(mean(allPeaks,1),sz),[])

%% pull out image of a mass
[~,knownMZidx]=min(abs(CMZ-542.5));
y=allPeaks(knownMZidx,:);
figure,imshow(reshape(y,sz),[])
%% Find worm masses
[w,beta,scoresPlot]=findWormMasses(sqPeaks,CMZ,542.5);
s=allPeaks(w,:);
s=s'*beta; %weight based on betas
figure,imshow(reshape(s,sz),[])
CMZ(w)
%% Divide worm in three
[out, worm]=segmentWorms(scoresPlot,1);
%% Pick the head
idx=pickTheHead(worm,out,1,{'00001.jpg'})
%% Make Y-vect for PLS specific to head/mid/tail
y=zeros(sz);
y=y-5;
y(out{1})=-10;
y(out{2})=-10;
y(out{3})=10;
figure,imshow(y,[])

%% Find masses related to sections
[XL,~,XS,~,~,PCTVAR] = plsregress(zscore(allPeaks'),reshape(y,[],1),10);
numPC=find(cumsum(PCTVAR(1,:))>.90,1,'first');
[XL,~,XS,~,BETA,PCTVAR] = plsregress(zscore(allPeaks'),reshape(y,[],1),numPC);
scoresPlot=reshape(sum(XS(:,1:numPC),2),sz);
figure,imagesc(scoresPlot), axis equal
colormap('redblue')
BETA=BETA(2:end);
[BETA_sorted,idx]=sort(BETA);
w=idx(zscore(BETA_sorted)>4);
BETAW=BETA(w);
%make image of those masses
s=allPeaks(w,:);
s=s'*BETAW; %weight based on betas
figure,imshow(reshape(s,sz),[])
CMZ(w)