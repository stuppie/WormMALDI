function [out, worm]=segmentWorms(in,showPlot)

% in         Grayscale image (2D matrix)
% showPlot   1 or 0

% outputs:
% out     The linear indices for the 3 segments
% worm    RGB image of the 3 segments (not in the original size). For viewing purposes only

if exist('showPlot','var')==0
    showPlot=1;
end
%in=scoresPlot;
in(isnan(in))=0;
mz=imresize(mat2gray(in),4);
mz2=mat2gray(imfilter(mz,fspecial('average',[8 8])));
%%
%mz3=mat2gray(graydist(mz2,mz2<=graythresh(mz2)));

bw1=mz2>=graythresh(mz2)*1.5;
bw2=imclose(bw1,strel('disk',6));
bw2=imfill(bw2,'holes');
bw2=bwmorph(bw2,'majority');
A=bwconncomp(bw2); %keep only the largest object
[~,idx]=max(cellfun(@length,A.PixelIdxList));
bw2(cat(1,A.PixelIdxList{setdiff([1:A.NumObjects],idx)}))=0;

%imoverlay(mz2,bwperim(bw2),[1 0 0]);
%%
bwSkel=skeleton(bw2);
cutoff=30;
keepGoing=1;
while cutoff<max(max(bwSkel)) && keepGoing
    bp=find(bwmorph(bwmorph(bwSkel>cutoff,'skel',Inf),'branchpoints'));
    if bp~=0
        cutoff=cutoff+1;
    else
        keepGoing=0;
    end
end
bw3=bwmorph(skeleton(bw2)>cutoff,'skel',Inf);
if showPlot, figure,imshow(imoverlay(mz2,logical(bw3+bwperim(bw2)),[1 0 0])), end

%% for removing pieces of the skeleton
if 0
    [skr,rad] = skeleton(bw2);
    imagesc(skr);
    colormap jet
    axis image off
end

%% Trace along the skeleton
[endpointsx,endpointsy]=find(bwmorph(bw3,'endpoints'));
pos=[endpointsx(1),endpointsy(1)];
S=pos;
bw3New=bw3;
bw3New(pos(1),pos(2))=0;
while ~(pos(1)==endpointsx(2) && pos(2)==endpointsy(2))
    local=bw3New(pos(1)-1:pos(1)+1,pos(2)-1:pos(2)+1);
    dx=sum(sum(local.*[-1,-1,-1;0,0,0;1,1,1]));
    dy=sum(sum(local.*[-1,0,1;-1,0,1;-1,0,1]));
    pos=[pos(1)+dx,pos(2)+dy];
    S=[S;pos];
    bw3New(pos(1),pos(2))=0;
end

%% Extend skeleton to endpoints
perim=bwperim(bw2);
perimPixels=bwdist(perim);
perimPixels = perimPixels<=1;

S = extendSkeleton(S, perimPixels);
S = flipud(S);
S = extendSkeleton(S, perimPixels);
S = flipud(S);

if showPlot, figure,imshow(perim), hold, plot(S(:,2),S(:,1),'x'), end

%% Calculate length of skeleton and divide into three
dist=sqrt((S(2:end,1)-S(1:end-1,1)).^2+ ...
    (S(2:end,2)-S(1:end-1,2)).^2);
dist=cumsum(dist);
cutoffs=[find(dist>(max(dist)/3),1,'first'), find(dist>(max(dist)/3)*2,1,'first')];

if showPlot
figure, imshow(perim), hold, plot(S(:,2),S(:,1),'y-')
plot(S(cutoffs,2), S(cutoffs,1), 'ro')
end

%% Segment

head=S(1:cutoffs(1),:);
a=zeros(size(bw3));
a(sub2ind(size(a),head(:,1),head(:,2)))=1;
head=(1-mat2gray(bwdist(a)));

mid=S(cutoffs(1):cutoffs(2),:);
a=zeros(size(bw3));
a(sub2ind(size(a),mid(:,1),mid(:,2)))=1;
mid=(1-mat2gray(bwdist(a)));

tail=S(cutoffs(2):end,:);
a=zeros(size(bw3));
a(sub2ind(size(a),tail(:,1),tail(:,2)))=1;
tail=(1-mat2gray(bwdist(a)));

worm=cat(3,head,mid,tail);
%normalize histograms
worm(:,:,2) = histeq(worm(:,:,2), imhist(worm(:,:,1)));
worm(:,:,3) = histeq(worm(:,:,3), imhist(worm(:,:,1)));
worm(worm==0)=1;
worm = imfilter(worm,fspecial('average',[6 6]));

for i=1:size(worm,1)
    for j=1:size(worm,2)
        [~,idx]=max(squeeze(worm(i,j,:)));
        worm(i,j,:)=NaN;
        worm(i,j,idx)=1;
    end
end
worm(:,:,1)=worm(:,:,1).*bw2;
worm(:,:,2)=worm(:,:,2).*bw2;
worm(:,:,3)=worm(:,:,3).*bw2;
if showPlot, figure,imshow(worm), end

%% output
head=abs(imresize(worm(:,:,1),0.25,'nearest'))>0;
out{1}=find(head);
mid=abs(imresize(worm(:,:,2),0.25,'nearest'))>0;
out{2}=find(mid);
tail=abs(imresize(worm(:,:,3),0.25,'nearest'))>0;
out{3}=find(abs(tail));
end
