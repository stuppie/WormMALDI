function spotIDs = readSpotIDs(filename)
% Use msconvert to convert the raw file into text format. The text format contains the
% spot coordinates which can be parsed out, while the mzXML doesn't.
% filename: path to raw file
% spotIDs: [number of scans x 4] matrix

system(['msconvert ' filename ' --text --filter "threshold count 1 most-intense']);
%%
[~,fname]=fileparts(filename);
fid=fopen([fname '.txt'],'r');
line=fgetl(fid);
scans=[];
spotIDs=cell(0);

while 1
    line=fgetl(fid);
    if ~ischar(line), break, end
    while isempty(regexp(line,'scan=(\d{1,8})','once'))
        line=fgetl(fid);
        if ~ischar(line), break, end
    end
    if ~ischar(line), break, end
    index=regexp(line,'scan=(\d{1,8})','tokens');
    scans=[scans,str2num(index{1}{1})];
    line=fgetl(fid);
    if isempty(strfind(line,'spotID:'))
        error('Something went wrong')
    end
    spotID=regexp(line,'([-]?\d{1,2})x([-]?\d{1,2}),(\d+)x(\d+)','tokens');
    spotID=spotID{1};
    spotIDs{end+1}=cellfun(@str2num,spotID);
end
fclose(fid);
spotIDs=cat(1,spotIDs{:});
%mins=min(spotIDs);
%spotIDs(:,1:2)=bsxfun(@minus,spotIDs(:,1:2),mins(1,2)-1);
end
