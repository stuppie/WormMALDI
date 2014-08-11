function varargout = pickTheHead(varargin)
% ***********
% Required arguments:
% segmentScans   cell array containing the linear indices for the scans that form each
%                segment
% worms          RGB image of the 3 segments (use the output from segmentWorms, it is 4x 
%                larger) than the values in segmentScans suggest it should be.
% name           path to image of worm. Used to figure out which segment is the head
% 
% Output: the segment that you clicked on
% ***********

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pickTheHead_OpeningFcn, ...
                   'gui_OutputFcn',  @pickTheHead_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pickTheHead is made visible.
function pickTheHead_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pickTheHead (see VARARGIN)

% Choose default command line output for pickTheHead
handles.output = hObject;

% Parse input args
segmentScans=varargin{1};
wormColors=varargin{2};
name=varargin{3};

% Read in the raw image of the worm
imshow(imresize(imread(name),.25),'parent',handles.axes2)

imshow(wormColors,'Parent',handles.axes1)
[handles.X,handles.Y] = ginput(1);
pickedIdx=sub2ind([size(wormColors,1)/4,size(wormColors,2)/4],round(handles.Y/4),round(handles.X/4));
segment=find(~cellfun(@isempty,cellfun(@(x) find(round(x)==round(pickedIdx)),segmentScans,'uniformoutput',0)));
if isempty(segment)
    segment=0;
end
handles.segment=segment;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pickTheHead wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pickTheHead_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.segment;
varargout{2} = handles.X;
varargout{3} = handles.Y;

close(handles.figure1)