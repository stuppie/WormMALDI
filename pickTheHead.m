function varargout = pickTheHead(varargin)
% PICKTHEHEAD MATLAB code for pickTheHead.fig
%      PICKTHEHEAD, by itself, creates a new PICKTHEHEAD or raises the existing
%      singleton*.
%
%      H = PICKTHEHEAD returns the handle to a new PICKTHEHEAD or the handle to
%      the existing singleton*.
%
%      PICKTHEHEAD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PICKTHEHEAD.M with the given input arguments.
%
%      PICKTHEHEAD('Property','Value',...) creates a new PICKTHEHEAD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pickTheHead_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pickTheHead_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pickTheHead

% Last Modified by GUIDE v2.5 22-Aug-2013 18:08:54

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
wormColors=varargin{1};
out=varargin{2};
sample=varargin{3};
names=varargin{4};

% get the raw image
imshow(imresize(imread(names{sample}),.25),'parent',handles.axes2)

imshow(wormColors,'Parent',handles.axes1)
[handles.X,handles.Y] = ginput(1);
pickedIdx=sub2ind([size(wormColors,1)/4,size(wormColors,2)/4],round(handles.Y/4),round(handles.X/4));
segment=find(~cellfun(@isempty,cellfun(@(x) find(round(x)==round(pickedIdx)),out,'uniformoutput',0)));
if isempty(segment)
    segment=0;
end
handles.segment=segment;
% Update handles structure
guidata(hObject, handles);
%assignin('caller','x',handles.X);
%assignin('caller','y',handles.Y);
%close(handles.figure1)

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
