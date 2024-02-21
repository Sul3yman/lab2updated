function varargout = Interface(varargin)
% INTERFACE MATLAB code for Interface.fig
%      INTERFACE, by itself, creates a new INTERFACE or raises the existing
%      singleton*.
%
%      H = INTERFACE returns the handle to a new INTERFACE or the handle to
%      the existing singleton*.
%
%      INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERFACE.M with the given input arguments.
%
%      INTERFACE('Property','Value',...) creates a new INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Interface

% Last Modified by GUIDE v2.5 13-Jan-2019 03:37:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Interface_OpeningFcn, ...
                   'gui_OutputFcn',  @Interface_OutputFcn, ...
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


% --- Executes just before Interface is made visible.
function Interface_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Interface (see VARARGIN)

clc

% Choose default command line output for Interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function tx_text_Callback(hObject, eventdata, handles)
% hObject    handle to tx_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Hints: get(hObject,'String') returns contents of tx_text as text
%        str2double(get(hObject,'String')) returns contents of tx_text as a double


% --- Executes during object creation, after setting all properties.
function tx_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_transmit.
function pushbutton_transmit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_transmit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_flag=1;

%Get data to transmit
tx_infobits=handles.data.tx_infobits;          


%Transmit 
tx_samples=transmit(tx_infobits,plot_flag);


%Send samples over Channel
noise_flag=get(handles.noise_yes,'Value');
channel=get(handles.channel1,'Value');
rx_samples=simchannel(tx_samples, noise_flag,channel);
 
%Receive 
rx_infobits=receive(rx_samples,plot_flag);

%convert into string of bits
rx_infobits_str=num2str(rx_infobits);
rx_infobits_str(isspace(rx_infobits_str)) = '';

% store bit string in data, and update tx_bits text box
handles.data.rx_infobits_str=rx_infobits_str;
guidata(hObject, handles);
set(handles.rx_bits,'String',rx_infobits_str);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.tx_bits,'String','');
set(handles.tx_text,'String','');
set(handles.rx_bits,'String','');
set(handles.rx_text,'String','');
set(handles.BER,'String','');



% --- Executes on button press in pushbutton_exit.
function pushbutton_exit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all;
clear all



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to tx_bits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tx_bits as text
%        str2double(get(hObject,'String')) returns contents of tx_bits as a double


% --- Executes during object creation, after setting all properties.
function tx_bits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_bits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_convert_to_text.
function pushbutton_convert_to_text_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_convert_to_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rx_infobits_str=handles.data.rx_infobits_str;

rx_infobits_str_tmp=reshape(rx_infobits_str,[],8);
rx_string=char(bin2dec(rx_infobits_str_tmp)');

set(handles.rx_text,'String',rx_string);


% --- Executes on button press in pushbutton_convert_to_bits.
function pushbutton_convert_to_bits_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_convert_to_bits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tx_string=get(handles.tx_text,'String');
tx_infobits_str=reshape(dec2bin(tx_string,8),1,[]);
tx_infobits=tx_infobits_str-'0';

%store data and update text box
handles.data.tx_infobits=tx_infobits;
guidata(hObject, handles);
set(handles.tx_bits,'String',tx_infobits_str);


% --- Executes on button press in pushbutton_Test_BER.
function pushbutton_Test_BER_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Test_BER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_flag=0;
%Generate test sequence of bits
BER_target=5*10E-3;
nRealizations=1/BER_target*1000;
%tx_infobits=randi(1,1,nRealizations)-1; %Error : corrected in line below
tx_infobits=randi(2,1,nRealizations)-1;

%Transmit and receive test sequence of bits
tx_samples=transmit(tx_infobits,plot_flag);

noise_flag=get(handles.noise_yes,'Value');
channel=get(handles.channel1,'Value');
rx_samples=simchannel(tx_samples, noise_flag,channel);
  
rx_infobits=receive(rx_samples,plot_flag);

%Calcualt bit error rate (BER)
BER=    sum(rx_infobits ~= tx_infobits)/length(tx_infobits);
set(handles.BER,'String',num2str(BER));
