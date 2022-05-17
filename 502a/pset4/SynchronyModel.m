function varargout = SynchronyModel(varargin)
% SYNCHRONYMODEL MATLAB code for SynchronyModel.fig
%      SYNCHRONYMODEL, by itself, creates a new SYNCHRONYMODEL or raises the existing
%      singleton*.
%
%      H = SYNCHRONYMODEL returns the handle to a new SYNCHRONYMODEL or the handle to
%      the existing singleton*.
%
%      SYNCHRONYMODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYNCHRONYMODEL.M with the given input arguments.
%
%      SYNCHRONYMODEL('Property','Value',...) creates a new SYNCHRONYMODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SynchronyModel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SynchronyModel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SynchronyModel

% Last Modified by GUIDE v2.5 24-Mar-2016 13:34:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SynchronyModel_OpeningFcn, ...
    'gui_OutputFcn',  @SynchronyModel_OutputFcn, ...
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


% --- Executes just before SynchronyModel is made visible.
function SynchronyModel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SynchronyModel (see VARARGIN)

% Choose default command line output for SynchronyModel
handles.output = hObject;

%Setup the network(s)
handles.Net.Ne=800;
handles.Net.Ni=200;
handles.SimTime = 1000; %in ms

% Update handles structure
guidata(hObject, handles);

% Set the ranges of the sliders
set(handles.AmpSlider, 'Min', 0, 'Max', 1, 'Value', 0);
set(handles.text5, 'String', {'Oscillation', 'Amplitude', sprintf('%1.2f', 0)});
set(handles.FreqSlider, 'Min', 1, 'Max', 80, 'Value', 20);
set(handles.text6, 'String', sprintf('Frequency (%2.0f Hz)', 20));
set(handles.RelWeightSlider, 'Min', -100, 'Max', 100, 'Value', 0);
set(handles.text8, 'String', {'Relative', 'Weight', sprintf('%2.0f%%', 0)});
axis(handles.ControlledAxes, 'off');
axis(handles.CompeteAxes, 'off');
set(handles.ControlledAxes, 'XTickLabel', '', 'YTickLabel', '');
set(handles.CompeteAxes, 'XTickLabel', '', 'YTickLabel', '');

handles = SimulateNetworks(handles);

% UIWAIT makes SynchronyModel wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles = SimulateNetworks(handles),

samp_f = 1000;
tran_I_X_strength = 0;
tran_ovr_strength = str2double(get(handles.edit1, 'String'));
tran_E_E_strength = str2double(get(handles.edit2, 'String'));
tran_E_I_strength = str2double(get(handles.edit3, 'String'));
down_internal_strength = str2double(get(handles.edit4, 'String'));
rel_inp_strength = get(handles.RelWeightSlider, 'Value');
tran_sparse = 0.2;

compete_inp = randn(handles.Net.Ne + handles.Net.Ni, handles.SimTime);
control_inp = get(handles.AmpSlider, 'Value').*sin(repmat([0:(handles.SimTime-1)], [handles.Net.Ne + handles.Net.Ni 1])./samp_f*get(handles.FreqSlider, 'Value')*2*pi) + ...
    (1-get(handles.AmpSlider, 'Value')).*randn(handles.Net.Ne + handles.Net.Ni, handles.SimTime);

re = rand(handles.Net.Ne, 1);
ri = rand(handles.Net.Ni, 1);
cntl_a=[0.02*ones(handles.Net.Ne,1); 0.02+0.08*ri];
cntl_b=[0.2*ones(handles.Net.Ne,1); 0.25-0.05*ri];
cntl_c=[-65+15*re.^2; -65*ones(handles.Net.Ni,1)];
cntl_d=[8-6*re.^2; 2*ones(handles.Net.Ni,1)];
cntl_S=0.5*[0.5*rand(handles.Net.Ne+handles.Net.Ni, handles.Net.Ne), -rand(handles.Net.Ne+handles.Net.Ni,handles.Net.Ni)]; %%% Synaptic weight matrix - PCNN
cntl_v=-65*ones(handles.Net.Ne+handles.Net.Ni,1); % Initial values of v
cntl_u=cntl_b.*cntl_v; % Initial values of u
cntl_firings=[]; % spike timings

re = rand(handles.Net.Ne, 1);
ri = rand(handles.Net.Ni, 1);
cmpt_a=[0.02*ones(handles.Net.Ne,1); 0.02+0.08*ri];
cmpt_b=[0.2*ones(handles.Net.Ne,1); 0.25-0.05*ri];
cmpt_c=[-65+15*re.^2; -65*ones(handles.Net.Ni,1)];
cmpt_d=[8-6*re.^2; 2*ones(handles.Net.Ni,1)];
cmpt_S=0.5*[0.5*rand(handles.Net.Ne+handles.Net.Ni, handles.Net.Ne), -rand(handles.Net.Ne+handles.Net.Ni,handles.Net.Ni)]; %%% Synaptic weight matrix - PCNN
cmpt_v=-65*ones(handles.Net.Ne+handles.Net.Ni,1); % Initial values of v
cmpt_u=cmpt_b.*cmpt_v; % Initial values of u
cmpt_firings=[]; % spike timings

re = rand(handles.Net.Ne, 1);
ri = rand(handles.Net.Ni, 1);
down_a=[0.02*ones(handles.Net.Ne,1); 0.02+0.08*ri];
down_b=[0.2*ones(handles.Net.Ne,1); 0.25-0.05*ri];
down_c=[-65+15*re.^2; -65*ones(handles.Net.Ni,1)];
down_d=[8-6*re.^2; 2*ones(handles.Net.Ni,1)];
cntl_down_Stran = cat(2, cat(1, tran_E_E_strength*rand(handles.Net.Ne, handles.Net.Ne), tran_E_I_strength*rand(handles.Net.Ni, handles.Net.Ne)), ...
    -tran_I_X_strength*rand(handles.Net.Ne+handles.Net.Ni,handles.Net.Ni));
cntl_down_Stran(rand(size(cntl_down_Stran)) <= (1-tran_sparse)) = 0;
cntl_down_Stran = tran_ovr_strength*((100 - rel_inp_strength)./100)*cntl_down_Stran;

cmpt_down_Stran = cat(2, cat(1, tran_E_E_strength*rand(handles.Net.Ne, handles.Net.Ne), tran_E_I_strength*rand(handles.Net.Ni, handles.Net.Ne)), ...
    -tran_I_X_strength*rand(handles.Net.Ne+handles.Net.Ni,handles.Net.Ni));
cmpt_down_Stran(rand(size(cmpt_down_Stran)) <= (1-tran_sparse)) = 0;
cmpt_down_Stran = tran_ovr_strength*((100 + rel_inp_strength)./100)*cmpt_down_Stran;

down_S = down_internal_strength*[0.5*rand(handles.Net.Ne+handles.Net.Ni, handles.Net.Ne), -rand(handles.Net.Ne+handles.Net.Ni,handles.Net.Ni)]; %%% Synaptic weight matrix - PCNN
down_v=-65*ones(handles.Net.Ne+handles.Net.Ni,1); % Initial values of v
down_u=down_b.*down_v; % Initial values of u
down_firings=[]; % spike timings

%Simulates the networks based on inputs
for t = 1:handles.SimTime, % in ms
    %The control input network gets its input from outside drive
    cntl_I = [5*ones(handles.Net.Ne, 1); 2*ones(handles.Net.Ni, 1)].*control_inp(:, t);
    cntl_fired = find(cntl_v(:,t) >= 30); % indices of spikes
    cntl_firings = [cntl_firings; t+0*cntl_fired, cntl_fired];
    cntl_v(cntl_fired,t)=cntl_c(cntl_fired);
    cntl_u(cntl_fired) = cntl_u(cntl_fired) + cntl_d(cntl_fired);
    cntl_I = cntl_I + sum(cntl_S(:,cntl_fired), 2);
    
    cntl_v(:,t+1) = cntl_v(:,t) + 0.5 .* (0.04.*cntl_v(:,t).^2 + 5 .* cntl_v(:,t) + 140 - cntl_u + cntl_I); % step 0.5 ms
    cntl_v(:,t+1) = cntl_v(:,t+1) + 0.5.*(0.04.*cntl_v(:,t+1).^2 + 5 .* cntl_v(:,t+1) + 140 - cntl_u + cntl_I); % for numerical
    cntl_u = cntl_u + cntl_a.*(cntl_b.*cntl_v(:,t+1) - cntl_u); % stability
    
    %The competing input network gets its input from outside drive
    cmpt_I = [5*ones(handles.Net.Ne, 1); 2*ones(handles.Net.Ni, 1)].*compete_inp(:, t);
    cmpt_fired = find(cmpt_v(:,t) >= 30); % indices of spikes
    cmpt_firings = [cmpt_firings; t+0*cmpt_fired, cmpt_fired];
    cmpt_v(cmpt_fired,t)=cmpt_c(cmpt_fired);
    cmpt_u(cmpt_fired) = cmpt_u(cmpt_fired) + cmpt_d(cmpt_fired);
    cmpt_I = cmpt_I + sum(cmpt_S(:,cmpt_fired), 2);
    
    cmpt_v(:,t+1) = cmpt_v(:,t) + 0.5 .* (0.04.*cmpt_v(:,t).^2 + 5 .* cmpt_v(:,t) + 140 - cmpt_u + cmpt_I); % step 0.5 ms
    cmpt_v(:,t+1) = cmpt_v(:,t+1) + 0.5.*(0.04.*cmpt_v(:,t+1).^2 + 5 .* cmpt_v(:,t+1) + 140 - cmpt_u + cmpt_I); % for numerical
    cmpt_u = cmpt_u + cmpt_a.*(cmpt_b.*cmpt_v(:,t+1) - cmpt_u); % stability
    
    %The downstream network gets its input from both input networks
    down_I = sum(cntl_down_Stran(:,cntl_fired), 2) + sum(cmpt_down_Stran(:,cmpt_fired), 2) + [5*ones(handles.Net.Ne, 1); 2*ones(handles.Net.Ni, 1)].*randn(handles.Net.Ne + handles.Net.Ni, 1);
    down_fired = find(down_v(:,t) >= 30); % indices of spikes
    down_firings = [down_firings; t+0*down_fired, down_fired];
    down_v(down_fired, t) = down_c(down_fired);
    down_u(down_fired) = down_u(down_fired) + down_d(down_fired);
    down_I = down_I + sum(down_S(:,down_fired), 2);
    
    down_v(:,t+1) = down_v(:,t) + 0.5 .* (0.04.*down_v(:,t).^2 + 5 .* down_v(:,t) + 140 - down_u + down_I); % step 0.5 ms
    down_v(:,t+1) = down_v(:,t+1) + 0.5.*(0.04.*down_v(:,t+1).^2 + 5 .* down_v(:,t+1) + 140 - down_u + down_I); % for numerical
    down_u = down_u + down_a.*(down_b.*down_v(:,t+1) - down_u); % stability
end

%Save our variables to our handles structure for use elsewhere
handles.ControlActivity = cntl_firings;
handles.CompeteActivity = cmpt_firings;
handles.DownActivity = down_firings;

%Update figures
axes(handles.ControlledAxes); cla;
plot([0:(t-1)], control_inp(randsample(size(control_inp, 1), 1), 1:t), 'k-'); hold on;
e_ind = randsample(handles.Net.Ne, 30);
for i = 1:length(e_ind),
    temp_ind = ismember(cntl_firings(:, 2), e_ind(i));
    if ~isempty(temp_ind),
        plot(cntl_firings(temp_ind, 1), 3+i/5*ones(sum(temp_ind), 1), 'k.');
    end
end
i_ind = handles.Net.Ne + randsample(handles.Net.Ni, 12);
for i = 1:length(i_ind),
    temp_ind = ismember(cntl_firings(:, 2), i_ind(i));
    if ~isempty(temp_ind),
        plot(cntl_firings(temp_ind, 1), 9+i/5*ones(sum(temp_ind), 1), 'r.');
    end
end

axes(handles.CompeteAxes); cla;
plot([0:(t-1)], compete_inp(randsample(size(compete_inp, 1), 1), 1:t), 'k-'); hold on;
e_ind = randsample(handles.Net.Ne, 30);
for i = 1:length(e_ind),
    temp_ind = ismember(cmpt_firings(:, 2), e_ind(i));
    if ~isempty(temp_ind),
        plot(cmpt_firings(temp_ind, 1), 3+i/5*ones(sum(temp_ind), 1), 'k.');
    end
end
i_ind = handles.Net.Ne + randsample(handles.Net.Ni, 12);
for i = 1:length(i_ind),
    temp_ind = ismember(cmpt_firings(:, 2), i_ind(i));
    if ~isempty(temp_ind),
        plot(cmpt_firings(temp_ind, 1), 9+i/5*ones(sum(temp_ind), 1), 'r.');
    end
end

axes(handles.OutputAxes); cla;
plot(down_firings(:, 1), down_firings(:, 2), '.');
ind = (down_firings(:, 2) <= handles.Net.Ne);
plot(down_firings(ind, 1), down_firings(ind, 2), 'k.'); hold all;
plot(down_firings(~ind, 1), down_firings(~ind, 2), 'r.');
set(gca, 'YLim', [0 (handles.Net.Ne+handles.Net.Ni)], 'XLim', [0 (t-1)]);
drawnow;


% --- Outputs from this function are returned to the command line.
function varargout = SynchronyModel_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function AmpSlider_Callback(hObject, eventdata, handles)
% hObject    handle to AmpSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.text5, 'String', {'Oscillation', 'Amplitude', sprintf('%1.2f', get(hObject, 'Value'))});
handles = SimulateNetworks(handles);

% --- Executes during object creation, after setting all properties.
function AmpSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AmpSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function FreqSlider_Callback(hObject, eventdata, handles)
% hObject    handle to FreqSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.text6, 'String', sprintf('Frequency (%2.0f Hz)', get(hObject, 'Value')));
handles = SimulateNetworks(handles);

% --- Executes during object creation, after setting all properties.
function FreqSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in RunSimButton.
function RunSimButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunSimButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = SimulateNetworks(handles);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Run a bunch of simulations, estimating the transfer between regions
amp_list = [0:0.05:1];
for cur_amp_ind = 1:length(amp_list), 
    cur_amp = amp_list(cur_amp_ind);
    set(handles.AmpSlider, 'Value', cur_amp);
    set(handles.text5, 'String', {'Oscillation', 'Amplitude', sprintf('%1.2f', cur_amp)});
    
    handles = SimulateNetworks(handles);
    drawnow;
    
    %Look at correlation between firing activity
    cntl_act = zeros(handles.SimTime, 1);
    cmpt_act = zeros(handles.SimTime, 1);
    down_act = zeros(handles.SimTime, 1);
    
    for t = 1:handles.SimTime,
%         cntl_act(t) = sum(handles.ControlActivity(handles.ControlActivity(:, 2) <= handles.Net.Ne, 1) == t);
%         cmpt_act(t) = sum(handles.CompeteActivity(handles.CompeteActivity(:, 2) <= handles.Net.Ne, 1) == t);
%         down_act(t) = sum(handles.DownActivity(handles.DownActivity(:, 2) <= handles.Net.Ne, 1) == t);
        
        cntl_act(t) = sum(handles.ControlActivity(:, 1) == t);
        cmpt_act(t) = sum(handles.CompeteActivity(:, 1) == t);
        down_act(t) = sum(handles.DownActivity(:, 1) == t);
    end
    
    %[c,lags] = xcorr(cntl_act, down_act);
    net_corr(cur_amp_ind, 1) = corr(cntl_act, down_act);
    %[c,lags] = xcorr(cmpt_act, down_act);
    net_corr(cur_amp_ind, 2) = corr(cmpt_act, down_act);
end

figure;
plot(amp_list, net_corr);
legend({'Controlled Neural Population', 'Competing Neural Population'});
xlabel('Strength of Oscillation in Input to Controlled Network');
ylabel('Correlation between Input Networks and Output Network');
title('Effect of Oscillatory Synchrony on Information Transfer');


% --- Executes on slider movement.
function RelWeightSlider_Callback(hObject, eventdata, handles)
% hObject    handle to RelWeightSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.text8, 'String', {'Relative', 'Weight', sprintf('%2.0f%%', get(hObject,'Value'))});

% --- Executes during object creation, after setting all properties.
function RelWeightSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RelWeightSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
