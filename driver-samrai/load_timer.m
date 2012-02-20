function varargout = load_timer(varargin)
% LOAD_TIMER M-file for load_timer.fig
%      LOAD_TIMER, by itself, creates a new LOAD_TIMER or raises the existing
%      singleton*.
%
%      H = LOAD_TIMER returns the handle to a new LOAD_TIMER or the handle to
%      the existing singleton*.
%
%      LOAD_TIMER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOAD_TIMER.M with the given input arguments.
%
%      LOAD_TIMER('Property','Value',...) creates a new LOAD_TIMER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before load_timer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to load_timer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help load_timer

% Last Modified by GUIDE v2.5 21-Sep-2010 13:19:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @load_timer_OpeningFcn, ...
                   'gui_OutputFcn',  @load_timer_OutputFcn, ...
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


% --- Executes just before load_timer is made visible.
function load_timer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to load_timer (see VARARGIN)
handles.output = hObject;
set(handles.load_plot,'Position',[0.05 0.05 0.9 0.18])
set(handles.load_text,'Position',[0.05 0.230 0.3 0.025])
set(handles.timer_table,'Position',[0.05 0.26 0.9 0.67])
set(handles.function_text,'Position',[0.05 0.93 0.8 0.025])
set(handles.load,'Position',[0.07 0.962 0.07 0.028])
set(handles.reset,'Position',[0.16 0.962 0.07 0.028])
set(handles.select_proc_text,'Position',[0.38 0.955 0.07 0.028],'HorizontalAlignment','right')
set(handles.select_proc,'Position',[0.47 0.96 0.08 0.028])
set(handles.select_proc,'Visible','off');
set(handles.select_proc_text,'Visible','off');
set(handles.load_plot,'Visible','off');
set(handles.load_text,'Visible','off');
set(handles.timer_table,'Position',[0.05 0.05 0.9 0.88])
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = load_timer_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function display_data(handles)
% This function displays all of the data to the figure window
% Load the data for the function selected
N_proc = handles.data(1).N_procs;
if isempty(handles.call)
    % Get the relavent timers for each processor
    for p = 1:N_proc
        timers = handles.data(p).timer;
        set(handles.function_text,'String','All Functions');
        timers_proc{p} = timers; %#ok<AGROW>
    end
else
    % Get the relavent timers for each processor
    for p = 1:N_proc
        timers = handles.data(p).timer;
        trace = handles.data(p).trace;
        % Find the timers in the selected call hierarcy
        index = [];
        for i = 1:size(handles.call,1)
            j = find(strcmp({timers.message},handles.call{i,1})&strcmp({timers.file},handles.call{i,2}));
            if isempty(j)
                index = -1;
                break;
            end
            index = [index,timers(j).id]; %#ok<AGROW>
        end
        if index(1) == -1
            timers_proc{p} = []; %#ok<AGROW>
        end
        index = unique(index);
        % Keep only those traces that have all of the selected call hierarchy
        for i = length(trace):-1:1
            for j = 1:length(index)
                if trace(i).id~=index(j) && all(trace(i).trace~=index(j))
                    trace(i) = [];
                    break;
                end
            end
        end
        % Keep only the timers that have traces remaining
        for i = length(timers):-1:1
            if all(timers(i).id~=[trace.id])
                timers(i) = [];
            end
        end
        % Add the information from the trace logs
        for j = 1:length(timers)
            i = find([trace.id]==timers(j).id&[trace.thread]==timers(j).thread);
            timers(j).N   = sum([trace(i).N]);
            timers(j).min = min([trace(i).min]);
            timers(j).max = max([trace(i).max]);
            timers(j).tot = sum([trace(i).tot]);
        end
        index = find([timers.tot]==0);
        timers(index) = []; %#ok<FNDSB,NASGU>
        timers_proc{p} = timers; %#ok<AGROW>
    end
end
% Create the timers to display in the table
if get(handles.select_proc,'Value')==1
    % We want to average each timer
    % First, get the timers for each processor, averaging the threads
    timers = timers_proc{p}(1); % Initialize the structure
    k = 1;
    for p = 1:N_proc
        id = unique([timers_proc{p}.id]);
        for i = 1:length(id)
            j = find([timers_proc{p}.id]==id(i));
            timers(k) = timers_proc{p}(j(1));
            timers(k).thread = [timers_proc{p}(j).thread];
            timers(k).N = mean([timers_proc{p}(j).N]);
            timers(k).min = mean([timers_proc{p}(j).min]);
            timers(k).max = mean([timers_proc{p}(j).max]);
            timers(k).tot = mean([timers_proc{p}(j).tot]);
            k = k+1;
        end
    end
    i = 1;
    while i <= length(timers)
        j = find([timers.id]==timers(i).id);
        timers(i).N = round(sum([timers(j).N])/N_proc); %#ok<AGROW>
        timers(i).min = sum([timers(j).min])/N_proc; %#ok<AGROW>
        timers(i).max = sum([timers(j).max])/N_proc; %#ok<AGROW>
        timers(i).tot = sum([timers(j).tot])/N_proc; %#ok<AGROW>
        timers(j(2:length(j))) = []; %#ok<AGROW>
        i = i+1;
    end
elseif get(handles.select_proc,'Value')==2
    % We want to take the minimum value for each processor
    timers = timers_proc{p}(1); % Initialize the structure
    k = 1;
    for p = 1:N_proc
        id = unique([timers_proc{p}.id]);
        for i = 1:length(id)
            j = find([timers_proc{p}.id]==id(i));
            timers(k) = timers_proc{p}(j(1));
            timers(k).thread = [timers_proc{p}(j).thread];
            timers(k).N = min([timers_proc{p}(j).N]);
            timers(k).min = min([timers_proc{p}(j).min]);
            timers(k).max = min([timers_proc{p}(j).max]);
            timers(k).tot = min([timers_proc{p}(j).tot]);
            k = k+1;
        end
    end
    i = 1;
    while i <= length(timers)
        j = find([timers.id]==timers(i).id);
        timers(i).N = min([timers(j).N]); %#ok<AGROW>
        timers(i).min = min([timers(j).min]); %#ok<AGROW>
        timers(i).max = min([timers(j).max]); %#ok<AGROW>
        timers(i).tot = min([timers(j).tot]); %#ok<AGROW>
        timers(j(2:length(j))) = []; %#ok<AGROW>
        i = i+1;
    end
elseif get(handles.select_proc,'Value')==3
    % We want to take the maximum value for each processor
    timers = timers_proc{p}(1); % Initialize the structure
    k = 1;
    for p = 1:N_proc
        id = unique([timers_proc{p}.id]);
        for i = 1:length(id)
            j = find([timers_proc{p}.id]==id(i));
            timers(k) = timers_proc{p}(j(1));
            timers(k).thread = [timers_proc{p}(j).thread];
            timers(k).N = mean([timers_proc{p}(j).N]);
            timers(k).min = max([timers_proc{p}(j).min]);
            timers(k).max = max([timers_proc{p}(j).max]);
            timers(k).tot = max([timers_proc{p}(j).tot]);
            k = k+1;
        end
    end
    i = 1;
    while i <= length(timers)
        j = find([timers.id]==timers(i).id);
        timers(i).N = max([timers(j).N]); %#ok<AGROW>
        timers(i).min = max([timers(j).min]); %#ok<AGROW>
        timers(i).max = max([timers(j).max]); %#ok<AGROW>
        timers(i).tot = max([timers(j).tot]); %#ok<AGROW>
        timers(j(2:length(j))) = []; %#ok<AGROW>
        i = i+1;
    end
else
    % We have selected a specific processor
    timers = timers_proc{get(handles.select_proc,'Value')-3};
end
tot_time = max([timers.tot]);
% Update the table
if isempty(handles.call)
    set(handles.function_text,'String','All Functions');
else
    tmp = handles.call(:,1);
    i = 1;
    while i < length(tmp)
        j = find(strcmp(tmp{i},tmp));
        j(1) = [];
        if ~isempty(j)
            tmp(j) = [];
        end
        i = i+1;
    end
    text = tmp{1};
    for i = 2:length(tmp)
        text = sprintf('%s  ->  %s',text,tmp{i});
    end
    set(handles.function_text,'String',text);
end
table_data = cell(length(handles.data(1).timer),10);
for i = 1:length(timers)
    table_data{i,1} = timers(i).message;
    table_data{i,2} = timers(i).file;
    if length(timers(i).thread)==1
        table_data{i,3} = timers(i).thread;
    else
        table_data{i,3} = 'multiple';
    end
    table_data{i,4} = timers(i).start;
    table_data{i,5} = timers(i).stop;
    table_data{i,6} = timers(i).N;
    table_data{i,7} = timers(i).min;
    table_data{i,8} = timers(i).max;
    table_data{i,9} = timers(i).tot;
    table_data{i,10} = 100*timers(i).tot/tot_time;
end
% Plot the data in the table
[Y,I] = sort([table_data{:,9}],2,'descend');
table_data = table_data(I,:);
set(handles.timer_table,'data',table_data);
% Plot the load balance
if N_proc > 1 
    [v,i] = max([timers.tot]);
    id = timers(i).id;
    total_time = zeros(1,N_proc);
    for i = 1:N_proc
        j = find([timers_proc{i}.id]==id);
        total_time(1:length(j),i) = [timers_proc{i}(j).tot];
    end
    axes(handles.load_plot);
    if size(total_time,2)*sqrt(size(total_time,1)) <= 128
        bar(total_time')
    elseif size(total_time,1)==1
        plot(total_time)
    else
        error('Not finished');
    end
    axis([0 N_proc+1 0 1.05*max(total_time(:))+1])
    hold on
    plot([0 N_proc+1],[mean(total_time(:)) mean(total_time(:))],'r--')
    hold off
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles) %#ok<INUSL,INUSD,DEFNU>
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'pathname')
    [FileName,PathName] = uigetfile('*.timer','Select the timer file',handles.pathname);
else
    [FileName,PathName] = uigetfile('*.timer','Select the timer file');
end
if FileName == 0
    return;
end
n = length(FileName);
if ~strcmp(FileName(max(n-5,1):n),'.timer')
    fprintf(1,'Not a valid filename\n');
    return;
end
[data.N_procs,data.timer,data.trace] = load_timer_file([PathName,FileName]);
if data(1).N_procs > 1 
    % Load the data for each timer file
    i = find(FileName=='.',2,'last');
    for j = 1:data(1).N_procs
        file = [FileName(1:i(1)),num2str(j),'.timer'];
        [data(j).N_procs,data(j).timer,data(j).trace] = load_timer_file([PathName,file]); %#ok<AGROW>
    end
    if ~all([data.N_procs]==data(1).N_procs)
        error('# of processors must match for all files');
    end
    % Change the ids of the timers
    if isnumeric(data(1).timer(1).id)
        % This is the old format, we need to re-number so the processors
        % share the same ids
        timers = [];
        for p = 1:data(1).N_procs
            timers = [timers,data(p).timer]; %#ok<AGROW>
        end
        i = 1;
        while i <= length(timers)
            timers(i).id = i-1; %#ok<AGROW>
            for j = length(timers):-1:i+1
                if ( strcmp(timers(j).message,timers(i).message) && strcmp(timers(j).file,timers(i).file) ...
                        && timers(j).start==timers(i).start && timers(j).stop==timers(i).stop )
                    timers(j) = []; %#ok<AGROW>
                end
            end
            i = i+1;
        end
        for p = 1:data(1).N_procs
            map = zeros(1,length(data(p).timer));
            for j = 1:length(data(p).timer)
                k = find(strcmp(data(p).timer(j).message,{timers.message}) & ...
                    strcmp(data(p).timer(j).file,{timers.file}) & ...
                    data(p).timer(j).start==[timers.start] & ...
                    data(p).timer(j).stop==[timers.stop] );
                map(j) = k-1;
                data(p).timer(j).id = k-1; %#ok<AGROW>
            end
            for j = 1:length(data(p).trace)
                data(p).trace(j).id = map(data(p).trace(j).id+1); %#ok<AGROW>
                if ~isempty(data(p).trace(j).trace)
                    data(p).trace(j).trace = map(data(p).trace(j).trace+1); %#ok<AGROW>
                end
            end
        end
    else
        % This is the new format, all processors use the same id, but it is
        % an alpha-numeric id.
        map = {};
        for i = 1:length(data)
            map = [map {data(i).timer.id}]; %#ok<AGROW>
        end
        map = unique(map);
        for i = 1:length(data)
            id1 = {data(i).timer.id};
            id2 = {data(i).trace.id};
            for j = 1:length(map);
                index = find(cellfun(@(x) all(x==map{j}),id1));
                for k = 1:length(index)
                    data(i).timer(index(k)).id = j; %#ok<AGROW>
                end
                index = find(cellfun(@(x) all(x==map{j}),id2));
                for k = 1:length(index)
                    data(i).trace(index(k)).id = j; %#ok<AGROW>
                    trace_id1 = data(i).trace(index(k)).trace;
                    trace_id2 = zeros(size(trace_id1));
                    for m = 1:length(trace_id1)
                        trace_id2(m) = find(cellfun(@(x) all(x==trace_id1{m}),map));
                    end
                    data(i).trace(index(k)).trace = trace_id2;
                end
            end
        end
    end
    % Set the timer selection button
    string_text{1,1} = 'Average';
    string_text{2,1} = 'Minimum';
    string_text{3,1} = 'Maximum';
    for i = 1:data(1).N_procs
        string_text{i+3,1} = ['Proc ',num2str(i)];
    end
    set(handles.select_proc,'String',string_text);
    set(handles.select_proc,'Visible','on');
    set(handles.select_proc_text,'Visible','on');
    set(handles.load_plot,'Visible','on');
    set(handles.load_text,'Visible','on');
    set(handles.load_plot,'Position',[0.05 0.05 0.9 0.18])
    set(handles.load_text,'Position',[0.05 0.230 0.3 0.025])
    set(handles.timer_table,'Position',[0.05 0.26 0.9 0.67])
else
    data.N_procs = 1;
    set(handles.select_proc,'Visible','off');
    set(handles.select_proc_text,'Visible','off');
    set(handles.load_plot,'Visible','off');
    set(handles.load_text,'Visible','off');
    set(handles.timer_table,'Position',[0.05 0.05 0.9 0.88])
end
set(handles.function_text,'Position',[0.05 0.93 0.8 0.025])
set(handles.load,'Position',[0.07 0.962 0.07 0.028])
set(handles.reset,'Position',[0.16 0.962 0.07 0.028])
set(handles.select_proc_text,'Position',[0.38 0.955 0.07 0.028],'HorizontalAlignment','right')
set(handles.select_proc,'Position',[0.47 0.96 0.08 0.028])
handles.data = data;
handles.call = [];
handles.pathname = PathName;
guidata(hObject,handles);
display_data(handles);


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.call = [];
guidata(hObject,handles);
display_data(handles);


% --- Executes when selected cell(s) is changed in timer_table.
function timer_table_CellSelectionCallback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to timer_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if isempty(eventdata.Indices)
    % Nothing is selected
elseif eventdata.Indices(2) == 1 
    % We want to dive into the function
    data = get(handles.timer_table,'Data');
    call{1} = data{eventdata.Indices(1),1};
    call{2} = data{eventdata.Indices(1),2};
    handles.call = [handles.call; call];
    guidata(hObject, handles);
    display_data(handles)
else
    % No options for this case
    % set(handles.FunctionTable
end


% --- Executes on selection change in select_proc.
function select_proc_Callback(hObject, eventdata, handles) %#ok<INUSL,INUSD,DEFNU>
% hObject    handle to select_proc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display_data(handles)



% --- Executes during object creation, after setting all properties.
function select_proc_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to select_proc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


