function [N_procs,timer,trace] = load_timer_file(file)
if ~exist('file','var')
    file = 'Z:\xraid\output_test_cond2\test_plasma.1.timer';
end


% Read in the file
fid = fopen(file);
text = char(fread(fid,inf,'char')');
fclose(fid);
pause(0);

% Find the number of processors
i1 = strfind(text,'<N_procs=');
i2 = find(text=='>');
i2 = i2(find(i2>i1,1,'first'));
tmp = str2num(strrep(text(i1+9:i2-1),',id=',' ')); %#ok<ST2NM>
if isempty(tmp)
    error('N_procs was not defined');
end
N_procs = tmp(1);
if length(tmp) > 1
    id = tmp(2); %#ok<NASGU>
end

% Find the timers
i1 = strfind(text,'<timer:');
i2 = find(text=='>');
timer = struct();
for i = 1:length(i1)
    i3 = i2(find(i2>i1(i),1,'first'));
    tmp = text(i1(i):i3);
    j1 = strfind(tmp,'id=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    id = tmp(j1+3:j2-1);
    if ~isempty(str2num(id)) %#ok<ST2NM>
        timer(i).id = str2double(tmp(j1+3:j2-1));
    else
        timer(i).id = tmp(j1+3:j2-1);
    end
    j1 = strfind(tmp,'message=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    timer(i).message = tmp(j1+8:j2-1);
    j1 = strfind(tmp,'file=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    timer(i).file = tmp(j1+5:j2-1);
    j1 = strfind(tmp,'thread=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    timer(i).thread = str2double(tmp(j1+7:j2-1));
    j1 = strfind(tmp,'start=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    timer(i).start = str2double(tmp(j1+6:j2-1));
    j1 = strfind(tmp,'stop=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    timer(i).stop = str2double(tmp(j1+5:j2-1));
    j1 = strfind(tmp,'N=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    timer(i).N = str2double(tmp(j1+2:j2-1));
    j1 = strfind(tmp,'min=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    timer(i).min = str2double(tmp(j1+4:j2-1));
    j1 = strfind(tmp,'max=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    timer(i).max = str2double(tmp(j1+4:j2-1));
    j1 = strfind(tmp,'tot=');
    j2 = find(tmp==','|tmp=='>');
    j2 = j2(find(j2>j1,1,'first'));
    timer(i).tot = str2num(tmp(j1+4:j2-1)); %#ok<ST2NM>
end

% Find the traces
if exist(strrep(file,'.timer','.trace'),'file') && false
    % The detailed trace data exists, use it
    trace = load_trace_file(strrep(file,'.timer','.trace'));
    % We need to create the traces and the min, max ,and total times
    error('Unfinished');
else
    % The detailed trace data does not exist, use the data in the timer file
    i1 = strfind(text,'<trace:');
    i2 = find(text=='>');
    trace = struct;
    for i = 1:length(i1)
        i3 = i2(find(i2>i1(i),1,'first'));
        tmp = text(i1(i):i3);
        j1 = strfind(tmp,'id=');
        j2 = find(tmp==','|tmp=='>');
        j2 = j2(find(j2>j1,1,'first'));
        id = tmp(j1+3:j2-1);
        if ~isempty(str2num(id)) %#ok<ST2NM>
            trace(i).id = str2double(id);
        else
            trace(i).id = id;
        end
        j1 = strfind(tmp,'thread=');
        j2 = find(tmp==','|tmp=='>');
        j2 = j2(find(j2>j1,1,'first'));
        trace(i).thread = str2double(tmp(j1+7:j2-1));
        j1 = strfind(tmp,'N=');
        j2 = find(tmp==','|tmp=='>');
        j2 = j2(find(j2>j1,1,'first'));
        trace(i).N = str2double(tmp(j1+2:j2-1));
        j1 = strfind(tmp,'min=');
        j2 = find(tmp==','|tmp=='>');
        j2 = j2(find(j2>j1,1,'first'));
        trace(i).min = str2double(tmp(j1+4:j2-1));
        j1 = strfind(tmp,'max=');
        j2 = find(tmp==','|tmp=='>');
        j2 = j2(find(j2>j1,1,'first'));
        trace(i).max = str2double(tmp(j1+4:j2-1));
        j1 = strfind(tmp,'tot=');
        j2 = find(tmp==','|tmp=='>');
        j2 = j2(find(j2>j1,1,'first'));
        trace(i).tot = str2double(tmp(j1+4:j2-1));
        j1 = strfind(tmp,'active=[');
        j2 = find(tmp==']');
        j2 = j2(find(j2>j1,1,'first'));
        tmp2 = strtrim(tmp(j1+8:j2-1));
        if isempty(tmp2)
            trace(i).trace = str2num(tmp2); %#ok<ST2NM>
        else
            if isempty(tmp2)
                trace(i).trace = {};
            else
                k = [0 find(tmp2==' ') length(tmp2)+1];
                for j = 1:length(k)-1
                    trace(i).trace{1,j} = tmp2(k(j)+1:k(j+1)-1);
                end
            end
        end
    end
    % Find all timers that have no traces and create a single trace for each of them
    no_trace = [];
    id_trace = unique([trace.id]);
    id_timer = unique([timer.id]);
    for i = 1:length(id_timer)
        if isempty(find(id_trace==id_timer(i),1))
            no_trace = [no_trace,id_timer(i)]; %#ok<AGROW>
        end
    end
    for k = 1:length(no_trace)
        i = find([timer.id]==no_trace(k));
        j = length(trace)+1;
        trace(j).id = timer(i).id;
        trace(j).N = timer(i).N;
        trace(j).min = timer(i).min;
        trace(j).max = timer(i).max;
        trace(j).tot = timer(i).tot;
        trace(j).trace = [];
    end
end
1;



function trace = load_trace_file(file)
% Read the file
fid = fopen(file);
trace = struct([]);
while 1
    pos = ftell(fid);
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    i1 = strfind(tline,'id='); i1=i1(1);
    i2 = find(tline==',',1,'first');
    id = tline(i1+3:i2-1); %#ok<ST2NM>
    i1 = strfind(tline,'thread=');
    i2 = find(tline==','|tline=='>');
    i2 = i2(find(i2>i1,1,'first'));
    thread = str2double(tline(i1+7:i2-1));
    i1 = strfind(tline,'N='); i1=i1(1);
    i2 = find(tline==':',1,'first');
    N = str2num(tline(i1+2:i2-1)); %#ok<ST2NM>
    fseek(fid,pos+i2,'bof');
    tmp = fread(fid,2*N,'double');
    fread(fid,1,'char');
    i = length(trace)+1;
    if ~isempty(str2num(id)) %#ok<ST2NM>
        trace(i).id = str2num(id); %#ok<ST2NM>
    else
        trace(i).id = id;
    end
    trace(i).thread = thread;
    trace(i).start = tmp(1:N)';
    trace(i).stop = tmp(N+1:2*N)';
end
fclose(fid);
% % Combine the trace results
% id1 = {trace.id};
% id2 = unique(id1);
% trace2 = struct();
% thread1 = [trace.thread];
% thread2 = unique(thread1);
% k = 1;
% for i = 1:length(id2)
%     for j = 1:length(thread2)
%         index = find(cellfun(@(x) all(x==id2{i}),id1)&thread1==thread2(j));
%         if isempty(index)
%             continue;
%         end
%         trace2(k).id = id2{i};
%         trace2(k).thread = thread2(j);
%         trace2(k).start = [trace(index).start];
%         trace2(k).stop = [trace(index).stop];
%         k = k+1;
%     end
% end
% trace = trace2;
% 1;

