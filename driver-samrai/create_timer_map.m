function map = create_timer_map(trace)
% This function creates a timer map

% Return if there is not work
if isempty(trace)
    map = [];
    return;
end

% Find all head traces (those with no sum traces)
index = [];
for i = 1:length(trace)
    if isempty(trace(i).trace)
        index = [index,i]; %#ok<AGROW>
    end
end

% Create the call map
map = struct();
for i = 1:length(index)
    map(i).id = trace(index(i)).id;
    map(i).N = trace(index(i)).N;
    map(i).min = trace(index(i)).min;
    map(i).max = trace(index(i)).max;
    map(i).tot = trace(index(i)).tot;
    map(i).children = [];
end
trace(index) = [];


% Find (and recursivly add) any traces that are only the children of one head
for i = 1:length(map)
    index = [];
    for j = 1:length(trace)
        if any(trace(j).trace==map(i).id)
            add_trace = true;
            for k = 1:length(map)
                if k==i
                    continue;
                end
                if any(trace(j).trace==map(k).id)
                    add_trace = false;
                    continue;
                end
            end
            if add_trace
                index = [index,j]; %#ok<AGROW>
            end
        end
    end
    trace2 = trace(index);
    for j = 1:length(trace2)
        trace2(j).trace(trace2(j).trace==map(i).id) = [];
    end
    trace(index) = [];
    map(i).children = create_timer_map(trace2);
end
if isempty(trace)
    return;
end
check_map(map);


% Process the cases where a trace is a child of multiple parents, but no other processes
index = [];
for i = 1:length(trace)
    N_parents = 0;
    for j = 1:length(map)
        if any(trace(i).trace==map(j).id)
            N_parents = N_parents+1;
        end
    end
    if length(trace(i).trace) == N_parents
        map = add_complex_trace(map,trace(i));
        index = [index,i]; %#ok<AGROW>
    end
end
trace(index) = [];

% Process all the remaining cases where a trace is a child of more than one parent
for i = 1:length(trace)
    map = add_complex_trace(map,trace(i));
end

% Check the map
check_map(map);



function map = add_complex_trace(map,trace)
% This function adds a single trace that is a potential child of multiple parents
if length(trace) ~= 1
    error('Only 1 trace is suported');
end
% Get a list of all potential paths to store the child
paths = {};
for i = trace.trace
    % Find the possible paths for teh given parent timer
    tmp = search_map(map,i);
    % Check if it is a valid path
    for j = 1:length(tmp)
        keep_path = true;
        for k = 1:length(tmp{j})
            if isempty(find(trace.trace==tmp{j}(k),1))
                keep_path = false;
            end
        end
        if keep_path
            paths{length(paths)+1,1} = tmp{j}; %#ok<AGROW>
        end
    end
end
% Choose the appropriate path(s) to save the trace
if isempty(paths) 
    error('No where to store the path');
end
% Choose the appropriate paths to save the data
if length(paths)>1
    % Eliminate any paths that are sub-paths of any other path
    i = 1;
    while i <= length(paths)
        index = [];
        for j = 1:length(paths)
            if i == j
                continue;
            elseif length(paths{j})>length(paths{i})
                continue;
            elseif all(paths{j}==paths{i}(1:length(paths{j})))
                index = [index,j]; %#ok<AGROW>
            end
        end
        paths(index) = []; %#ok<AGROW>
        if isempty(index)
            i = i+1;
        else
            i = 1;
        end
    end
    1;
end
% Save the trace in the given paths
if length(paths) == 1
    if length(paths{1}) ~= length(trace.trace)
        error('There are more sub traces to account for');
    end
    map = recursively_add(map,trace,paths{1});
else
    if any(unique([paths{:}])~=sort(trace.trace))
        error('There are more sub traces to account for');
    end
    for i = 1:length(paths)
        map = recursively_add(map,trace,paths{i});
    end
end
1;   


function map = recursively_add(map,trace,path)
% This function will recursively add the trace entry
i = find([map.id]==path(1));
if length(path) > 1
    map(i).children = recursively_add(map(i).children,trace,path(2:length(path)));
    return;
end
if ~isempty(map(i).children)
    j = find([map(i).children.id]==trace.id);
    if ~isempty(j)
        % An entry for the trace already exists, add the two together
        map(i).children(j).id  = map(i).children(j).id;
        map(i).children(j).N   = map(i).children(j).N   + trace.N;
        map(i).children(j).min = map(i).children(j).min + trace.min;
        map(i).children(j).max = map(i).children(j).max + trace.max;
        map(i).children(j).tot = map(i).children(j).tot + trace.tot;
        return;
    end
end
j = length(map(i).children)+1;
map(i).children(j).id = trace.id;
map(i).children(j).N  = trace.N;
map(i).children(j).min = trace.min;
map(i).children(j).max = trace.max;
map(i).children(j).tot = trace.tot;
map(i).children(j).children = [];
1;


function path = search_map(map,id)
% This function searches a map for a given id
% If the id is found arrays of the paths to the id are returned
path = {};
for i = 1:length(map)
    if map(i).id == id
        path{length(path)+1,1} = id; %#ok<AGROW>
        continue;
    end
    if isempty(map(i).children)
        continue;
    end
    tmp = search_map(map(i).children,id);
    if isempty(tmp)
        continue;
    end
    for j = 1:length(tmp)
        path{length(path)+1,1} = [map(i).id tmp{j}]; %#ok<AGROW>
    end
end


function check_map(map)
if isempty(map)
    return;
end
for i = 1:length(map)
    for j = 1:length(map(i).children)
        if sum([map(i).children.id]==map(i).children(j).id)>1
            error('Duplicate children detected');
        end
    end
    check_map(map(i).children);
end
