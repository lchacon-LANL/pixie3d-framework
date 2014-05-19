function data = loadDebugData(filename)
% This function loads the debug information for the pixie3dApplication

fid = fopen(filename);

i = 1;
data = struct([]);
while 1
    % Read the header block
    tline = [];
    while isempty(tline)
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
    end
    if ~ischar(tline)
        break
    end
    iteration = str2num(tline(12:length(tline))); %#ok<ST2NM>
    if length(iteration)~=1
        error('Problem reading data');
    end
    tline = fgetl(fid);
    time  = str2num(tline(7:length(tline))); %#ok<ST2NM>
    tline = fgetl(fid);
    lower = str2num(tline(10:length(tline)-1)); %#ok<ST2NM>
    tline = fgetl(fid);
    upper = str2num(tline(10:length(tline)-1)); %#ok<ST2NM>
    tline = fgetl(fid);
    nbox  = str2num(tline(9:length(tline)-1)); %#ok<ST2NM>
    tline = fgetl(fid);
    N_levels = str2num(tline(11:length(tline))); %#ok<ST2NM>
    tline = fgetl(fid);
    N_vars = str2num(tline(9:length(tline))); %#ok<ST2NM>
    data(i).iteration = iteration;
    data(i).time = time;
    data(i).lower = lower;
    data(i).upper = upper;
    data(i).nbox = nbox;
    data(i).N_levels = N_levels;
    data(i).N_vars = N_vars;
    % Read each variable
    data(i).var = struct([]);
    for j = 1:data(i).N_vars
        tline = fgetl(fid);
        data(i).var(j).var_name = tline(12:length(tline));
        if data(i).N_levels > 0
            % Reading in debug type 1
            for k = 1:data(i).N_levels
                tline = fgetl(fid);
                index = find(tline=='=');
                level = str2num(tline(index(1)+1:index(2)-9))+1; %#ok<ST2NM>
                data(i).ratio{i} = str2num(tline(index(2)+3:index(3)-12)); %#ok<ST2NM>
                n_patch = str2num(tline(index(3)+1:length(tline))); %#ok<ST2NM>
                for m = 1:n_patch
                    t0 = ftell(fid);
                    tline = fgetl(fid);
                    fseek(fid,t0+length(tline)+1,'bof');
                    header = convert_patch_header(tline);
                    N = [header.ilast-header.ifirst+1+2*header.gcw, header.depth];
                    patch = header.patch_num+1;
                    if strcmp(header.type,'cell')
                        % Reading cell-centered data
                        data_read = fread(fid,prod(N),'double');
                        tline = fgetl(fid);
                        if ~isempty(tline)
                            fseek(fid,-length(tline),'cof');
                        end
                        data_read = reshape(data_read,N);
                        data(i).var(j).gcw = header.gcw;
                        data(i).var(j).depth = header.depth;
                        data(i).var(j).ifirst{level,patch} = header.ifirst;
                        data(i).var(j).ilast{level,patch} = header.ilast;
                        data(i).var(j).data{level,patch} = data_read;
                    elseif strfind(header.type,'side')==1
                        % Reading cell-centered data
                        data(i).var(j).gcw = header.gcw;
                        data(i).var(j).depth = header.depth;
                        data(i).var(j).ifirst{level,patch} = header.ifirst;
                        data(i).var(j).ilast{level,patch} = header.ilast;
                        for d = 1:3
                            N2 = N;
                            N2(d) = N2(d)+1;
                            data_read = fread(fid,prod(N2),'double');
                            tline = fgetl(fid);
                            if ~isempty(tline)
                                fseek(fid,-length(tline),'cof');
                            end
                            data_read = reshape(data_read,N2);
                            data(i).var(j).data{level,patch}{d} = data_read;
                        end
                    else
                        error('Unknown data format');
                    end
                end
            end
        else
            % Reading in debug type 2
            tline = fgetl(fid);
            index = find(tline=='=');
            data(i).var(j).gcw = [0 0 0];
            depth = str2num(tline(index+1:length(tline))); %#ok<ST2NM>
            N = prod(data(i).nbox)*depth;
            data_read = fread(fid,N,'double');
            data_read = reshape(data_read,[data(i).nbox,depth]);
            data(i).var(j).ifirst{1} = [0 0 0];
            data(i).var(j).ilast{1} = data(i).nbox-1;
            data(i).var(j).depth = depth;
            data(i).var(j).data{1} = data_read;
            tline = fgetl(fid); %#ok<NASGU>
            if ~isempty(tline)
                fseek(fid,-length(tline),'cof');
            end
        end
    end
    if data(i).N_levels == -1
        data(i).N_levels = 1;
    end
    i = i+1;
end
fclose(fid);


function header = convert_patch_header(tline)
% Get the patch number
index = strfind(tline,'patch_num');
tline2 = tline(index:length(tline));
tline2 = tline2(find(tline=='=',1)+1:length(tline2));
header.patch_num = str2double(tline2(1:find(tline2==',',1)-1));
% Get ifirst
index = strfind(tline,'ifirst');
tline2 = tline(index:length(tline));
tline2 = tline2(find(tline2=='(',1)+1:find(tline2==')',1)-1);
header.ifirst = str2num(tline2); %#ok<ST2NM>
% Get ilast
index = strfind(tline,'ilast');
tline2 = tline(index:length(tline));
tline2 = tline2(find(tline2=='(',1)+1:find(tline2==')',1)-1);
header.ilast = str2num(tline2); %#ok<ST2NM>
% Get gcw
index = strfind(tline,'gcw');
tline2 = tline(index:length(tline));
tline2 = tline2(find(tline2=='(',1)+1:find(tline2==')',1)-1);
header.gcw = str2num(tline2); %#ok<ST2NM>
% Get depth
index = strfind(tline,'depth');
tline2 = tline(index:length(tline));
tline2 = tline2(find(tline2=='=',1)+1:length(tline2));
header.depth = str2double(tline2(1:min([find(tline2==',',1)-1,length(tline2)])));
% Get type
tline2 = tline(find(tline==',',1,'last')+1:length(tline));
header.type = strtrim(tline2(1:min([find(tline2==',',1)-1,length(tline2)])));

