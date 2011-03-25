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
                level = str2num(tline(index(1)+1:index(2)-9)); %#ok<ST2NM>
                data(i).ratio{i} = str2num(tline(index(2)+3:index(3)-12)); %#ok<ST2NM>
                n_patch = str2num(tline(index(3)+1:length(tline))); %#ok<ST2NM>
                for m = 1:n_patch
                    t0 = ftell(fid);
                    tline = fgetl(fid);
                    fseek(fid,t0+length(tline)+1,'bof');
                    index = find(tline=='=');
                    patch_num = str2num(tline(index(1)+1:index(2)-10)); %#ok<ST2NM>
                    ifirst = str2num(tline(index(2)+3:index(3)-10)); %#ok<ST2NM>
                    ilast = str2num(tline(index(3)+3:index(4)-8)); %#ok<ST2NM>
                    gcw = str2num(tline(index(4)+3:index(5)-10)); %#ok<ST2NM>
                    depth = str2num(tline(index(5)+1:length(tline))); %#ok<ST2NM>
                    N = prod(ilast-ifirst+1+2*gcw)*depth;
                    data_read = fread(fid,N,'double');
                    tline = fgetl(fid);
                    if ~isempty(tline)
                        fseek(fid,-length(tline),'cof');
                    end
                    data_read = reshape(data_read,[ilast-ifirst+1+2*gcw,depth]);
                    data(i).var(j).gcw = gcw;
                    data(i).var(j).depth = depth;
                    data(i).var(j).ifirst{level+1,patch_num+1} = ifirst;
                    data(i).var(j).ilast{level+1,patch_num+1} = ilast;
                    data(i).var(j).data{level+1,patch_num+1} = data_read;
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
            data(i).var(j).ifirst{1} = [1 1 1];
            data(i).var(j).ilast{1} = data(i).nbox;
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

