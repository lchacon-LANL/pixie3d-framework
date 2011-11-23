% This is a MATLAB program that will check if all of the ghost cells in a 
% simulation match any other cells on the same level

file = '/projects/pixie3d-v3/bin/debugFile4';
data = loadDebugData(file);

error = zeros(length(data),data(1).N_vars);
for i = 1:length(data)
    for j = 1:data(i).N_vars
        level_data = nan([data(i).nbox+2*data(i).var(j).gcw,data(i).var(j).depth]);
        for k = 1:size(data(i).var(j).data,2)
            if any(isnan(data(i).var(j).data{1,k}(:)))
                error('This program relies on the data having no NaNs');
            end
            i1 = data(i).var(j).ifirst{1,k}(1)+1:data(i).var(j).ilast{1,k}(1)+1+2*data(i).var(j).gcw(1);
            i2 = data(i).var(j).ifirst{1,k}(2)+1:data(i).var(j).ilast{1,k}(2)+1+2*data(i).var(j).gcw(2);
            i3 = data(i).var(j).ifirst{1,k}(3)+1:data(i).var(j).ilast{1,k}(3)+1+2*data(i).var(j).gcw(3);
            patch_data = data(i).var(j).data{1,k};
            tmp = level_data(i1,i2,i3,:);
            diff = tmp-patch_data;
            diff(isnan(tmp)) = 0;
            patch_error = max(abs(diff(:)));
            if patch_error>0
                imagesc(diff(:,:,2,1))
                1;
            end
            error(i,j) = max(error(i,j),patch_error);
            level_data(i1,i2,i3,:) = patch_data;
        end
    end
end
if max(max(error))==0
    fprintf(1,'No errors detected\n');
else
    fprintf(1,'Errors detected\n');
end
