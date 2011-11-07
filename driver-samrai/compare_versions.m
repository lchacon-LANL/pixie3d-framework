% This is a MATLAB program that will compare two debug files to search for
% differences.  The absolute error is recorded in "error".  It was created
% to compare a single level single patch run with a single level multiple
% patch run.

file1 = '/projects/pixie3d/bin/debugFile';
file2 = '/projects/tmp/pixie3d/bin/debugFile';
% file2 = '/projects/pixie3d/bin/debugFile.good';


data1 = loadDebugData(file1);
data2 = loadDebugData(file2);

var_name1 = {data1(1).var.var_name}';
var_name2 = {data2(1).var.var_name}';

found_diff = false;
h1 = figure;
%h2 = figure;
diff_max = zeros(length(data1),1);
for n = 1:length(data1)
    for i2 = 1:length(var_name2)
        i1 = find(strcmp(var_name1,var_name2{i2}));
        if isempty(i1)
            continue;
        end
        tmp1 = data1(n).var(i1).data{1};
        tmp2 = data2(n).var(i2).data{1};
        % diff = (tmp1-tmp2)./max(abs(tmp1(:)));
        diff = tmp1-tmp2;
        diff(isnan(diff)) = 0;
        diff_max(n) = max(diff_max(n),max(abs(diff(:))));
        if max(abs(diff(:))) == 0
            continue;
        elseif max(abs(diff(:))) < 1e-6
            fprintf(1,'Small difference detected in %s (%e)\n',var_name1{i1},max(abs(diff(:))));
            continue;
        end
        figure(h1), imagesc(diff(:,:,2,1)), colorbar
        %figure(h2), imagesc(tmp2(:,:,2,1)), colorbar
        fprintf(1,'Difference detected in %s  (%e)\n',var_name1{i1},max(abs(diff(:))));
        found_diff = true;
        % pause
    end
    if found_diff
        fprintf(1,'Diferenced detected in iteration %i\n',n);
        break;
    end
end
close(h1)
% close(h2)
% fprintf('No differences detected\n');
