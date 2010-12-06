% This is a MATLAB program that will compare two debug files to search for
% differences.  The absolute error is recorded in "error".  It was created
% to compare a single level single patch run with a single level multiple
% patch run.

file1 = '/projects/pixie3d/bin/debugFile1';
file2 = '/projects/pixie3d/bin/debugFile9';

data1 = loadDebugData(file1);
data2 = loadDebugData(file2);

error = zeros(length(data1),data1(1).N_vars);
for i = 1:length(data1)
    for j = 1:data1(i).N_vars
        d1 = zeros([data1(i).nbox+data1(i).var(j).gcw,data1(i).var(j).depth]);
        for k = 1:size(data1(i).var(j).data,2)
            i1 = data1(i).var(j).ifirst{1,k}(1)+1:data1(i).var(j).ilast{1,k}(1)+1+2*data1(i).var(j).gcw(1);
            i2 = data1(i).var(j).ifirst{1,k}(2)+1:data1(i).var(j).ilast{1,k}(2)+1+2*data1(i).var(j).gcw(2);
            i3 = data1(i).var(j).ifirst{1,k}(3)+1:data1(i).var(j).ilast{1,k}(3)+1+2*data1(i).var(j).gcw(3);
            d1(i1,i2,i3,:) = data1(i).var(j).data{1,k};
        end
        d2 = zeros([data2(i).nbox+data2(i).var(j).gcw,data2(i).var(j).depth]);
        for k = 1:size(data2(i).var(j).data,2)
            i1 = data2(i).var(j).ifirst{1,k}(1)+1:data2(i).var(j).ilast{1,k}(1)+1+2*data2(i).var(j).gcw(1);
            i2 = data2(i).var(j).ifirst{1,k}(2)+1:data2(i).var(j).ilast{1,k}(2)+1+2*data2(i).var(j).gcw(2);
            i3 = data2(i).var(j).ifirst{1,k}(3)+1:data2(i).var(j).ilast{1,k}(3)+1+2*data2(i).var(j).gcw(3);
            d2(i1,i2,i3,:) = data2(i).var(j).data{1,k};
        end
        diff = d1-d2;
        error(i,j) = max(abs(diff(:)));
    end
end
if max(max(error))==0
    fprintf(1,'No errors detected\n');
else
    fprintf(1,'Errors detected\n');
end


% i = 5;
% j = 23;
% d1 = zeros([data1(i).nbox+data1(i).var(j).gcw,data1(i).var(j).depth]);
% for k = 1:size(data1(i).var(j).data,2)
%     i1 = data1(i).var(j).ifirst{1,k}(1)+1:data1(i).var(j).ilast{1,k}(1)+1+2*data1(i).var(j).gcw(1);
%     i2 = data1(i).var(j).ifirst{1,k}(2)+1:data1(i).var(j).ilast{1,k}(2)+1+2*data1(i).var(j).gcw(2);
%     i3 = data1(i).var(j).ifirst{1,k}(3)+1:data1(i).var(j).ilast{1,k}(3)+1+2*data1(i).var(j).gcw(3);
%     d1(i1,i2,i3,:) = data1(i).var(j).data{1,k};
% end
% d2 = zeros([data2(i).nbox+data2(i).var(j).gcw,data2(i).var(j).depth]);
% for k = 1:size(data2(i).var(j).data,2)
%     i1 = data2(i).var(j).ifirst{1,k}(1)+1:data2(i).var(j).ilast{1,k}(1)+1+2*data2(i).var(j).gcw(1);
%     i2 = data2(i).var(j).ifirst{1,k}(2)+1:data2(i).var(j).ilast{1,k}(2)+1+2*data2(i).var(j).gcw(2);
%     i3 = data2(i).var(j).ifirst{1,k}(3)+1:data2(i).var(j).ilast{1,k}(3)+1+2*data2(i).var(j).gcw(3);
%     d2(i1,i2,i3,:) = data2(i).var(j).data{1,k};
% end
% diff = d1-d2;
% imagesc(diff(:,:,2,2)'),colorbar

