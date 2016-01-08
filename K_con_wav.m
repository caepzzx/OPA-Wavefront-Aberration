function [ K_conwav] = K_con_wav( K_con )
%K_CON_wav 给出参量过程中三波的耦合项的参数在指定波长时的取值[zzx]
%   K_con为元胞数组 
%   K_con{1}为S光中耦合项参数随波长的变化
%   K_con{2}为I光中耦合项参数随波长的变化
%   K_con{3}为P光中耦合项参数随波长的变化
%   size(K_con(j),2)==1 OR numel(t)
for j=1:numel(K_con)
    if size(K_con{j},2)==1
        K_conwav(j)=K_con{j};
    else
        temp_K_con=K_con{j};
        idx=numel(temp_K_con)/2;
        if(idx~=0)
            idx=ceil(idx);
        K_conwav(j)=temp_K_con(idx);
        else
        K_conwav(j)=mean([temp_K_con(idx),temp_K_con(idx+1)]);
        end
    end
end

