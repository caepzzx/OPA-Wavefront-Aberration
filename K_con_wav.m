function [ K_conwav] = K_con_wav( K_con,j)
%K_CON_wav 给出参量过程中三波的耦合项的参数在指定波长时的取值[zzx]
%   K_con为元胞数组 
%   K_con{1}为S光中耦合项参数随波长的变化
%   K_con{2}为I光中耦合项参数随波长的变化
%   K_con{3}为P光中耦合项参数随波长的变化
%   i 耦合项的索引
%   j 波长的索引
%   size(K_con(j),2)==1 OR numel(t)

 
for i=1:numel(K_con)
    if size(K_con{i},2)==1
        K_conwav(i)=K_con{i};
    else
        temp_K_con=K_con{i};
        K_conwav(i)=temp_K_con(j);
        end
    end
end

