%
% modifying paths
%
if(not(isdeployed))

 mypaths =[{sprintf('%s%stransient%s'                           ,pwd ,filesep, filesep)}
           {sprintf('%s%slibrary%s'                             ,pwd ,filesep, filesep)}
           {sprintf('%s%slibrary%smatlab_general%s'             ,pwd ,filesep, filesep)}
           {sprintf('%s%slibrary%sid_simulation%s'              ,pwd ,filesep, filesep)}];

 for pathidx =1:length(mypaths)
   addpath(mypaths{pathidx});
 end

 clear mypaths pathidx
end
