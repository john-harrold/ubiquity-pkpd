function []=build_exe(exe_name)
%
% Creates a stand alone executable from model_gui.m in the directory above
% this one. To specify the name of the executable call this function with a
% string. 
%
% Examples:
%  Create model_gui.exe in directory above:
%  >> build_exe 
% 
%  Create my_model.exe in directory above:
%  >> build_exe('my_model')
% 
%  Create my_model-yyyy_mm_dd.exe with the current date in place of the
%  placeholders in directory above:
%  >> build_exe(sprintf('my_model-%s', datestr(now, 'yyyy_mm_dd')))
%


try
    delete pdm_gui_log.txt;
end

try 
    delete gui_state.mat;
end

mcc -e model_gui.m -a . -d ..;

if((exist('exe_name', 'var')))
  movefile(sprintf('..%smodel_gui.exe',filesep), sprintf('..%s%s.exe',filesep, exe_name));
end
