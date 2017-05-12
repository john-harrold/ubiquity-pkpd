function []=build_system(usersysfile)

sysfile = 'system.txt';

if(exist('usersysfile', 'var'))
   sysfile=usersysfile;
end

% generating the system
disp('-----------------------------------')
disp('Building the system to    ')
disp('generate the model targets');
disp('   -> Matlab: C/Simulink   ');
disp('   -> Matlab: m-file       ');
perl_output = perl('build_system.pl', sysfile);

if(not(strcmp(perl_output, '')))
    disp('Build reported errors and');
    disp('may have failed, see below:');
    disp(perl_output);
end

% If there is a simulink license, we'll compile it
myversion = ver;
disp('-----------------------------------')
if(any(strcmp('Simulink', {myversion.Name})))
  disp('Simulink found');
  disp('mex-ing the model file');
  try
  mex ode_model.c;
  catch
  disp('mex failed, you will only be able');
  disp('to use the m-file format');
  disp('To debug try the following');
  disp('mex ode_model.c');
  end
else
  disp('Simulink _NOT_ found');
  disp('C-matlab targets will not work');
end
disp('-----------------------------------')

clear myversion, perl_output;
