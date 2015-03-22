
% generating the system
disp('-----------------------------------')
disp('Running the perl script to')
disp('generate the model targets');
disp('   -> Adapt 5');
disp('   -> Berkeley Madonna');
disp('   -> Matlab: C/Simulink   ');
disp('   -> Matlab: m-file       ');
disp('   -> MONOLIX              ');
%disp('   -> NONMEM               ');
disp('   -> PottersWheel ');
perl_output = perl('build_system.pl');

if(not(strcmp(perl_output, '')))
    disp('perl script may have failed, see below:');
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
  end
else
  disp('Simulink _NOT_ found');
  disp('C-matlab targets will not work');
end
disp('-----------------------------------')

clear myversion, perl_output;
