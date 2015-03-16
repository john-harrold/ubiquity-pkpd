function str = pad_string(str, maxlength)
%  function str = padstring(str, maxlength)
%
%  adds spaces to the beginning of the string 'str' until it is lenth
%  'maxlength'
%

if(length(str) < maxlength)
  evalstr = strcat('str = sprintf(''%', mat2str(maxlength),'s'',str);');
  eval(evalstr);
end
