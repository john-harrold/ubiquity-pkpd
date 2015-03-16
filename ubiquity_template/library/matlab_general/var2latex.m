function str = var2latex(var, format) 
%  function str = var2latex(var, format) 
% takes 'var' and returns string according to 'format' and then latexify the
% resulting format. 
%
% examples:
% call                        |   returns
%-----------------------------|----------------------
% var2latex(20,'%.2e')        |   2.0\times 10^{1}
% var2latex(2.45e-2,'%.2e')   |   2.45\times 10^{-2}
% var2latex(pi,'%.2e')        |   3.14
% var2latex(sin(pi/2),'%.2e') |   1
%                             |   


% author:  john harrold <jmh@member.fsf.org>
% date:    2003.12.09
% version: 1

myerror = 0;

if ~isnumeric(var)
    disp('var must be a number');
    myerror = 1;
end

if ~ischar(format)
    disp('var must be a number');
    exit;
    myerror = 1;
end

if ~myerror 
    eval_str = sprintf('str = sprintf(''%s'', var );',format);


    eval(eval_str);
    str = latexify(str);
end


function str = latexify(str)
% exponents can have the form:
% xxe-0xx
% xxe-xx
% xxe+0xx
% xxe+xx

warning off

if regexp(str, 'e-0*\d*')
  str = regexprep(str, 'e-0*(\d*)', '\\times 10^{-$1}', 'tokenize');
end
if regexp(str,'e-\d*')
  str = regexprep(str, 'e-(\d*)',  '\\times 10^{-$1}', 'tokenize');
end
if regexp(str,'e\+0*$')
  str = regexprep(str, 'e\+0*$',     '',  'tokenize');
end
if regexp(str,'e\+0*\d*')
  str = regexprep(str, 'e\+0*(\d*)', '\\times 10^{$1}',  'tokenize');
end
if regexp(str,'e\+\d*')
  str = regexprep(str, 'e\+(\d*)',  '\\times 10^{$1}',  'tokenize');
end

if regexp(str,'\.0*$')
  str = regexprep(str, '\.0*$',  '',  'tokenize');
end

if regexp(str,'\.00*\\')
  str = regexprep(str, '\.00*\\',  '.0\',  'tokenize');
end

warning on

