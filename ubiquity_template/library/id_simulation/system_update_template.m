function []=system_update_template()
% function []=system_update_template()
%
% To update the template files to the latest version simply run this script
% from the top level directory of the template.

disp(perl('library/ubiquity/update_tempalte.pl'))
