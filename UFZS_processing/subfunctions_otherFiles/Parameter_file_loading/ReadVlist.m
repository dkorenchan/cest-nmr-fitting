% ReadVlist: Reads in values from a Bruker variable list file (e.g. valist, 
% vdlist, etc) and outputs them as a numerical array
%
%   INPUTS:
%       path        -   String containing path to list file
%       filename    -   String containing name of list file 
%
%   OUTPUTS:
%       n           -   Array of type double containing values read from
%                       list file. Time values are returned in units of 
%                       seconds.
%
function n=ReadVlist(path,filename)
%Read in the sf, offset, sw_p, xdim and si pars from fname
%xwinnmr parameter file
%DK edits 4/5/21: Fixed bug - if there was a space at end of line, it
%would return NaN 
fname = fullfile(path,filename);
fin = fopen (fname, 'r');
n = [];
while (~feof(fin)) 
    line = fgetl(fin);
    % Search for any ending time unit (m, u) for numbers in list, convert to s
    % DK edit 4/5/21: Change endsWith to contains - fixes problem if line
    % ends in a space
    if contains(line,'m')
        num = str2double(extractBefore(line,'m'));
        num = num / 1000;
    elseif contains(line,'u')
        num = str2double(extractBefore(line,'u'));
        num = num / 1000000;  
    % end DK edits 4/5/21
    else
        num = str2double(line);
    end
    n = [n num];  
end
fclose(fin);
end