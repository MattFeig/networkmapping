function [Ci] = run_infomap_on_pajekfile_SSRDE(pajekfilename,reps)
%[Ci] = run_infomap_on_pajekfile(pajekfilename,reps)
%
%
% This script runs infomap on a pajekfile with some number of
% repetitions. It then returns the community assignments found.
%
% Updated to run on SSRDE - DVD Dec 2021

% infomapfolder = '/data/nil-bluearc/GMT/Evan/Scripts';
infomapfolder = '/sphere/greene-lab/Shared_Tools';

% this will be the relevant output of infomap
[pathstr,cluname,ext] = filenamefinder(pajekfilename,'dotsout');
clufile = [ pathstr '/' cluname '.clu' ];

% obtain seed #
clear randnum;
randnum=ceil(rand*1000000);

disp('running infomap...')
% run infomap
[failed, message] = system([infomapfolder '/Infomap-0.15.7/Infomap --clu -2 -s' num2str(randnum) ' -N' num2str(reps) ' ' pajekfilename ' ' pathstr]);
if logical(failed)
    disp(message)
end

% does this hang up the job???
% So parfor doesn't crap out
isclufile = exist(clufile);
while isclufile == 0
    pause(60)
    isclufile = exist(clufile);
end

Ci = textread(clufile,'%d','headerlines',1);
end
