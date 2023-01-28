function [ NetworkMapping ] = Infomap_Assigning(template,ciftifile,outputfile)
% Infomap template matching
%   Template = string - path to template
%   ciftifile = string - path to regularized ciftifile
%   outputfile = string - path to output file

% updated for SSRDE -DVD Dec 2021

addpath(genpath('./COMBINED_UTILS'));
addpath(genpath('/sphere/greene-lab/lab_members/matt/Allvisit_Infomap_matt/COMBINED_UTILS'));

% addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/'))
% addpath(genpath('/data/nil-bluearc/GMT/Scott/MSC_Subcortical/Scripts'))
% Read in template

Template = ft_read_cifti_mod(template);
Template = Template.data(1:59412,1);
%Template(Template==6) = 10;
% Read in regularized ciftifile
Data = ft_read_cifti_mod(ciftifile);
Data = Data.data(1:59412,:);
NetworkMapping = zeros(size(Data));

% Loop through thresholds
for thr = 1:size(Data,2)
   ThisThr = Data(:,thr); 
   % In a threshold, loop through communtiy assignment values
   for p = 1:max(ThisThr)
      
      % For a community assignment value, record all indicides of its occurence in a threshold
      Idx = find(ThisThr == p);
      
      % Look at the corresponding indicies in the template map for that community assignment value
      % Take the mode / most occuring index 
      if Idx
          TemplateVerts = Template(Idx,1);
          Table = tabulate(TemplateVerts);
          [a b] = max(Table(:,2));  
          NetworkMapping(Idx,thr) = Table(b,1);
         %  NetworkMapping(Idx,thr) = mode(TemplateVerts);

      end
   end
   
end


Data = ft_read_cifti_mod(ciftifile);
Data.data(1:59412,:) = NetworkMapping;
ft_write_cifti_mod(outputfile,Data)
end

