function remove_islands(regularized_ciftifile, out_map_colored_single_path)
% clear
addpath(genpath('./COMBINED_UTILS'));
% regularized_ciftifile='../data/sub-MSCPI07_CONCAT0.2FD_rawassn_minsize400_regularized.dtseries.nii'
% out_map_colored_single_path='../results/matt_results/sub-MSCPI07_minsize400_final.dtseries.nii'

minsize=20;


cifti_data = ft_read_cifti_mod(regularized_ciftifile); 
neighbors = cifti_neighbors(regularized_ciftifile);

out_map_colored_single_v = ft_read_cifti_mod(out_map_colored_single_path); 
out_map_colored_single = out_map_colored_single_v.data;
allcolors = unique(out_map_colored_single);


for color = allcolors(:)'
    clusteredmetric = zeros(size(out_map_colored_single));
    thiscolorverts = find(out_map_colored_single==color);
    for vertex = thiscolorverts'
        
        %find the neighbors of this vertex
        vertexneighbors = neighbors(vertex,:);
        vertexneighbors(isnan(vertexneighbors)) = [];
        
        %find which of those neighbors also pass the thresholds
        vertexneighbors_thiscolor = intersect(thiscolorverts,vertexneighbors);
        
        %find if those neighbors have already been assigned different cluster values
        uniqueneighborvals = unique(clusteredmetric(vertexneighbors_thiscolor));
        uniqueneighborvals(uniqueneighborvals==0) = [];
        
        %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
        if isempty(uniqueneighborvals)
            clusteredmetric(vertexneighbors_thiscolor) = vertex;
            %if there is only one previous cluster identifier present, make all the neighbors that value
        elseif length(uniqueneighborvals)==1
            clusteredmetric(vertexneighbors_thiscolor) = uniqueneighborvals;
            %if there are multiple cluster identifier values in the neighborhood, merge them into one
        else
            for valuenum = 2:length(uniqueneighborvals)
                clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
            end
        end
    end
    uniqueclustervals = unique(clusteredmetric);
    uniqueclustervals(uniqueclustervals==0) = [];
    
    for clusternum = uniqueclustervals'
        if nnz(clusteredmetric==clusternum) < minsize
            neighborverts = unique(neighbors((clusteredmetric==clusternum),2:end));
            neighborverts(isnan(neighborverts)) = [];
            borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
            borderverts(out_map_colored_single(borderverts)<1) = [];
            mode_neighborval = mode(out_map_colored_single(borderverts));
            out_map_colored_single(clusteredmetric==clusternum) = mode_neighborval;
        end
    end
end


out = out_map_colored_single;


% consensusmap = out;
% cifti_data.data = consensusmap;
% if ~exist('cifti_data.mapname')
%     cifti_data.mapname = {'Column number ' num2str(1)};
%     cifti_data.dimord = 'scalar_pos';
% else
%     cifti_data.mapname = cifti_data.mapname(1);
% end
% 
% dotsloc = strfind(regularized_ciftifile,'.');
% basename = regularized_ciftifile(1:(dotsloc(end-1)-1));
% outname = [basename '_recolored'];
% ft_write_cifti_mod(outname,cifti_data);
% set_cifti_powercolors([outname '.dscalar.nii'])
% 
% 
% 


consensusmap = out;
out_map_colored_single_v.data = consensusmap;
if ~exist('cifti_data.mapname')
    out_map_colored_single_v.mapname = {'Column number ' num2str(1)};
    out_map_colored_single_v.dimord = 'scalar_pos';
else
    out_map_colored_single_v.mapname = out_map_colored_single_v.mapname(1);
end

dotsloc = strfind(out_map_colored_single_path,'.');
basename = out_map_colored_single_path(1:(dotsloc(end-1)-1));
outname = [basename '_island_getaway'];
ft_write_cifti_mod(outname,out_map_colored_single_v);
set_cifti_powercolors([outname '.dscalar.nii'])
