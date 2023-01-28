#!/usr/bin/env python3

__doc__ = """This script creates a consensus InfoMap solution. It combines InfoMaps ran on cortical resting state fMRI 
data, done at multiple sparsity thresholds, in a semi-intelligent way. It then matches the consensus communtiues to the
nearest match on a given template network map.
"""

import os, subprocess, sys, argparse
import numpy as np
import matplotlib.pyplot as plt

pwd = os.getcwd()
matlab_tools_path = ('./matlab_tools')

def load_nii(nii_path, purge = True):
    ''' Load in a .nii file to work with. 
    
    In lieu of a proper package to interface with .nii's in python, 
    this function, calls workbench to convert the .nii to text, and then
    uses numpy to load in the text file'''

    path_to_wb = '/Applications/workbench/bin_macosx64/wb_command'

    end = nii_path.find('.')
    outname = nii_path[:end] + '.txt'
    wb_comm = ' '.join([path_to_wb, '-cifti-convert -to-text', nii_path, outname])
    subprocess.call(wb_comm, shell=True)
    nii_data = np.loadtxt(outname)
    if purge == True:
        os.remove(outname)
    return nii_data

def save_nii(array, output_name, output_dir_path, wb_required_template_path, purge = True):
    '''Save a numpy array as a .nii file, utilizing workbench as an intermediary
    Argument: 
        wb_required_template_path is a dtseries of the same dimension used a template to write over'''

    path_to_wb = '/Applications/workbench/bin_macosx64/wb_command'
    
    if not os.path.isdir(output_dir_path):
        raise Exception(f"The output folder {output_dir_path} does not exist")
    out_path = os.path.join(output_dir_path,output_name)
    outnamecifti = out_path+'.dtseries.nii'
    if os.path.isfile(outnamecifti):
        print('-WARNING: Overwriting')
    np.savetxt(out_path, array)
    wb_comm = ' '.join([path_to_wb, '-cifti-convert -from-text', out_path, wb_required_template_path, outnamecifti])
    subprocess.call(wb_comm, shell=True)
    if purge == True:
        os.remove(out_path)

def sparsest_template_match(regularized_dtseries_cortex, template_cortex, man_edit_array = np.array([[-99,-99]]), 
                            force = False):

    '''To do: Refactor and add docstring'''

    out_map_colored_single = np.zeros(91282)

    for i, thr in enumerate(regularized_dtseries_cortex.T[::-1]):    
        Idx_missing = np.where(thr < 1)[0]
        for idx in Idx_missing:
            thisassignments = regularized_dtseries_cortex[idx,:]
            nomissing = np.array([x for x in thisassignments if x >=1])
            if nomissing.shape[0]>=1:
                thr[idx] = nomissing[0]
        # Loop through infomap community assignments in the current threshold
        for p in np.arange(1,max(thr)+1):
            if p in thr:
                # Indicies of the current community assignment
                A = (thr == p)
                Idx = np.where(thr == p)[0]
                D_list = []
                # Loop through all the networks in a given template
                for templatenet in np.unique(template_cortex.astype(int)):
                    # Indicies of the current template network
                    B = (template_cortex.astype(int) == templatenet)
                    networkinds = np.where(template_cortex.astype(int) == templatenet)[0]
                    # Calculate Dice index/overlap of the corrent community assigntment and the current template network 
                    D = np.logical_and(A,B).sum()/np.logical_or(A,B).sum()
                    D_list.append(D)
                # If the dice index of the current community assignemnt is large w.r.t any template network 
                # Choose the template network with the largest dice index 
                # Set the final output vertices to be equal to that template network, at the indicies 
                # of the current community assignment             
                # D_list[6] = 0
                potential_matches = np.where(np.array(D_list)>.1)[0]
                if (len(potential_matches) > 1) and (6 in potential_matches):
                    D_list[6] = 0
                if p not in man_edit_array[:,0]:
                    if np.max(D_list) > .1:
                        out_map_colored_single[Idx] = np.argmax(D_list)
                else:
                    if force == True:    
                        p_row_ind = np.where(man_edit_array[:,0]==p)[0][0]
                        man_assignment = man_edit_array[p_row_ind,1]
                        out_map_colored_single[Idx] = man_assignment
                    else:
                        potential_matches = np.where(np.array(D_list)>.1)[0]
                        if len(potential_matches) == 0:
                            pass
                        elif len(potential_matches) == 1:
                            out_map_colored_single[Idx] = potential_matches[0]
                            print(potential_matches[0] == np.argmax(D_list))
                        elif len(potential_matches) > 1:
                            p_row_ind = np.where(man_edit_array[:,0]==p)[0][0]
                            man_assignment = man_edit_array[p_row_ind,1]
                            if man_assignment in potential_matches:
                                out_map_colored_single[Idx] = man_assignment
                            else:
                                out_map_colored_single[Idx] = np.argmax(D_list)
                                print(np.argmax(D_list) in potential_matches)
                        else: 
                            print("error")

    return out_map_colored_single

def main():

    arg_parser = argparse.ArgumentParser()
    if len(sys.argv[1:])==0:
        print('\nArguments required. Use -h option to print FULL usage.\n')
    arg_parser.add_argument('regularized_dt_path', type=os.path.abspath, 
                            help='Path to the regularized dt series')
    arg_parser.add_argument('output_name', type = str, help = 'Naming of output file')
    arg_parser.add_argument('-output_dir', type = os.path.abspath, required=False, default = 'results', help = 'Naming of output file', dest = 'output_dir')
    arg_parser.add_argument('-t', action='store', default = 'data/Networks_template.dscalar.nii', type=os.path.abspath,
                            required=False, help='Path to the network template', dest = 'net_template_path')
    arg_parser.add_argument('-nocleanup', action = 'store_false', default = True, required=False,
                            help = 'Skip removal of islands from final consensus map', dest = 'island_cleanup')
    arg_parser.add_argument('-w', action='store', default = 'data/surfaces/92ktemplate.dtseries.nii', type=os.path.abspath,
                            required=False, help='Path to the same dimension .niii to use as template to overwrite, for saving purposes',
                            dest = 'wb_required_template_path')
    args = arg_parser.parse_args()

    regularized_dt_path = args.regularized_dt_path
    net_template_path = args.net_template_path
    wb_required_template_path = args.wb_required_template_path
    island_cleanup = args.island_cleanup
    output_name = args.output_name
    output_dir = args.output_dir

    final_path = os.path.join(output_dir,output_name)+'.dtseries.nii'
    net_template_data = load_nii(net_template_path)
    regularized_dtseries = load_nii(regularized_dt_path)
    template_cortex = net_template_data[0:59412]
    regularized_dtseries_cortex = regularized_dtseries[0:59412,:]
    out_map_colored_single = sparsest_template_match(regularized_dtseries_cortex, template_cortex)
    save_nii(out_map_colored_single, output_name, output_dir, wb_required_template_path)

    if island_cleanup == True:
        regularized_dt_path_mat_input = os.path.join('..',regularized_dt_path)
        regularized_dt_path_mat_input = f"'{regularized_dt_path_mat_input}'"
        outmap_single_mat_input = os.path.join('..',final_path)
        outmap_single_mat_input = f"'{outmap_single_mat_input}'"
        os.chdir(matlab_tools_path)
        matlab_call = f""" matlab -nodesktop -nodisplay -nosplash -r "remove_islands({regularized_dt_path_mat_input}, {outmap_single_mat_input});exit" """
        subprocess.call(matlab_call, shell=True)
        # os.chdir('..')

if __name__ == '__main__':
    sys.exit(main())