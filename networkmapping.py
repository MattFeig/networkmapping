#!/usr/bin/env python3

__doc__ = """This script creates a consensus InfoMap solution. It combines InfoMaps ran on cortical resting state fMRI 
data, done at multiple sparsity thresholds, in a semi-intelligent way. It then matches the consensus communities to the
nearest match on a given template network map.
"""

import os, subprocess, sys, argparse
import numpy as np

pwd = os.getcwd()

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

    for i, thr in enumerate(regularized_dtseries_cortex.T[::-1]): # Loop through community solutions as different sparsity thresholds
        Idx_missing = np.where(thr < 1)[0]
        
        for idx in Idx_missing:
            thisassignments = regularized_dtseries_cortex[idx,:]
            nomissing = np.array([x for x in thisassignments if x >=1])
            if nomissing.shape[0]>=1:
                thr[idx] = nomissing[0]
        for p in np.arange(1,max(thr)+1): # Loop through community assignments in the current threshold
            if p in thr:
                A = (thr == p)
                Idx = np.where(thr == p)[0]
                D_list = []

                for templatenet in np.unique(template_cortex.astype(int)): # Loop through all the networks in a given template
                    B = (template_cortex.astype(int) == templatenet)
                    D = np.logical_and(A,B).sum()/np.logical_or(A,B).sum()  
                    # Calculate Dice overlap of the current community and current template network 
                    D_list.append(D)
                potential_matches = np.where(np.array(D_list)>.1)[0]
                if (len(potential_matches) > 1) and (6 in potential_matches):
                    # If there are multiple potential matches, and premotor is in this list, set it to 0 for now.
                    D_list[6] = 0
                if p not in man_edit_array[:,0]:
                    # If there are no manual edits required, assign the largest overlap >.1 dice to the final output
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

def matlab_cleanup(regularized_dt_path, conensus_communties_path):
    '''Calls matlab tools remove_islands.m, requries the orginal path of communty solutions across thresholds,
    along with the final consensus path'''

    matlab_tools_path = ('./matlab_tools')

    regularized_dt_path_mat_input = f"'{regularized_dt_path}'"
    conensus_communties_path_mat_input = f"'{conensus_communties_path}'"
    os.chdir(matlab_tools_path)
    matlab_call = f""" matlab -nodesktop -nodisplay -nosplash -r "remove_islands({regularized_dt_path_mat_input}, {conensus_communties_path_mat_input});exit" """
    subprocess.call(matlab_call, shell=True)

def main():

    arg_parser = argparse.ArgumentParser()
    if len(sys.argv[1:])==0:
        print('\nArguments required. Use -h option to print FULL usage.\n')
    
    arg_parser.add_argument('regularized_dt_path', type=os.path.abspath,
                            help = '''the path to the nii file with multiple community detections solution 
                            to be fed into the consensus algorithm''')
    arg_parser.add_argument('output_name', type = str, help = 'the desired name of the output file')
    arg_parser.add_argument('-output_dir', default = 'results', type = os.path.abspath, required= False,
                            help = 'the desired output directory, the default value is ./results', dest = 'output_dir')
    arg_parser.add_argument('-t', action='store', default = 'data/Networks_template.dscalar.nii', type=os.path.abspath, required=False,
                            help= '''the path the desired network organizition for template matching, the default path
                            is data/Networks_template.dscalar.nii''', dest = 'net_template_path')
    arg_parser.add_argument('-nocleanup', action = 'store_true', default = False, required=False,
                            help = 'use this flag to skip the island removal step, helpful if you dont have matlab', dest = 'skip_cleanup')
    arg_parser.add_argument('-w', action='store', default = 'data/surfaces/92ktemplate.dtseries.nii', type=os.path.abspath,
                            required=False, help='''the path to the same dimension .nii to use as template for saving,
                            default is data/surfaces/92ktemplate.dtseries.nii''', dest = 'wb_required_template_path')
    args = arg_parser.parse_args()

    regularized_dt_path = args.regularized_dt_path
    net_template_path = args.net_template_path
    wb_required_template_path = args.wb_required_template_path
    skip_cleanup = args.skip_cleanup
    output_name = args.output_name
    output_dir = args.output_dir

    final_path = os.path.join(output_dir,output_name)+'.dtseries.nii'
    net_template_data = load_nii(net_template_path)
    regularized_dtseries = load_nii(regularized_dt_path)
    template_cortex = net_template_data[0:59412]
    regularized_dtseries_cortex = regularized_dtseries[0:59412,:]
    out_map_colored_single = sparsest_template_match(regularized_dtseries_cortex, template_cortex)
    save_nii(out_map_colored_single, output_name, output_dir, wb_required_template_path)

    if skip_cleanup == False:
        regularized_dt_path_mat_input = os.path.join('..',regularized_dt_path)
        outmap_single_mat_input = os.path.join('..',final_path)
        matlab_cleanup(regularized_dt_path_mat_input, outmap_single_mat_input)

if __name__ == '__main__':
    sys.exit(main())