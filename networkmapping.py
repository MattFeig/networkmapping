#!/usr/bin/env python3

__doc__ = """This script creates a consensus InfoMap solution. It combines InfoMaps ran on cortical resting state fMRI 
data, done at multiple sparsity thresholds, in a semi-intelligent way. It then matches the consensus communities to the
nearest match on a given template network map.
"""

import os, subprocess, sys, argparse
import numpy as np

pwd = os.getcwd()

def load_nii(nii_path, purge=True):
    '''
    Load in a .nii file to work with.
    
    This function converts the .nii to a text file using workbench, and then
    uses numpy to load in the text file.
    '''

    path_to_wb = '/Applications/workbench/bin_macosx64/wb_command'

    end = nii_path.find('.')
    outname = nii_path[:end] + '.txt'
    wb_comm = ' '.join([path_to_wb, '-cifti-convert -to-text', nii_path, outname])
    subprocess.call(wb_comm, shell=True)
    nii_data = np.loadtxt(outname)
    
    if purge:
        os.remove(outname)
    
    return nii_data

def save_nii(array, output_name, output_dir_path, wb_required_template_path, purge = True):
    '''
    Save a numpy array as a .nii file, utilizing workbench as an intermediary.

    Arguments:
        array: the numpy array to save
        output_name: the name to use for the output file
        output_dir_path: the directory path to save the output file in
        wb_required_template_path: a dtseries of the same dimension used as a template to write over
        purge: if True, delete the intermediate text file after creating the .nii file
    '''

    # Define the path to workbench
    path_to_wb = '/Applications/workbench/bin_macosx64/wb_command'
    
    # Check if the output directory exists
    if not os.path.isdir(output_dir_path):
        raise Exception(f"The output folder {output_dir_path} does not exist")
    
    # Define the output path and the name for the intermediate .txt file
    out_path = os.path.join(output_dir_path, output_name)
    outnamecifti = out_path+'.dtseries.nii'
    
    # Check if the output file already exists and print a warning if it does
    if os.path.isfile(outnamecifti):
        print('-WARNING: Overwriting')
    
    # Save the numpy array as a text file
    np.savetxt(out_path, array)
    
    # Use workbench to convert the text file to a .nii file using the provided template
    wb_comm = ' '.join([path_to_wb, '-cifti-convert -from-text', out_path, wb_required_template_path, outnamecifti])
    subprocess.call(wb_comm, shell=True)
    
    # Delete the intermediate text file if purge is True
    if purge:
        os.remove(out_path)

def sparsest_template_match(regularized_dtseries_cortex, template_cortex, dthr = .1, man_edit_array = np.array([[-99,-99]]), 
                            force = False):

    '''To do: Refactor and add docstring'''

    out_map_colored_single = np.zeros(91282)

    for i, thr in enumerate(regularized_dtseries_cortex.T[::-1]): # Loop through community solutions as different sparsity thresholds
        
        # Idx_missing = np.where(thr < 1)[0]
        
        # for idx in Idx_missing:
        #     thisassignments = regularized_dtseries_cortex[idx,:]
        #     nomissing = np.array([x for x in thisassignments if x >=1])
        #     if nomissing.shape[0]>=1:
        #         thr[idx] = nomissing[0]

        for p in np.arange(-1,max(thr)+1): # Loop through community assignments in the current threshold
            if p in thr:
                A = (thr == p)
                Idx = np.where(thr == p)[0]
                D_list = []

                net_colors = range(1,18)
                for templatenet in net_colors: # Loop through all the networks in a given template and calculate dice overlap   
                    B = (template_cortex.astype(int) == templatenet)
                    D = np.logical_and(A,B).sum()/np.logical_or(A,B).sum()  
                    D_list.append(D)
                potential_matches = np.where(np.array(D_list)>dthr)[0]

                # if (len(potential_matches) > 1) and (6 in potential_matches):
                #     # If there are multiple potential matches, and premotor is in this list, set it to 0 for now.
                #     D_list[6] = 0

                if p not in man_edit_array[:,0]:
                    # If there are no manual edits required, assign the largest overlap >.1 dice to the final output
                    if np.max(D_list) > dthr:
                        out_map_colored_single[Idx] = net_colors[np.argmax(D_list)]
                else:
                    if force == True:    
                        p_row_ind = np.where(man_edit_array[:,0]==p)[0][0]
                        man_assignment = man_edit_array[p_row_ind,1]
                        out_map_colored_single[Idx] = man_assignment
                    else:
                        potential_matches = np.where(np.array(D_list)>dthr)[0]
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
    arg_parser = argparse.ArgumentParser(description="""this script creates a consensus InfoMap solution and then performs template matching. 
                example call: python3 networkmapping.py data/infomap_threshold_maps.dtseries.nii consensusmatchedoutput""")
    if len(sys.argv[1:])==0:
        print('\nArguments required. Use -h option to print FULL usage.\n')
    
    arg_parser.add_argument('regularized_dt_path', type=os.path.abspath,
                            help = '''the path to the nii file with multiple community detections solution 
                            to be fed into the consensus algorithm''')
    arg_parser.add_argument('output_name', type = str, help = 'the desired name of the output file')
    arg_parser.add_argument('-output_dir', default = 'results', type = os.path.abspath, required= False,
                            help = 'the desired output directory, the default value is ./results', dest = 'output_dir')
    arg_parser.add_argument('-d', default = .1 , type = float, required= False,
                           help = 'dice threshold, default = .1', dest = 'dthr')
    arg_parser.add_argument('-s', default = 0 , type = int, required= False,
                           help = 'thresholds to skip', dest = 'skip')
    arg_parser.add_argument('-t', action='store', default = 'data/Networks_template3.dscalar.nii', type=os.path.abspath, required=False,
                            help= '''the path the desired network organizition for template matching, the default path
                            is data/Networks_template3.dscalar.nii''', dest = 'net_template_path')
    arg_parser.add_argument('-w', action='store', default = 'data/surfaces/92ktemplate.dtseries.nii', type=os.path.abspath,
                            required=False, help='''the path to the same dimension .nii to use as template for saving,
                            default is data/surfaces/92ktemplate.dtseries.nii''', dest = 'wb_required_template_path')
    arg_parser.add_argument('--nocleanup', action = 'store_true', default = False, required=False,
                            help = 'use this flag to skip the island removal step, helpful if you dont have matlab', dest = 'skip_cleanup')
    args = arg_parser.parse_args()

    final_path = os.path.join(args.output_dir,args.output_name)+'.dtseries.nii'
    net_template_data = load_nii(args.net_template_path)
    regularized_dtseries = load_nii(args.regularized_dt_path)
    template_cortex = net_template_data[0:59412]
    regularized_dtseries_cortex = regularized_dtseries[0:59412, args.skip:]
    out_map_colored_single = sparsest_template_match(regularized_dtseries_cortex,template_cortex, dthr=args.dthr)
    save_nii(out_map_colored_single, args.output_name, args.output_dir, args.wb_required_template_path)

    if args.skip_cleanup == False:
        regularized_dt_path_mat_input = os.path.join('..',args.regularized_dt_path)
        outmap_single_mat_input = os.path.join('..',final_path)
        matlab_cleanup(regularized_dt_path_mat_input, outmap_single_mat_input)
        os.remove(final_path)

if __name__ == '__main__':
    sys.exit(main())