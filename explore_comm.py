#!/usr/bin/env python3


import os, subprocess, sys, argparse
import numpy as np

pwd = os.getcwd()
path_to_wb = '/Applications/workbench/bin_macosx64/wb_command'
template_path = '/Users/matt/workspace/code/code_projects/child_networks/surfaces/92ktemplate.dtseries.nii'

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
        

def main():
    arg_parser = argparse.ArgumentParser()
    if len(sys.argv[1:])==0:
        print('\nArguments required. Use -h option to print FULL usage.\n')
    
    arg_parser.add_argument('regularized_dt_path', type=os.path.abspath,
                            help = 'the path to the nii file with multiple community detections solution')
    arg_parser.add_argument('col_num', type=int, help = 'the column number')
    arg_parser.add_argument('output_dir',  type = os.path.abspath, help = 'the desired output directory')
    args = arg_parser.parse_args()


    regularized_dtseries = load_nii(args.regularized_dt_path)
    data_thr = regularized_dtseries[:,1+args.col_num]

    base_name = os.path.basename(args.regularized_dt_path) 
    outname = f'Column_{args.col_num}_{base_name[:10]}'
    save_nii(data_thr, outname, args.output_dir, template_path)


    values = np.unique(data_thr)
    values = values[values != 0]
    data = np.random.rand(data_thr.shape[0], 18)
    save_nii(data, f'{outname}_randvals', args.output_dir, args.regularized_dt_path)

    x = os.path.join(args.output_dir, f'{outname}.dtseries.nii')
    y = os.path.join(args.output_dir, f'{outname}.dlabel.nii')
    z = os.path.join(args.output_dir, f'{outname}_randvals.dtseries.nii')
    w = os.path.join(args.output_dir, f'{outname}.ptseries.nii')
    o = os.path.join(args.output_dir, f'{outname}.pconn.nii')


    wb_comm = ' '.join([path_to_wb, '-cifti-label-import', x, "''", y])
    subprocess.call(wb_comm, shell=True)

    wb_comm = ' '.join([path_to_wb, '-cifti-parcellate', z, y, 'COLUMN', w])
    subprocess.call(wb_comm, shell=True)

    wb_comm = ' '.join([path_to_wb, '-cifti-correlation', w, o, '-fisher-z'])
    subprocess.call(wb_comm, shell=True)

    os.remove(x)
    os.remove(y)
    os.remove(z)
    os.remove(w)


if __name__ == '__main__':
    sys.exit(main())