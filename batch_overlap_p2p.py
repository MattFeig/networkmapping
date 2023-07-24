import os, subprocess, sys, argparse
import numpy as np

# mscpi_sublist = ['MSCPI05', 'MSCPI06', 'MSCPI07', 'MSCPI08', 'MSCPI10', 'MSCPI11', 'MSCPI12', 'MSCPI13',
#                   'MSCPI14', 'MSCPI15', 'MSCPI17', 'MSCPI18', 'MSCPI19']

# msc_sublist = ['MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC08', 'MSC09', 'MSC10']

from_list = ['MSCPI05']

to_list = ['MSCPI06', 'MSCPI11', 'MSCPI12', 'MSCPI13', 'MSCPI14', 'MSCPI15', 'MSCPI17', 'MSCPI18', 'MSCPI19','MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC08', 'MSC09', 'MSC10']

def main():
    output_dir = '/Users/matt/workspace/code/code_projects/child_networks/data/dice_manipulation_data/mscpi_to_other'
    
    for sub1 in from_list:

        sub1_communitites_path = f'/Users/matt/workspace/code/code_projects/child_networks/data/communities_data/mscpi_regularized/sub-{sub1}_CONCAT_6.0mm_SMOOTHED_0.2FD_rawassn_minsize10_regularized.dtseries.nii'
        if os.path.exists(sub1_communitites_path):
             pass
        else: sub1_communitites_path = f'/Users/matt/workspace/code/code_projects/child_networks/data/communities_data/msc_regularized/sub-{sub1}_CONCAT_6.0mm_SMOOTHED_0.2FD_rawassn_minsize10_regularized.dtseries.nii'

        for sub2 in to_list:

            outname = f'{sub1}_to_{sub2}'

            sub2_template_path = f'/Users/matt/workspace/code/code_projects/child_networks/data/point_one_networks/MSCPI_results/abcd_template/{sub2}_island_getaway.dscalar.nii'
            if os.path.exists(sub2_template_path):
                pass
            else: sub2_template_path = f'/Users/matt/workspace/code/code_projects/child_networks/data/point_one_networks/MSC_results/template3/{sub2}_island_getaway.dscalar.nii'

            run_command=f'python3 overlap_manipulation.py {sub1_communitites_path} {outname} -output_dir {output_dir} -t {sub2_template_path}'
            print(f'Running command: {run_command} \n')
            subprocess.call(run_command.split())
        
if __name__ == '__main__':
    sys.exit(main())


