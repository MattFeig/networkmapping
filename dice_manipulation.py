import os, subprocess, sys, argparse
import numpy as np

values = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]

def main():
    

    arg_parser = argparse.ArgumentParser(description="""this script creates a consensus InfoMap solution and then performs template matching. 
                example call: python3 networkmapping.py data/infomap_threshold_maps.dtseries.nii consensusmatchedoutput""")
    if len(sys.argv[1:])==0:
        print('\nArguments required. Use -h option to print FULL usage.\n')
    
    arg_parser.add_argument('regularized_dt_path', type=os.path.abspath,
                            help = '''the path to the nii file with multiple community detections solution 
                            to be fed into the consensus algorithm''')
    
    arg_parser.add_argument('output_name', type = str, help = 'the desired name of the output file')
    arg_parser.add_argument('-output_dir', default = 'results', type = os.path.abspath, required= True,
                            help = 'the desired output directory, the default value is ./results', dest = 'output_dir')
    arg_parser.add_argument('-t', action='store', default = 'data/Networks_template3.dscalar.nii', type=os.path.abspath, required=False,
                            help= '''the path the desired network organizition for template matching, the default path
                            is data/Networks_template3.dscalar.nii''', dest = 'net_template_path')
    arg_parser.add_argument('--nocleanup', action = 'store_true', default = False, required=False,
                            help = 'use this flag to skip the island removal step, helpful if you dont have matlab', dest = 'skip_cleanup')
    args = arg_parser.parse_args()

    for value in values:
        output_name2 = f'{args.output_name}_{value}dthr'
        run_command=f'python3 networkmapping.py {args.regularized_dt_path} {output_name2} -output_dir {args.output_dir} -t {args.net_template_path} -d {value}'
        print(f'Running command: {run_command} \n')
        subprocess.call(run_command.split())
        

if __name__ == '__main__':
    sys.exit(main())