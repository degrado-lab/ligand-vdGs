import os
import argparse
import subprocess

def parse_args():
    '''
    Parses command-line arguments.
    '''
    argp = argparse.ArgumentParser()
    argp.add_argument('--probe-path', help="Path to probe executable.")
    argp.add_argument('--input-pdb', help="Path to input PDB to run Probe on.")
    argp.add_argument('--outdir', help="Directory path for outputting probe files.")
    return argp.parse_args()

def main():
    args = parse_args()
    _probe = args.probe_path
    input_pdb = args.input_pdb
    outdir = args.outdir

    basename = os.path.basename(input_pdb)
    for ext in (".pdb.gz", ".pdb"):
        if basename.endswith(ext):
            basename = basename[:-len(ext)]
            break
    pdb_code = basename[-4:]
    out_filename = pdb_code + '.probe'
    out_subdir = os.path.join(outdir, pdb_code[1:3])
    out_path = os.path.join(out_subdir, out_filename)
    if not os.path.exists(out_subdir):
        os.makedirs(out_subdir)
          
    cmd = f'{_probe} -U -CON -WEAKH -DE32 -4 -Both ALL ALL {input_pdb} > {out_path}'
    print('-'*62)
    print('Probe command:')
    print(cmd)
    
    result = subprocess.run(cmd, shell=True)

    if os.path.exists(out_path):
        gzip_cmd = ['gzip', out_path]
        print('-'*62)
        print('gzip command:')
        print(' '.join(gzip_cmd))
        subprocess.run(gzip_cmd, check=True)
        print('-'*62)
        print('Job completed.')
    else:
        print('-'*62)
        print('The expected pdb file was not output: ', out_path)

if __name__ == '__main__':
    main()