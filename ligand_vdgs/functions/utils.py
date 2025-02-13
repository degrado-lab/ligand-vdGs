import os
import sys
import shutil

def handle_existing_files(out_dir):
   os.makedirs(out_dir, exist_ok=True)
   if len(os.listdir(out_dir)) > 0:
      raise ValueError(f'The output dir {out_dir} is not empty. Remove files or define a new '
            'output dir name to prevent accidental overwriting.')

def set_up_outdir(outdir, overwrite):
    '''Create outdir if it does not exist, and require the overwrite flag if it does,
    so that there are no stale files.'''

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise ValueError(f'The filename you designated as the output directory, {outdir}, already '
                             'exists and is not a directory.')
        if overwrite: 
            shutil.rmtree(outdir) # removes the dir and everything in it
            os.makedirs(outdir)
    
    else: # dir doesn't exist, so create it
        parent_dir = os.path.dirname(outdir)
        if parent_dir: # skip cases where parent_dir is "" for the current working dir.
            if not os.path.exists(parent_dir):
                os.makedirs(parent_dir) 
            os.makedirs(outdir)

def valid_database_subdir_format(input_dir):
    ''' Ensure that the pdb database dir has subdirs formatted similarly to the RCSB
    mirror format (see docs/database_generation.txt). '''
    
    def error_message():
        print('Database structure must be similar to the RCSB mirror format; '
              'see docs/database_generation.txt')
    
    # Input path must be a dir
    if not os.path.exists(input_dir):
        error_message()
        return False
    # Input dir must not be empty
    if not len(os.listdir(input_dir)) > 0: 
        error_message()
        return False
    # Ensure that there is at least one correctly formatted subdir
    for p in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, p)
        if os.path.isdir(subdir_path):
            # Subdirs should be the inner 2 characters of a 4-char pdb file name
            if len(p) == 2: 
                for pdb in os.listdir(subdir_path):
                    if pdb[1:3] == p:
                        return True
    # No valid subdir was found
    error_message()
    return False