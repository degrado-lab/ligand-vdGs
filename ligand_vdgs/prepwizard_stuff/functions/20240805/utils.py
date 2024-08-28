import os
import shutil

def set_up_outdir(outdir, overwrite):
    '''Create outdir if it does not exist, and require the overwrite flag if it does,
    so that there are no stale files.
    NOTE: TAKE OUT RESUME FLAG because iteration over the existing files, and check_output_file_exists, 
    is handled downstream of here. '''

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise ValueError(f'The filename you designated as the output directory, {outdir}, already '
                             'exists and is not a directory.')
        if overwrite: 
            shutil.rmtree(outdir) # removes the dir and everything in it
            os.makedirs(outdir)
        else:
            # Are there files?
            files = os.listdir(outdir)
            if files:
                raise FileExistsError(
                f'The dir "{outdir}" already exists and has files. Please remove them and re-run this script, '
                'pass an overwrite flag, or modify code to handle how to resume usage of files generated from '
                'a previous run.')
                #if resume is None: 
                #    raise FileExistsError(
                #    f'The dir "{outdir}" already exists and has files. Please remove them and re-run this script, '
                #    'or pass an overwrite or resume flag.')
                #if resume==False:  # resume is an existing flag in the script
                #    raise FileExistsError(
                #    f'The dir "{outdir}" already exists and has files. Please remove them and re-run this script, '
                #    'or pass the overwrite flag, to prevent files from being forcefully overwritten. You may '
                #    'also pass the resume flag if you would like to use files generated from a previous run.')

    else: # dir doesn't exist, so create it
        parent_dir = os.path.dirname(outdir)
        if parent_dir: # skip cases where parent_dir is "" for the current working dir.
            if not os.path.exists(parent_dir):
                os.makedirs(parent_dir) 
            os.makedirs(outdir)