#!/usr/bin/env python3
"""
Validate the paths in a location file. If the paths are relative paths 
from the root of a directory in the resources dolder of this pipeline,
the paths are updated to relative paths from the pipeline root directory.
"""

from sys import argv, exit
import os

def validate_paths(location_file):
    # read the locations file
    with open(location_file, 'r') as file:
        data = file.readlines()
    
    # validate the paths
    update = False
    for i, line in enumerate(data):
          values = line.rstrip().split()
          if len(values) > 2:
            sys.stderr.write('Error: invalid line: {}\n'.format(line))
            sys.exit(-1)
          path = values[-1]
          if os.path.exists(path) and not os.path.isdir(path):
              continue
          new_path = os.path.join('resources', path)
          if os.path.exists(new_path) and not os.path.isdir(new_path):
              new_line = "{}\n".format(new_path)
              if len(values) == 2:
                  new_line = "{} {}".format(values[0], new_line)
              data[i] = new_line
              update = True
          else:
            sys.stderr.write('Error: path not found: {}\n'.format(path))
            sys.exit(-1)
            
    # update the file with new paths
    if update:
        with open(location_file, 'w') as file:
            file.writelines(data)
    

    
validate_paths(snakemake.input[0])