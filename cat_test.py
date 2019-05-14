import subprocess
import os
import sys
from sys import argv
import re


def cat_and_move(input_loc, output_loc):
    """Combines .fsa_nt files in given dir and sends to specified location"""
    # need to check naming conventions

    output_name = (output_loc.split("/"))[-1]  # name of the outfile
    output_parent = os.path.abspath(os.path.join(output_loc, os.pardir))  # outfile's parent directory

    print("Merging all files in the folder " + input_loc + " ending in .fsa_nt into "\
     + output_name + " at " + output_loc)

    list_to_cat = []
    for i in os.listdir(input_loc):
        if (i.find(".fsa_nt") != -1): # requires .fsa_nt (NCBI) convention
            list_to_cat.append(i)

    my_cmd = (['cat'] + list_to_cat)

    with open(f"{output_loc}", "w") as outfile:
        print("Merging " + str(len(list_to_cat)) + " files...")
        try:
            subprocess.run(my_cmd,
                        stdout=outfile,
                        cwd=f"{input_loc}",
                        check=True)
            print("Success!")
        except subprocess.CalledProcessError as e:
            sys.exit(e)


cat_and_move(argv[1], argv[2])
