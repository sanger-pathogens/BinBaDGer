# tree_builder.py
import numpy as np
import dendropy
import subprocess
import os
import sys

def buildRapidNJ(phylip_path, meta_ID):
    """Use rapidNJ for more rapid tree building
    Takes a path to phylip, system calls to rapidnj executable, loads tree as
    dendropy object (cleaning quotes in node names), removes temporary files.
    """

    # construct tree
    rapidnj_cmd = "rapidnj " + phylip_path + " -n -i pd -o t -x " + meta_ID + ".raw"
    
    print(rapidnj_cmd)
    
    try:
        subprocess.run(rapidnj_cmd, shell=True, check=True)
        with open(meta_ID + ".raw", 'r') as f, open(f"{meta_ID}.nwk", 'w') as fo:
            for line in f:
                fo.write(line.replace("'", ''))
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Could not run command " + rapidnj_cmd + "; returned code: " + str(e.returncode) + "\n")
        sys.exit(1)
    
    tree = dendropy.Tree.get(path=f"{meta_ID}.nwk", schema="newick")
    return tree

def generate_phylogeny(phylip_path, meta_ID, tree_suffix, overwrite):
    """Generate phylogeny using dendropy or RapidNJ"""

    tree_filename = f"{meta_ID}.nwk"
    if overwrite or not os.path.isfile(tree_filename):
        
        sys.stderr.write("Building phylogeny\n")
        
        tree = buildRapidNJ(phylip_path, meta_ID)

        tree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
        
        tree.write(path=tree_filename, schema="newick", suppress_rooting=True, unquoted_underscores=True)
    
    else:
        sys.stderr.write("NJ phylogeny already exists; add --overwrite to replace\n")
