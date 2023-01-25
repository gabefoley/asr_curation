import os
import shutil
def test_snakemake_pipeline():

    # Move the snakefile into the test directory
    snakefile_dest = "./snakefile"
    shutil.copyfile("../snakefile", snakefile_dest)

    os.system("snakemake --cores 1 --configfile files/config/test_config.yaml -R validate_ids")


    # Remove snakefile
    # if os.path.exists(snakefile_dest):
    #     os.remove(snakefile_dest)