import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_all_images():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("./all_images/data")
        expected_path = PurePosixPath("./all_images/expected")
        config_path = PurePosixPath("./config")
        additional_data_path = PurePosixPath("./additional_data")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree(config_path, workdir / "config")
        shutil.copytree(additional_data_path, workdir / "additional_data")
        # dbg
        print("all_images", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "all_images",
            "-F", 
            "-j1",
            "--keep-target-files",
            "--configfile",
            "/Users/gabefoley/Dropbox/Code/Python_Workspace/ASR_Curation/asr_curation/config/test_config.yaml",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
