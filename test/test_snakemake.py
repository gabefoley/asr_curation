import os
import shutil
import pytest


@pytest.fixture
def test_dir(tmpdir):
    # Create a temporary directory for testing
    test_dir = tmpdir.mkdir("test_workflow")
    # Copy the Snakefile and test/files directory to the temporary directory
    snakefile_dest = os.path.join(test_dir, "snakefile")
    scripts_dest = os.path.join(test_dir, "scripts")
    test_dest = os.path.join(test_dir, "test")

    shutil.copyfile("snakefile", snakefile_dest)
    shutil.copytree("scripts", scripts_dest)
    shutil.copytree("test", test_dest)

    os.chdir(test_dir)

    return test_dir


def test_snakemake_pipeline(test_dir):
    print("/n Current workdir")
    print(os.getcwd())
    print(os.listdir())

    subsets_dir = os.path.join(
        test_dir,
        "test/files/test_workflow/example_workflow/datasets/test_ec_1_1_1_86/subsets"
    )

    print (subsets_dir)


    # Set the coverage directory to be the temporary directory defined in test_dir()
    os.system("snakemake --cores 1 --configfile test/files/config/test_config.yaml")

    assert os.path.exists(
        os.path.join(
            test_dir,
            "test/files/test_workflow/example_workflow/datasets/test_ec_1_1_1_86/subsets/all_kari/1/grasp_results/GRASP_ancestors.fa",
        )
    )

    # Delete the subsets directory after the test
    if os.path.exists(subsets_dir):
        shutil.rmtree(subsets_dir)

    print (subsets_dir)
