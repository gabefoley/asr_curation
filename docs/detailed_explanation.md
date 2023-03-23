## **Snakemake diagram**
The following diagram is the overall steps that have to be run by `snakemake` to generate the example data.

!!! note "Diagram"
    Note that in this case, the KARI data had three subsets defined, so `create_subsets` and all of the other steps downstream from this step are run three times, while the ALS data only has one subset defined and so `create_subsets` is only run once.

![snakemake directed acyclic graph](images/dag.png)


