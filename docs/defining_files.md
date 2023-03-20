# Defining input files


## **Defining the config files**

## **Defining the subset files**

Subset files need to be in the format 

`<subset_name> = <column_name> : <string> $ <other_column_name> : NOT <string>`

Where `subset_name` can be any name you wish to call the subset, `column_name` and `other_column_name` are columns that appear in the final annotation file (see /csv/custom/*.csv) and `string` is a term that you wish to include or exclude.

For each row / sequence in your annotation file, if the given column contains the string (or substring) then it will be included.

If you use the `NOT` modifier, for each row / sequence in your annotation file, if the given column contains the string (or substring) then it will not be included.

You can define multiple subsets within the one file.



```
# THIS IS A COMMENT LINE.

eukaryotic_kari = lineage_superkingdom : Eukaryota $ protein_name : NOT Deleted, Merged $ Non_AA_Character : FALSE

classI_kari = KARI_Class : Class_1 $ protein_name : NOT Deleted, Merged $ Non_AA_Character : FALSE

all_kari = *
```


## Including all sequences in a subset
The default way to include every row / sequence is 

`<subset_name> = *`

For example,

`all = *`  

