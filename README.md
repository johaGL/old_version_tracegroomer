# t4DIMet

Transform -tracer metabolome- xlsx file to csv files ready for DIMet analysis


## Requirements

The requirements are the same as in DIMet https://github.com/cbib/DIMet
First clone DIMet and create dimet conda environment (see https://github.com/cbib/DIMet), then clone t4DIMet. Both work with the same conda environment
```
conda activate dimet
```

## convert

Its execution takes only few seconds! Below we explain how to proceed.

### Input files

Three types of Tracer data are accepted by our `convert` module:
1. IsoCor results (.xlsx measurments file).
2. Results provided by the VIB Metabolomics Expertise Center (El-Maven results are shaped by VIB MEC team into a multi-sheet .xlsx file).  
3. A 'generic' .xlsx measurments file.

Your data is expected to correspond to one of the three types above. We provide you 'toy' examples of each one of them  [here](examples/readme_examples.md). Pick the example most suited to your data.

For example, let's take [toy2 example](examples/toy2/) , which contains:

 - 'data/' folder, with:
    * the tracer metabolomics data as one .xlsx file. 
    * the samples description that we call the "metadata", see below
    
 - 'analysis001/' folder, with:
     * your configuration file (extension .yml), see below.


Regarding the **metadata**, we explain it in detail in the section [Metadata](#the-metadata).


Regarding the .yml file, we supply examples that you can use as template, such  [this config.yml template](examples/toy1/analysis001/config-1-001.yml). To double-check your modifications there exist online editors, such as https://yamlchecker.com/, just copy-paste and edit!


### Execute `convert` 

The examples serve to demonstrate how fast this module can be. Take toy2 example, copy and paste the entire toy2 folder in your 'home/' folder, then from terminal:
```
python -m t4DIMet.convert toy2/analysis01/config-3-01.yml
```

--------------------
# Details

### Advanced prepare options

We provide advanced options for `convert` module, check the help:
```
python -m t4DIMet.convert --help
```
they appear as 'optional arguments' in the help menu.


You can:

- normalize by the amount of material (number of cells, tissue weight): setting the path to the file in your **.yml** configuration. The file must be like [this csv file](examples/toy2/data/nbcells-or-amountOfMaterial.csv), and the first column must contain the same names as in metadata 'former\_name'.
- normalize by an internal standard (present in your data) at choice: using the advanced option `--use_internal_standard`.

However we have some indications that can slightly differ for [users having VIB results as input](#users-having-vib-results), [users having IsoCor results](#users-having-isocor-results) or users having ['generic' type of data](#users-having-generic-data).


 
### Users having IsoCor results

A typical IsoCor results table is described in: https://isocor.readthedocs.io/en/latest/tutorials.html#output-files
 
 Our pipeline transforms its columns into tables, so here the 'Isocor column : DIMet table' correspondence:
    - corrected_area : isotopologuesAbsolute  
    - isotopologue_fraction : isotopologuesProportions
    - mean\_enrichment :  meanEnrich\_or_fracContrib
    - (nothing)   : Abundance 

Abundance table is the sum (per metabolite) of IsotopologuesAbsolute, we also perform it automatically, as this column is not present in the input data.     
    
Please stick to the example [toy1](examples/toy1/) for the names of the tables in the **.yml** file for isocorOutput
    
        
Options regarding to detection limit (LOD) and blanks will not have any effect on the IsoCor type of data: LOD is not provided in the data, and the same is true for blanks. 

All the other options do have effect: those related to internal standard, amount of material, and isotopologues.
 
 
 
### Users having VIB results

As shown in the example 'toy2' [here](examples/toy2/), give the names of the sheets that are present in your excel file coherently. 
 
Our pipeline performs, by default:
- the subtraction of the means of the blanks across all metabolites' abundance for each sample.
- seting to NaN the values of abundance that are under the limit of detection (LOD).
- excluding metabolites whose abundance values across all samples are under LOD (excluded then from all tables by default).
- stomping fractions values to be comprised between 0 and 1 (some negative and some superior to 1 values can occur after )

You can modify all those options depending on your needs, they appear as 'optional arguments' in the help menu. 


### Users having generic data

We have created this option for those formats that are not the other two scenarios, so your data is expectd to be in the form of a .xlsx file with sheets similar as in VIB results:
- sheets corresponding to isotopologue Proportions (when available) and isotopologue Absolute values must have isotopologues as columns and samples as rows.
- sheets corresponding to abundance and mean enrichment  must have metabolites as columns and samples as rows.

As in example [toy3](examples/toy3) if you only have isotopologue Absolute values, but not the other tables: put them as a single named sheet in your .xlsx file, and we automatically generate all the other types of tables for you ! 




## The "Metadata"

Here the first lines of the required metadata table, which must be a .csv (comma delimited) file : 

| name_to_plot   | timepoint | condition | timenum | short_comp  |  original_name |
|-------------------|-----------|-----------|-------|------------|--------------- |
| Control\_cell\_T0-1 | T0        | Control   | 0     | cell       | MCF001089_TD01 |
| Control\_cell\_T0-2 | T0        | Control   | 0     | cell       | MCF001089_TD02 |
| Control\_cell\_T0-3 | T0        | Control   | 0     | cell       |  MCF001089_TD03|

You can create it with any spreadsheet program such as Excel or Google Sheets or LibreOfice Calc. At the moment of saving your file you specify that the delimiter must be a comma, see https://support.microsoft.com/en-us/office/save-a-workbook-to-text-format-txt-or-csv-3e9a9d6c-70da-4255-aa28-fcacf1f081e6. 

Column names in metadata must be exactly: 
 - original\_name
 - name_to_plot
 - timepoint
 - timenum
 - condition
 - short\_comp

 
The column 'original\_name' must have the names of the samples **as given in your data**. 
  
 
 The column 'name_to_plot' must have the names as you want them to be (or set identical to original\_name if you prefer). To set  names that are meaningful is a better choice, as we will take them for all the results.
 
 
 The column 'timenum' must contain only the numberic part of the timepoint, for example 2,0, 10, 100  (this means, without letters ("T", "t", "s", "h" etc) nor any other symbol). Make sure these time numbers are in the same units (but do not write  the units here!).
  

The column 'short\_comp' is an abbreviation, coined by you, for the compartments. This will be used for the results' files names: the longer the compartments names are, the longer the output files' names! Please pick short and clear abbreviations to fill this column.
