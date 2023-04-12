# Tracegroomer

Processes **Trace**r metabolomics given file(s), producing the .csv files which are ready for DIMet analysis.


## Requirements

The requirements are the same as in DIMet https://github.com/cbib/DIMet
First clone DIMet and create dimet conda environment (see https://github.com/cbib/DIMet), then clone Tracegroomer. Both work with the same conda environment
```
conda activate dimet
```

## Tracegroomer.tidy

Its execution takes only few seconds! Below we explain how to proceed.

### Input files


Three types of Tracer data are accepted by our `Tracegroomer.tidy` module:

1. IsoCor results (.xlsx measurments file).
2. Results provided by the VIB Metabolomics Expertise Center (El-Maven results are shaped by VIB MEC team into a multi-sheet .xlsx file).  
3. A 'generic' .xlsx measurments file.


Your data is expected to correspond to one of the three types above. We provide you "toy groom examples" of each one of them  [here](groomexamples/). Pick the example most suited to your data:

1. 'toyp1' : IsoCor output (tsv file)
2. 'toyp2' : VIB MEC xlsx file
3. 'toyp3' : a generic type of xlsx file

Each of these example folders contains:

 - Compulsory files:
     * the tracer metabolome file (tsv, csv -tab delimited-, or xlsx)
     * your configuration file (extension .yml), see below.

 - 'data/' folder, with the **metadata** file. This file provides the samples description. This is compulsory. See "metadata" section below. 
 
 - Facultative files:
   * the amount of material by sample (csv file)
   * a file with metabolites to exclude (csv file)
    
    
For example, in [toyp2 example](groomexamples/toyp2/), the tracer metabolome file is a xlsx file.     
It depends on the platform/software used before us.

Regarding the **metadata**, we explain it in detail in the section [Metadata](#the-metadata).


Regarding the .yml file, you can use as template the ones we provide in the groomexamples, such  [this config.yml template](groomexamples/toyp1/config-1-groom.yml). To double-check your modifications there exist online editors, such as https://yamlchecker.com/, just copy-paste and edit!

Note that this pipeline does not correct for naturally ocurring isotopologues. Your data must be already processed by another software that performs such correction.


### Execute `Tracegroomer.tidy` 

The groomexamples serve also to demonstrate how fast this module can be. To run all the examples at once, copy-paste the 'groomexamples' folder in your $HOME, then run:

```
python3 -m Tracegroomer.tidy --targetedMetabo_path ~/groomexamples/toyp1/TRACER_IsoCor_out_example.tsv --type_of_file IsoCor_out_tsv ~/groomexamples/toyp1/config-1-groom.yml


python3 -m Tracegroomer.tidy --targetedMetabo_path ~/groomexamples/toyp2/TRACER_metabo_toy2.xlsx --type_of_file VIBMEC_xlsx --amountMaterial_path ~/groomexamples/toyp2/nbcells-or-amountOfMaterial.csv ~/groomexamples/toyp2/config-2-groom.yml


python3 -m Tracegroomer.tidy --targetedMetabo_path ~/groomexamples/toyp3/TRACER_generic_toy3.xlsx --type_of_file generic_xlsx ~/groomexamples/toyp3/config-3-groom.yml
```



## The output

The four output tables (three if absolute isotopologues are not provided) are saved in the  folder that  you specified in the config .yml file ("out_path" field). 

The examples we provide here use the `data/` folder (inside each given "toyp*" directory) as the output location. In this way we simply copy the entire `data/` folder to the tracer metabolomics project we want to run with DIMet !.

The format of these output files is tab-delimited .csv.


--------------------
## Details about Tracegroomer.tidy

### Advanced Tracegroomer.tidy options

We provide advanced options for `Tracegroomer.tidy` module, check the help:
```
python -m Tracegroomer.tidy --help
```
they appear as 'optional arguments' in the help menu.


You can:

- normalize by the amount of material (number of cells, tissue weight): setting the path to the file in your **.yml** configuration. The file must be like [this csv file](groomexamples/toyp2/nbcells-or-amountOfMaterial.csv), and the first column must contain the same names as in metadata 'former\_name'.
- normalize by an internal standard (present in your data) at choice: using the advanced option `--use_internal_standard`.

However we have some indications that can slightly differ for [users having VIB results as input](#users-having-vib-results), [users having IsoCor results](#users-having-isocor-results) or users having ['generic' type of data](#users-having-generic-data).


 
### Users having IsoCor results

A typical IsoCor results table is described in: https://isocor.readthedocs.io/en/latest/tutorials.html#output-files
 
 Our pipeline transforms its columns into tables, so here the 'Isocor column : DIMet table' correspondences:
 
- corrected_area : isotopologuesAbsolute  
- isotopologue_fraction : isotopologuesProportions
- mean\_enrichment :  meanEnrich\_or_fracContrib
- (nothing)   : Abundance 

Abundance table is the sum (per metabolite) of IsotopologuesAbsolute, we also perform it automatically, as this column is not present in the input data.     
    
Please stick to the example [toyp1](groomexamples/toyp1) for the names of the tables in the **.yml** file for isocorOutput
    
        
Options regarding to detection limit (LOD) and blanks will not have any effect on the IsoCor type of data: LOD is not provided in the data, and the same is true for blanks. 

All the other options do have effect: those related to internal standard, amount of material, and isotopologues.
 
 
 
### Users having VIB results

As shown in the example 'toyp2' [here](groomexamples/toyp2/), give the names of the sheets that are present in your excel file coherently. 
 
Our pipeline performs, by default:
- the subtraction of the means of the blanks across all metabolites' abundance for each sample.
- seting to NaN the values of abundance that are under the limit of detection (LOD).
- excluding metabolites whose abundance values across all samples are under LOD (excluded then from all tables by default).
- stomping fractions values to be comprised between 0 and 1 (some negative and some superior to 1 values can occur after correction of naturally occurring isotopologues when using software dedicated to such corrections)

You can modify all those options depending on your needs, they appear as 'optional arguments' in the help menu. 


### Users having generic data

We have created this option for those formats that are not the other two scenarios, so your data is expectd to be in the form of a .xlsx file with sheets similar as in the provided 'toyp3'  :
- sheets corresponding to isotopologue Proportions (when available) and isotopologue Absolute values must have isotopologues as columns and samples as rows.
- sheets corresponding to abundance and mean enrichment  must have metabolites as columns and samples as rows.

As in example [toy3](groomexamples/toy3) if you only have isotopologue Absolute values, but not the other tables: put them as a single named sheet in your .xlsx file, and we automatically generate all the other types of tables for you ! 




## The "Metadata"

Here the first lines of the required metadata table, which must be a .csv (comma delimited) file : 

| name_to_plot   | condition |  timepoint | timenum | short_comp  |  original_name |
|-------------------|-----------|-----------|-------|------------|--------------- |
| Control\_cell\_T0-1 | Control  | T0h   | 0     | cell       | MCF001089_TD01 |
| Control\_cell\_T0-2 | Control  | T0h   | 0     | cell       | MCF001089_TD02 |
| Control\_cell\_T0-3 | Control  | T0h   | 0     | cell       |  MCF001089_TD03|

You can create it with any spreadsheet program such as Excel or Google Sheets or LibreOfice Calc. At the moment of saving your file you specify that the delimiter must be a comma, see https://support.microsoft.com/en-us/office/save-a-workbook-to-text-format-txt-or-csv-3e9a9d6c-70da-4255-aa28-fcacf1f081e6. 

Column names in metadata must be exactly: 
 - original\_name
 - name\_to\_plot
 - timepoint
 - timenum
 - condition
 - short\_comp

 
The column 'original\_name' must have the names of the samples **as given in your data**. 
  
 
 The column 'name\_to\_plot' must have the names as you want them to be (or set identical to original\_name if you prefer). To set  names that are meaningful is a better choice.
 
 
 The column 'timenum' must contain only the numberic part of the timepoint, for example 2,0, 10, 100  (this means, without letters ("T", "t", "s", "h" etc) nor any other symbol). Make sure these time numbers are in the same units (but do not write  the units here!).
  

The column 'short\_comp' is an abbreviation, coined by you, for the compartments. This will be used for the results' files names: the longer the compartments names are, the longer the output files' names! Please pick short and clear abbreviations to fill this column.
