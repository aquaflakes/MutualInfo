---
title: "Mutual Infomation Calculation for SELEX"
output: 
  html_document: 
    keep_md: yes
---



## Introduction

To footprint the binding of transcription factors on DNA, the traditional motif mapping analysis requires previously curated motifs, and is limited to profile only 1 motif at a time. To overcome such issues, this package provides the implementation of the mutual-infomation based footprint analysis for transcription factors. 

Check also <https://www.nature.com/articles/s41586-018-0549-5> for details on how to use.


## Prerequisit
To use the functions in this package. Please note that the package **fjComm** needs to be installed before hand. There is an `install_github` function to install R packages hosted on GitHub in the **devtools** package. It requests developer’s name.:


```r
install_github("aquaflakes/fjComm")
```

## Input data requirement
The input data needs to be stored in a dataframe with only 1 column, each row corresponds to one sequence from the illumina sequencing. The data frame look like follows:

```
##                                      .
## 1 ACTACTATCTATTCTATCATTACTACTACATATCAT
## 2 ACTACTATCTATTCTATCATTACTACTACATATCAT
## 3 ACTACTATCTATTCTATCATTACTACTACATATCAT
## 4 ACTACTATCTATTCTATCATTACTACTACATATCAT
## 5 ACTACTATCTATTCTATCATTACTACTACATATCAT
```


## Calculate Mutual Infomation

You can also embed plots, for example:

![](Readme_files/figure-html/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
