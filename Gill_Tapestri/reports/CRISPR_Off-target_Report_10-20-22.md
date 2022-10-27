---
title: "CRISPR Off-Target Sites for CD45"
author: "Joe Jessee"
date: "2022-10-18"
output:
  html_document:
    toc: TRUE
    fig_caption: TRUE
  
---



## Acquiring Off-Target Sites

In this report I will show how to get the off-target sites for CD45 from "rgenome.net", annotate the target data, and convert the annotated target data for compliance with MissionBio's Mosaic Template Editor requirements.

### CRISPR RGEN Tools at regenome.net

<br>


![from rgenome.net](~/Gill_Tapestri/images/Screenshot from 2022-10-18 11-37-07.png)

<br>


|Bulge.type |crRNA                   |DNA                     |Chromosome | Position|Direction | Mismatches| Bulge.Size|
|:----------|:-----------------------|:-----------------------|:----------|--------:|:---------|----------:|----------:|
|X          |AAAATATGCAAACATCACTGNGG |cAAATtctCAAACAcCACTGTGG |chr8       |   413937|-         |          5|          0|
|X          |AAAATATGCAAACATCACTGNGG |AAAATAaGgAAAacTCcCTGAGG |chr8       |  1691247|+         |          5|          0|
|X          |AAAATATGCAAACATCACTGNGG |AAtgaATGatAACATCACTGAGG |chr8       |  2193196|+         |          5|          0|
|X          |AAAATATGCAAACATCACTGNGG |AAAATAgtCcAACAgCcCTGAGG |chr8       |  2634858|+         |          5|          0|
|X          |AAAATATGCAAACATCACTGNGG |AAtATtTGaAAtCAgCACTGAGG |chr8       |  4968618|-         |          5|          0|
|X          |AAAATATGCAAACATCACTGNGG |AggtTtTGCAAgCATCACTGTGG |chr8       |  5530766|+         |          5|          0|

