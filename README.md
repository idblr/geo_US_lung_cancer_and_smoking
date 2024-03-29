Geographic Patterns in U.S. Lung Cancer Mortality and Cigarette Smoking <img src="hex/hex.png" width="120" align="right"/>
===================================================

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
![GitHub last commit](https://img.shields.io/github/last-commit/idblr/geo_US_lung_cancer_and_smoking)

**Date repository last updated**: June 10, 2023

### Authors

* **Alaina H. Shreves**<sup>1,2</sup> - [ORCID](https://orcid.org/0000-0002-0127-4391)
* **Ian D. Buller**<sup>3,4</sup> - [ORCID](https://orcid.org/0000-0001-9477-8582)
* **Elizabeth Chase**<sup>5,6</sup> - [ORCID](https://orcid.org/0000-0003-0452-2976)
* **Hannah Creutzfeldt**<sup>3,7</sup> - [ORCID](https://orcid.org/0000-0001-7828-8735)
* **Jared A. Fisher**<sup>3</sup> - [ORCID](https://orcid.org/0000-0001-9203-5742)
* **Barry I. Graubard**<sup>6</sup> - [ORCID](https://orcid.org/0000-0002-6787-1105)
* **Robert N. Hoover**<sup>8</sup> - [ORCID](https://orcid.org/0000-0003-4329-4371)
* **Debra T. Silverman**<sup>3</sup> - [ORCID](https://orcid.org/0000-0001-8894-0301) 
* **Susan S. Devesa**<sup>5</sup> - *Co-Senior Author* - [ORCID](https://orcid.org/0000-0003-3235-4148)
* **Rena R. Jones**<sup>3</sup> - *Co-Senior Author* & *Corresponding Author* - [ORCID](https://orcid.org/0000-0003-1294-1679)

1.	Department of Epidemiology, Harvard T.H. Chan School of Public Health, Harvard University, Boston, MA, 02115, USA
2.	Trans-Divisional Research Program, Division of Cancer Epidemiology and Genetics (DCEG), National Cancer Institute (NCI), National Institutes of Health (NIH), Rockville, MD, 20850, USA
3.	Occupational and Environmental Epidemiology Branch, DCEG, NCI, Rockville, MD, 20850, USA
4.  Cancer Prevention Fellowship Program, Division of Cancer Prevention, NCI, Rockville, MD, 20850, USA
5.	Infections and Immunology Branch, DCEG, NCI, NIH, Rockville, MD, 20850, USA
6.	Department of Biostatistics, University of Michigan School of Public Health, University of Michigan, Rockville, MD, 20850, USA
7.	Fielding School of Public Health, University of California Los Angeles, Los Angeles, CA, 90095, USA
8.	Office of the Director, DCEG, NCI, NIH, Rockville, MD, 20850, USA

### Project Details
Lung cancer is the leading cause of cancer death in the United States (US) and variations in lung cancer mortality and smoking behavior are evident by sex and region. We apply geospatial statistical methods to describe patterns in lung cancer mortality rates (2005-2018) in relation to patterns in cigarette smoking prevalences (1997-2003) by sex at the US county level. Our findings identify counties where lung carcinogens other than smoking may be driving lung cancer mortality and where further study is needed. 

#### Project Timeframe

<table>
<colgroup>
<col width="20%"/>
<col width="80%"/>
</colgroup>
<thead>
<tr class="header">
<th>Time</th>
<th>Event</th>
</tr>
</thead>
<tbody>
<td><p align="center">1997-2003</p></td>
<td>NCI Model-based Small Area Estimates of Cancer-Related Measures smoking prevalences for persons aged 18+ years (see data availability section below)</td>
</tr>
<td><p align="center">2005-2018</p></td>
<td>Lung and bronchus cancer mortality rates among persons aged 20+ years from the National Vital Statistics System data from the National Center for Health Statistics (see data availability section below)</td>
</tr>
<td><p align="center">July 2020</p></td>
<td>Project Initiation</td>
</tr>
<td><p align="center">March 2022</p></td>
<td>Initial manuscript submission to <a href="https://cebp.aacrjournals.org/">Cancer Epidemiology, Biomarkers & Prevention</a> for peer-review</td>
</tr>
<td><p align="center">November 2022</p></td>
<td>Manuscript accepted by <a href="https://cebp.aacrjournals.org/">Cancer Epidemiology, Biomarkers & Prevention</a></td>
</tr>
<td><p align="center">February 2023</p></td>
<td>Manuscript published in <a href="https://doi.org/10.1158/1055-9965.EPI-22-0253">Cancer Epidemiology, Biomarkers & Prevention</a></td>
</tr>
<td><p align="center">June 2023</p></td>
<td>Update to the False Discovery Rate <a href="https://doi.org/10.1111/j.2517-6161.1995.tb02031.x">(Benjamini & Hochberg, 1995)</a> calculation for multiple testing correction that now orders the p-values in ascending order instead of in descending order.</td>
</tr>
</tbody>
<table>

### R Scripts Included In This Repository

This repository includes R scripts used to calculate the Lee's L statistic and render the geographic visualizations found in the following peer-reviewed manuscript:

Shreves AH, Buller ID, Chase E, Creutzfeld H, Fisher JA, Graubard BI, Hoover RN, Silverman DT, Devesa SS, Jones RR. (2023) Geographic Patterns in U.S. Lung Cancer Mortality and Cigarette Smoking. _Cancer Epidemiology, Biomarkers & Prevention_, 32(2):193-201. DOI:[10.1158/1055-9965.EPI-22-0253](https://doi.org/10.1158/1055-9965.EPI-22-0253) PMID:[36413442](https://pubmed.ncbi.nlm.nih.gov/36413442/).

<table>
<colgroup>
<col width="20%"/>
<col width="80%"/>
</colgroup>
<thead>
<tr class="header">
<th>R Script</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<td><p align="center"><code>functions.R</code></td>
<td>Custom functions to calculate the local Lee's L statistic with correction for multiple testing</td>
</tr>
<td><p align="center"><code>preparation.R</code></td>
<td>Calculate the local Lee's L statistics for the four comparisons. Requires a data set to run (not included; see notes within). </td>
</tr>
<td><p align="center"><code>figure1.R</code></p></td>
<td>Generate Figure 1</td>
</tr>
<td><p align="center"><code>figure2.R</code></p></td>
<td>Generate Figure 2</td>
</tr>
<td><p align="center"><code>supplemental1.R</code></p></td>
<td>Generate Supplemental Figure 1</td>
</tr>
<td><p align="center"><code>supplemental2.R</code></p></td>
<td>Generate Supplemental Figure 2</td>
</tr>
</tbody>
<table>

The repository also includes the code to create the project hexagon sticker.

### Getting Started

* Step 1: You must download the data (see Data Availability section)
* Step 2: Save the data set to the data directory in this repository. Currently specified as a CSV file, but modify the path on Line 58 of the `preparation.R` file based on data location and file name
* Step 3: Run R scripts for figures. The `preparation.R` file will source the `functions.R` file.

### Data Availability

County-level U.S. lung cancer mortality rates and smoking prevalences are downloadable from [Model-based Small Area Estimates of Cancer-Related Measures](https://sae.cancer.gov/nhis-brfss/) from the [Surveillance Research Program](https://surveillance.cancer.gov/) within the [Division of Cancer Control and Population Sciences](https://cancercontrol.cancer.gov/) of the [National Cancer Institute](https://www.cancer.gov/) and the [National Vital Statistics System](https://www.cdc.gov/nchs/nvss/index.htm) from the [National Center for Health Statistics](https://www.cdc.gov/nchs/index.htm) of the [Centers for Disease Control and Prevention](https://www.cdc.gov/).

### Questions?

For questions about the manuscript please e-mail the corresponding author [Dr. Rena R. Jones](mailto:rena.jones@nih.gov).
