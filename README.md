## Geographic Patterns in U.S. Lung Cancer Mortality and Cigarette Smoking

![license](https://img.shields.io/badge/license-apache-yellow)

**Date repository last updated**: January 15, 2022

### Authors

* **Alaina H. Shreves**<sup>1,2</sup> [ORCID](https://orcid.org/0000-0002-0127-4391)
* **Ian D. Buller**<sup>3,4</sup> [ORCID](https://orcid.org/0000-0001-9477-8582)
* **Elizabeth Chase**<sup>5,6</sup> - [ORCID](https://orcid.org/0000-0003-0452-2976)
* **Hannah Creutzfeldt**<sup>3,7</sup>
* **Jared A. Fisher**<sup>3</sup> [ORCID](https://orcid.org/0000-0001-9203-5742)
* **Barry I. Graubard**<sup>6</sup> [ORCID](https://orcid.org/0000-0002-6787-1105)
* **Robert N. Hoover**<sup>8</sup>
* **Debra T. Silverman**<sup>3</sup>
* **Susan S. Devesa**<sup>5</sup> - *Co-Senior Author*
* **Rena R. Jones**<sup>3</sup> - *Co-Senior Author* & *Corresponding Author* - [ORCID](https://orcid.org/0000-0003-1294-1679)

1.	Department of Epidemiology, Harvard T.H. Chan School of Public Health, Harvard University, Boston, MA, 02115, USA
2.	Trans-Divisional Research Program, Division of Cancer Epidemiology and Genetics (DCEG), National Cancer Institute (NCI), National Institutes of Health (NIH), Rockville, MD, 20850, USA
3.	Occupational and Environmental Epidemiology Branch, Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, MD, 20850, USA, Rockville, MD, 20850, USA
4.  Cancer Prevention Fellowship Program, Division of Cancer Prevention, National Cancer Institute, Rockville, MD, 20850, USA
5.	Infections and Immunology Branch, DCEG, NCI, NIH, Rockville, MD, 20850, USA
6.	Department of Biostatistics, University of Michigan School of Public Health, University of Michigan, Rockville, MD, 20850, USA
7.	Fielding School of Public Health, University of California Los Angeles, Los Angeles, CA, 90095, USA
8.	Office of the Director, DCEG, NCI, NIH, Rockville, MD, 20850, USA

### Project details
Lung cancer is the leading cause of cancer death in the United States (US) and variations in lung cancer mortality and smoking behavior are evident by sex and region. We apply geospatial statistical methods to describe patterns in lung cancer mortality rates (2005-2018) in relation to patterns in cigarette smoking prevalences (1997-2003) by sex at the US county level. Our findings identify counties where lung carcinogens other than smoking may be driving lung cancer mortality and where further study is needed. 

#### Project timeframe

<table>
<colgroup>
<col width="20%" />
<col width="80%" />
</colgroup>
<thead>
<tr class="header">
<th>Time</th>
<th>Event</th>
</tr>
</thead>
<tbody>
<td><p align="center">1997-2003</p></td>
<td>NCI Model-based Small Area Estimates of Cancer-Related Measures smoking prevalences for persons aged 18+ years</td>
</tr>
<td><p align="center">2005-2018</p></td>
<td>Lung and bronchus cancer mortality rates among persons aged 20+ years from the National Vital Statistics System data from the National Center for Health Statistics</td>
</tr>
<td><p align="center">July 2020</p></td>
<td>Project Initiation</td>
</tr>
<td><p align="center">TBD</p></td>
<td>Initial manuscript submission for peer-review</td>
</tr>
<td><p align="center">TBD</p></td>
<td>Manuscript accepted in a peer-reviewed journal</td>
</tr>
<td><p align="center">TBD</p></td>
<td>Manuscript published in a peer-reviewed journal</td>
</tr>

</tbody>
<table>

### R-scripts included in this repository

This repository includes R-scripts use to calculate the geospatial techniques and render the figures found in the following peer-reviewed manuscript:

[INSERT SUGGESTED CITATION HERE]

<table>
<colgroup>
<col width="20%" />
<col width="80%" />
</colgroup>
<thead>
<tr class="header">
<th>R-Script</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<td><p align="center"><code>functions.R</code></td>
<td>Custom functions to calculate the local Lee's L statistic with correction for multiple testing</td>
</tr>
<td><p align="center"><code>preparation.R</code></td>
<td>Calculate the local Lee's L statistics for the four comparisons. Requires a data set to run (not included). See notes within. </td>
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

### Questions?

For questions about the manuscript please e-mail the corresponding author [Dr. Rena R. Jones](mailto:rena.jones@nih.gov)
