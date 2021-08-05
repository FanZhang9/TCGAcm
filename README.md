
# TCGAcm is An R package for TCGA clinical data mining.

<!-- badges: start -->
<!-- badges: end -->
Since the TCGA clinical data were collected in a hug unfriendly and unreadable xml files, it's become a hug problem for us to properly dig deep into the TCGA potency. The TCGA XML files contained all kind of information, by it cannot be properly used due to their tree like structures. One person can have multiple drugs records in the same time, and of course got multiple follow up data, not surprised. So this package not only provide the tools for integrating the clinical data, but also write down some through about the cancer and data mining.

## Installation

Install from github.
```
library(devtools)
install_github("FanZhang9/TCGAmc")
```


## 1. Vital data including OS and PFS

The function `vitalParse()` can integrate TCGA clinical data from TCGA XML files, including the old vital data and updated latest follow up vital data. It can also retrieve minimum days to new tumor event after initial treatment.

First, you may download the XML data from TCGA data portal.

### input prepare
The only thing is to prepare the `filenames` as a input. It can contain a path `"data/nationwidechildrens.org_clinical.TCGA-A2-A3XY.xml"` or just a `"nationwidechildrens.org_clinical.TCGA-A2-A3XY.xml"` file name, but you must put it in your working directory. 

```{r echo=TRUE}
filenames <- list.files(path = ".", pattern = "*.xml",recursive = TRUE)
filenames
```

### function
The out put including (bcr_patient_barcode, tumor_tissue_site, vital_status, days_to_death, days_to_last_followup, new_tumor_event_after_initial_treatment, days_to_new_tumor_event_after_initial_treatment, days_tumor_free_till_last_followup).

1. OS
This three (vital_status, days_to_death, days_to_last_followup) can use as COX survival analysis.

2. PFS
This three (new_tumor_event_after_initial_treatment, days_to_new_tumor_event_after_initial_treatment, days_tumor_free_till_last_followup) can use as COX tumor free survival analysis.


```{r echo=TRUE, warning=FALSE}
library(TCGAcm)
vital_data <- vitalParse(filenames)
```


