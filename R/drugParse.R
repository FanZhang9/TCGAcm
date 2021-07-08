#' Parse TCGA clinical data drug information
#'
#'
#'
#' @description `drugParse()` extract the drug record information from TCGA clinical XML files.
#'
#' For each value in `filenames` `drugParse()` will extract the drug record information from clinical XML files downloaded from the TCGA database.
#'
#'
#'
#'
#'
#' @param filenames A vector contain the names of TCGA clinical XML files directorys and names (tcga_files_path/xml_file_names), which can easily retrieve by `list.files(path = ".", pattern = "*.xml",recursive = TRUE)`. The working directorys is the key, please check the names (which contain tcga_files_path/XML_file_names) is correct.
#'
#' @param simplify logical, when simplify = TRUE return a tibble with patient ID information (bcr_patient_barcode, bcr_patient_uuid, tumor_tissue_site, gender), and a column contain
#'                 drug names. (it's a list column, thank to the tibble freature, and it contain the patient drug names has undergone).
#'                 When it's simplify = FALSE, it will return return a tibble with patient ID information (bcr_patient_barcode, bcr_patient_uuid, tumor_tissue_site, gender), and a column contain
#'                 full drug record.
#'
#' @author Fan Zhang
#'
#' @return When simplify = TRUE, return a tibble with patients id and drug names list column. When simplify = FALSE, return a tibble contain full drugs record information.
drugParse <- function(filenames,simplify = TRUE){

  if(!all(grepl(".xml",filenames) == TRUE)){
    stop("filenames must be a vector of TCGA XML file names")}

  #parse XML_filenames extract drug information
  drug_list <- lapply(filenames,function(x){

    file <- XML::xmlParse(x,encoding="UTF-8")

    file <- XML::xmlToList(file)

    patient <- file[[2]]
    rm(file)


    checkNA <- function(x){
      out_put <- tryCatch(x,error=function(c){NA})
      out_put <- ifelse(is.null(out_put),NA,out_put)
      return(out_put)
    }

    #extract patient id information, creat tibble
    patient_infor <- tibble::tibble(bcr_patient_barcode = checkNA(patient[["bcr_patient_barcode"]][["text"]]),
                                    bcr_patient_uuid = checkNA(patient[["bcr_patient_uuid"]][["text"]]),
                                    age_at_initial_pathologic_diagnosis = checkNA(patient[["age_at_initial_pathologic_diagnosis"]][["text"]]),
                                    tumor_tissue_site = checkNA(patient[["tumor_tissue_site"]][["text"]]),
                                    gender = checkNA(patient[["gender"]][["text"]]),
                                    stage_event_system_version = checkNA(patient[["stage_event"]][["system_version"]][["text"]]),
                                    stage_pathologic_stage = checkNA(patient[["stage_event"]][["pathologic_stage"]][["text"]])
                                    )


    drug <- patient[["drugs"]]

    rm(patient)

    if(!is.null(drug)){

            #extract the xml text node
            drug <- lapply(drug,function(drug){

              #unlist the list
              x <- unlist(drug)

              #select the xml text in the list by names
              x <- x[grep(".text$",names(x))]

              names(x) <- stringr::str_sub(names(x),1,-6)

              return(x)
            })

            #out put all the drug information into data.frame

            ###merge the list, which contained single names vector, by its names
            drug <- lapply(drug,function(drug_record){
              drug_record <- as.data.frame(drug_record)
              drug_record <- tibble::rownames_to_column(drug_record,var = "term")
              return(drug_record)
            })

            #perform t() to
            drug <- suppressWarnings(Reduce(function(x,y) merge(x,y,by="term",all=TRUE),drug, accumulate = FALSE))

            drug <- as.data.frame(t(drug))

            colnames(drug) <- drug[1,]

            drug <- drug[-1,]

            drug_name <- drug$drug_name

            #drug list write into tibble
            #if drug is null write NA
            if(length(drug_name)==0){
              patient_infor$drug_list <- list(drug=NA)
            } else {
                if(simplify == TRUE){
                patient_infor$drug_list <- list(drug=drug$drug_name)
                }
                if(simplify == FALSE){
                patient_infor$drug_list <- list(drug=drug)
                }
            }

    }  else {
      patient_infor$drug_list <- list(NA)
     }

    return(patient_infor)

  })

    drug_list <- do.call(rbind, drug_list)
    return(drug_list)

}



