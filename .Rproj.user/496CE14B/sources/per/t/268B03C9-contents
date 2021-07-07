#' Parse TCGA clinical data patients stage information
#'
#' @param filenames A vector contain the names of TCGA clinical XML files directorys and names (tcga_files_path/xml_file_names), which can easily retrieve by `list.files(path = ".", pattern = "*.xml",recursive = TRUE)`. The working directorys is the key, please check the names (which contain tcga_files_path/XML_file_names) is correct.
#' @param simplify should return be simplified ?
#'
#'
#'
#'
#'
#' @author Fan Zhang
#'
#'
#'
#'
#'
#'
#'
#'
#'
stageParse <- function(filenames, simplify = TRUE){

  if(!all(grepl(".xml",filenames) == TRUE)){
    stop("filenames must be a vector of TCGA XML file names")}

  stage_infor <- lapply(filenames,function(filename){

        file <- XML::xmlParse(filename,encoding="UTF-8")

        file <- XML::xmlToList(file)

        patient <- file[[2]]
        rm(file)


        checkNA <- function(x){
          out_put <- tryCatch(x,error=function(c){NA})
          out_put <- ifelse(is.null(out_put),NA,out_put)
          return(out_put)
        }

        if(simplify == FALSE){
          #extract patient stage information, creat tibble
          patient_infor <- tibble::tibble(bcr_patient_barcode = checkNA(patient[["bcr_patient_barcode"]][["text"]]),
                                          bcr_patient_uuid = checkNA(patient[["bcr_patient_uuid"]][["text"]]),
                                          age_at_initial_pathologic_diagnosis = checkNA(patient[["age_at_initial_pathologic_diagnosis"]][["text"]]),
                                          tumor_tissue_site = checkNA(patient[["tumor_tissue_site"]][["text"]]),
                                          gender = checkNA(patient[["gender"]][["text"]]),
                                          race = checkNA(patient[["race_list"]][["race"]][["text"]]),
                                          stage_system_version = checkNA(patient[["stage_event"]][["system_version"]][["text"]]),
                                          stage_clinical_stage = checkNA(patient[["stage_event"]][["clinical_stage"]][["text"]]),
                                          stage_pathologic_stage = checkNA(patient[["stage_event"]][["pathologic_stage"]][["text"]]),
                                          stage_tnm_clinical_T = checkNA(patient[["stage_event"]][["tnm_categories"]][["clinical_categories"]][["clinical_T"]][["text"]]),
                                          stage_tnm_clinical_N = checkNA(patient[["stage_event"]][["tnm_categories"]][["clinical_categories"]][["clinical_N"]][["text"]]),
                                          stage_tnm_clinical_M = checkNA(patient[["stage_event"]][["tnm_categories"]][["clinical_categories"]][["clinical_M"]][["text"]]),
                                          stage_tnm_pathologic_T = checkNA(patient[["stage_event"]][["tnm_categories"]][["pathologic_categories"]][["pathologic_T"]][["text"]]),
                                          stage_tnm_pathologic_M = checkNA(patient[["stage_event"]][["tnm_categories"]][["pathologic_categories"]][["pathologic_N"]][["text"]]),
                                          stage_tnm_pathologic_N = checkNA(patient[["stage_event"]][["tnm_categories"]][["pathologic_categories"]][["pathologic_M"]][["text"]]),
                                          stage_psa_value = checkNA(patient[["stage_event"]][["psa"]][["psa_value"]][["text"]]),
                                          stage_days_to_psa = checkNA(patient[["stage_event"]][["psa"]][["days_to_psa"]][["text"]]),
                                          stage_gleason_gleason_score = checkNA(patient[["stage_event"]][["gleason_grading"]][["gleason_score"]][["text"]]),
                                          stage_gleason_primary_pattern = checkNA(patient[["stage_event"]][["gleason_grading"]][["primary_pattern"]][["text"]]),
                                          stage_gleason_secondary_pattern = checkNA(patient[["stage_event"]][["gleason_grading"]][["secondary_pattern"]][["text"]]),
                                          stage_gleason_tertiary_pattern = checkNA(patient[["stage_event"]][["gleason_grading"]][["tertiary_pattern"]][["text"]]),
                                          stage_ann_arbor_b_symptoms = checkNA(patient[["stage_event"]][["ann_arbor"]][["b_symptoms"]][["text"]]),
                                          stage_ann_arbor_extranodal_involvement = checkNA(patient[["stage_event"]][["ann_arbor"]][["extranodal_involvement"]][["text"]]),
                                          stage_serum_markers = checkNA(patient[["stage_event"]][["serum_markers"]][["text"]]),
                                          stage_igcccg_stage = checkNA(patient[["stage_event"]][["igcccg_stage"]][["text"]]),
                                          stage_masaoka_stage = checkNA(patient[["stage_event"]][["masaoka_stage"]][["text"]])
          )
          return(patient_infor)
        }


        if(simplify == TRUE){
          patient_infor <- tibble::tibble(bcr_patient_barcode = checkNA(patient[["bcr_patient_barcode"]][["text"]]),
                                          bcr_patient_uuid = checkNA(patient[["bcr_patient_uuid"]][["text"]]),
                                          age_at_initial_pathologic_diagnosis = checkNA(patient[["age_at_initial_pathologic_diagnosis"]][["text"]]),
                                          tumor_tissue_site = checkNA(patient[["tumor_tissue_site"]][["text"]]),
                                          gender = checkNA(patient[["gender"]][["text"]]),
                                          race = checkNA(patient[["race_list"]][["race"]][["text"]]),
                                          stage_system_version = checkNA(patient[["stage_event"]][["system_version"]][["text"]]),
                                          stage_pathologic_stage = checkNA(patient[["stage_event"]][["pathologic_stage"]][["text"]]),
                                          stage_tnm_pathologic_T = checkNA(patient[["stage_event"]][["tnm_categories"]][["pathologic_categories"]][["pathologic_T"]][["text"]]),
                                          stage_tnm_pathologic_M = checkNA(patient[["stage_event"]][["tnm_categories"]][["pathologic_categories"]][["pathologic_N"]][["text"]]),
                                          stage_tnm_pathologic_N = checkNA(patient[["stage_event"]][["tnm_categories"]][["pathologic_categories"]][["pathologic_M"]][["text"]])
                           )
          return(patient_infor)
        }
  })
  stage_infor <- do.call(rbind, stage_infor)
}
