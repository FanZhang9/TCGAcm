
#' Parse TCGA clinical follow up data
#'
#' @param filenames A vector contain the names of TCGA clinical XML files directorys and names (tcga_files_path/xml_file_names), which can easily retrieve by `list.files(path = ".", pattern = "*.xml",recursive = TRUE)`. The working directorys is the key, please check the names (which contain tcga_files_path/XML_file_names) is correct.
#' @param dir file directory contained the XML files, default is current directory "."
#' @param simplify should the output be simplified
#'
#'
#' @author Fan Zhang <fanzhang95@outlook.com>
#'
#'
follow_upParse <- function(filenames, dir=".", simplify = TRUE){

  if(!all(grepl(".xml",filenames))){
    stop("filenames must be a vector of TCGA XML file names")}

  if(!file.exists(file.path(dir,filenames[1]))){
    stop(paste("In", file.path(dir,filenames[1]),"\n cann't not find the xml files.\n Please input a correct directory (parameter dir) in this function, so i can find the files."))
  }


  one_case <- lapply(filenames,function(filenames){

      xfile <- xml2::read_xml(file.path(dir,filenames))

      patient <- xml2::as_list(xfile)[[1]][[2]]

      checkNA <- function(x){
        out_put <- tryCatch(x,error=function(c){NA})
        out_put <- ifelse(is.null(out_put),NA,out_put)
        return(out_put)
      }

      checkNAnull <- function(x){
        out_put <- tryCatch(x,error=function(c){NA})
        out_put <- ifelse(is.null(out_put[[1]]),NA,out_put)
        return(out_put)
      }

      #extract patient id information, creat tibble
      patient_infor <- tibble::tibble(bcr_patient_barcode = checkNA(patient[["bcr_patient_barcode"]][[1]]),
                                      bcr_patient_uuid = checkNA(patient[["bcr_patient_uuid"]][[1]]),
                                      age_at_initial_pathologic_diagnosis = checkNA(patient[["age_at_initial_pathologic_diagnosis"]][[1]]),
                                      tumor_tissue_site = checkNA(patient[["tumor_tissue_site"]][[1]]),
                                      gender = checkNA(patient[["gender"]][[1]]),
                                      stage_event_system_version = checkNA(patient[["stage_event"]][["system_version"]][[1]]),
                                      stage_pathologic_stage = checkNA(patient[["stage_event"]][["pathologic_stage"]][[1]]),
                                      stage_clinical_stage = checkNA(patient[["stage_event"]][["clinical_stage"]][[1]]),
                                      vital_status = checkNA(patient[["vital_status"]][[1]]),
                                      days_to_last_followup = checkNA(patient[["days_to_last_followup"]][[1]])
      )



      followup <- patient[["follow_ups"]]
      rm(patient)

      if(!length(followup)==0){

      followup <- lapply(followup,unlist)

      ###merge the list, which contained single names vector, by its names
      followup <- lapply(followup,function(each_named_vector){
        each_named_vector <- as.data.frame(each_named_vector)
        each_named_vector <- tibble::rownames_to_column(each_named_vector,var = "term")
        return(each_named_vector)
      })

      #perform t() to
      followup <- Reduce(function(x,y) merge(x,y,by="term",all=TRUE),followup, accumulate = FALSE)

      followup <- tibble::column_to_rownames(followup,var="term")

      followup <- tibble::as_tibble(t(followup))

      } else {
        followup <- NA
      }

      if(simplify==FALSE){
        patient_infor$followup_data <- list(followup)
        return(patient_infor)
      }
      if(simplify==TRUE){

        #change the colnames so the names can be matched
        if(!is.na(followup)){
          colnam <- stringr::str_split(colnames(followup),"\\.")
          colnames(followup)  <- lapply(colnam,tail,n=1)

        }

        patient_infor$followupInfor_days_to_last_followup <- checkNAnull(list(followup$days_to_last_followup))

        patient_infor$followupInfor_days_to_death <- checkNAnull(list(followup$days_to_death))

        patient_infor$followupInfor_vital_status <- checkNAnull(list(followup$vital_status))

        patient_infor$followupInfor_days_to_new_tumor_event_after_initial_treatment <- checkNAnull(list(followup$days_to_new_tumor_event_after_initial_treatment))

        patient_infor$followupInfor_new_tumor_neoplasm_event_type <- checkNAnull(list(followup$new_neoplasm_event_type))

        patient_infor$followupInfor_new_tumor_occurrence_anatomic_site <- checkNAnull(list(followup$new_neoplasm_event_occurrence_anatomic_site))

        patient_infor$followupInfor_new_tumor_event_after_initial_treatment <- checkNAnull(list(followup$new_tumor_event_after_initial_treatment))



        return(patient_infor)
      }

     })
  do.call(rbind, one_case)

}










