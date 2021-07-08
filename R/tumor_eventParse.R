#' Title Parse TCGA clinical tumor event data from XML files
#'
#'
#' @param filenames A vector contain the names of TCGA clinical XML files directorys and names (tcga_files_path/xml_file_names), which can easily retrieve by `list.files(path = ".", pattern = "*.xml",recursive = TRUE)`. The working directorys is the key, please check the names (which contain tcga_files_path/XML_file_names) is correct.
#' @param simplify When simplify = TRUE, it will return the patients basic ID information, new_tumor_neoplasm_event_type, new_tumor_occurrence_anatomic_site, new_neoplasm_event_type.
#'                 When simplify = FALSE, it will return the patients basic ID information, and a tibble column, new_tumor_event_list, contain all the tumor event of the patient.
#' @param dir file directory contained the XML files, default is current directory "."
#'
#' tumor_eventParse() will extract the TCGA xml files, as you frome the follow up record of XML files.
#' Some patient may have multiple record of tumor event, others have none or
#'
#' When simplify = TRUE, the returned tibble will contain column as following (bcr_patient_barcode, bcr_patient_uuid, age_at_initial_pathologic_diagnosis, tumor_tissue_site, gender, stage_event_system_version, stage_pathologic_stage, new_tumor_event_list).
#' last 3 column (new_tumor_neoplasm_event_type, new_tumor_occurrence_anatomic_site, new_neoplasm_event_type), is a list in a tibble column, due to one patient often have several new tumor metastasis event.
#'
#' When simplify = FALSE, the returned tibble will contain bcr_patient_barcode, bcr_patient_uuid, age_at_initial_pathologic_diagnosis, tumor_tissue_site, gender, stage_event_system_version, stage_pathologic_stage, new_tumor_event_list.
#' The last one new_tumor_event_list is the full tumor event data, store in the tibble column, and each element is a table
#'
#'
#' @return a tibble with tumor event data
#'
#'
#' @author Fan Zhang
#'
#' @seealso the vitalParse function for extract the vital data, also have the days to new tumor event, which can be used for PFS analysis.
#'  \code{link{vitalParse}}
tumor_eventParse <- function(filenames, dir=".", simplify = TRUE){

  if(!all(grepl(".xml",filenames))){
    stop("filenames must be a vector of TCGA XML file names")}

  if(!file.exists(file.path(dir,filenames[1]))){
    stop(paste("In", file.path(dir,filenames[1]),"\n cann't not find the xml files.\n Please input a correct directory (parameter dir) in this function, so i can find the files."))
  }

  one_case <- lapply(filenames,function(filenames){

    file <- XML::xmlParse(file.path(dir,filenames),encoding="UTF-8")

    file <- XML::xmlToList(file)

    patient <- file[[2]]
    rm(file)


    checkNA <- function(x){
      out_put <- tryCatch(x,error=function(c){NA})
      out_put <- ifelse(is.null(out_put),NA,out_put)
      return(out_put)
    }


    ###patient ID infor

    #extract patient id information, creat tibble
    patient_infor <- tibble::tibble(bcr_patient_barcode = checkNA(patient[["bcr_patient_barcode"]][["text"]]),
                                    bcr_patient_uuid = checkNA(patient[["bcr_patient_uuid"]][["text"]]),
                                    age_at_initial_pathologic_diagnosis = checkNA(patient[["age_at_initial_pathologic_diagnosis"]][["text"]]),
                                    tumor_tissue_site = checkNA(patient[["tumor_tissue_site"]][["text"]]),
                                    gender = checkNA(patient[["gender"]][["text"]]),
                                    stage_event_system_version = checkNA(patient[["stage_event"]][["system_version"]][["text"]]),
                                    stage_pathologic_stage = checkNA(patient[["stage_event"]][["pathologic_stage"]][["text"]]),
                                    stage_clinical_stage = checkNA(patient[["stage_event"]][["clinical_stage"]][["text"]])
    )
    #if without follow up data
    if(!is.null(patient[["follow_ups"]])){
      #extrect the new_tumor_events in the follow up
      #some cases has multiple follow up
      followup_all_tumor_event <- lapply(patient[["follow_ups"]],function(followup){
        #follow up infor

        #
        if(length(followup[["new_tumor_events"]])>=2){
          #extract from the second new_tumor_events data to the last one
          tumor_event_in_one_followup <- lapply(followup[["new_tumor_events"]][2:length(followup[["new_tumor_events"]])],function(one_tumor_event_infor){

            #unlist the list
            x <- unlist(one_tumor_event_infor)

            #select the xml text in the list by names
            x <- x[grep(".text$",names(x))]

            names(x) <- stringr::str_sub(names(x),1,-6)
            #return each from the second new_tumor_events data to the last one
            return(x)
          })

          ###merge the list, which contained single names vector, by its names
          one_followup_tumor_event <- lapply(tumor_event_in_one_followup,function(each_named_vector){
            each_named_vector <- as.data.frame(each_named_vector)
            each_named_vector <- tibble::rownames_to_column(each_named_vector,var = "term")
            return(each_named_vector)
          })

          #perform t() to
          one_followup_tumor_event <- Reduce(function(x,y) merge(x,y,by="term",all=TRUE),one_followup_tumor_event, accumulate = FALSE)

          one_followup_tumor_event <- tibble::as_tibble(one_followup_tumor_event)



          return(one_followup_tumor_event)

        } else {
          tumor_event_in_one_followup <- NA
          return(tumor_event_in_one_followup)
        }

      })

    } else {
      #if no follow up data
      followup_all_tumor_event <- NA
    }

    #rm NA
    followup_all_tumor_event <- followup_all_tumor_event[!is.na(followup_all_tumor_event)]

    norow <- sapply(followup_all_tumor_event,nrow)
    followup_all_tumor_event <- followup_all_tumor_event[norow!=0]

    #not null,
    if(!is.null(followup_all_tumor_event)){
      #combind multiple follow up new tumor event
      all_new_tumor_event <- Reduce(function(x,y) merge(x,y,by="term",all=TRUE),followup_all_tumor_event, accumulate = FALSE)

      #rename from 1 to colnames(all_new_tumor_event)
      #colnames(all_new_tumor_event)[2:length(colnames(all_new_tumor_event))] <- paste("new_tumor_event",1:( length(colnames(all_new_tumor_event))-1), sep = "_" )

    } else {
      #is null, return NA
      all_new_tumor_event <- NA
    }

    if(simplify == FALSE){

      ###three condition in the same time: not null, more than one row, at lease one not NA
      if(all(!is.null(all_new_tumor_event), ncol(all_new_tumor_event)!=0, any(!is.na(all_new_tumor_event)))){

       patient_infor$new_tumor_event_list <- list(all_new_tumor_event)
      } else {
        patient_infor$new_tumor_event_list <- list(NA)
      }
      return(patient_infor)
    }
    if(simplify == TRUE){
      if(all(!is.null(all_new_tumor_event), ncol(all_new_tumor_event)!=0, any(!is.na(all_new_tumor_event)))){
        all_new_tumor_event <- tibble::column_to_rownames(all_new_tumor_event,"term")

        all_new_tumor_event <- tibble::as_tibble(t(all_new_tumor_event))

          patient_infor$days_to_new_tumor_event_after_initial_treatment <- checkNA(list(all_new_tumor_event$days_to_new_tumor_event_after_initial_treatment))

          patient_infor$new_tumor_neoplasm_event_type <- checkNA(list(all_new_tumor_event$new_neoplasm_event_type))

          patient_infor$new_tumor_occurrence_anatomic_site <- checkNA(list(all_new_tumor_event$new_neoplasm_event_occurrence_anatomic_site))

          patient_infor$new_neoplasm_event_type <- checkNA(list(all_new_tumor_event$new_neoplasm_event_type))
      } else {
        patient_infor$days_to_new_tumor_event_after_initial_treatment <- list(NA)
        patient_infor$new_tumor_neoplasm_event_type <- list(NA)
        patient_infor$new_tumor_occurrence_anatomic_site <- list(NA)
        patient_infor$new_neoplasm_event_type <- list(NA)
      }
      return(patient_infor)
    }
  })
  do.call(rbind, one_case)
}
