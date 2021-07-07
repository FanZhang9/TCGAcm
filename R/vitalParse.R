# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#'  Parse TCGA clinical vital data from XML to OS and PFS for COX analysis
#'
#' @description `vitalParse()` integate TCGA clinical data from XML including the old vital data and updated latest follow up vital data. And the minimum days to new tumor event after initial treatment.
#'
#' For each value in `list` (in fact, it's a vector) `vitalParse()` extract the vital data from clinical XML files downloaded from the TCGA database
#' including vital_status and days to last follow up (OS), the latest follow_up vital data, which is updataed vital data and it can be easily omited. Meanwhile the days to new tumor event (PFS) is the earliest days to tumor event, after_initial_treatment.
#' So we can use this OS, PFS to perform COX analysis. When simplify=TRUE,if the patient without new tumor event, it's pregression free time will be last follow up time.
#'
#'
#' @param list A vector names of TCGA clinical XML files directorys and names (tcga_files_path/xml_file_names), which can easily retrieve by `list.files(path = ".", pattern = "*.xml",recursive = TRUE)`. The working directorys is the key, please check the names (which contain tcga_files_path/XML_file_names) is correct.
#'
#' @param simplify logical, when it's FLASE return TCGA vital information (including bcr_patient_barcode, vital_status, days_to_death, followupInfor_days_to_death,
#' days_to_last_followup, followupInfor_vital_status, followupInfor_days_to_last_followup).
#' when it's TURE, it return only (bcr_patient_barcode, days_to_last_followup, vital_status,days_to_new_tumor_event, tumor_event_status, tumor_event_free_days)
#' which is latest vital data (because the latest follow up vital data was integrated if avaliable) ready for survival analysis.
#'
#' @author Fan Zhang
#'
#' @return If simplify=TRUE a data frame ready for survival analysis of OS and PFS, with vital information and new tumor event free times, where the new follow up vital data is integrated. As simplify=FALSE, a data frame contain all vital, follow up times and tumor event days will return.
#' @export
#'
#' @examples
#' #list <- list.files(".","*.xml",recursive = TRUE)
#' #vitalParse(list,simplify=TRUE)
vitalParse <- function(list,simplify=TRUE){


  if(!all(grepl(".xml",list) == TRUE)){
    stop("list must be a vector of TCGA XML file names")}

  #extract the TCGA xml files one by one

  vital_infor <- lapply(list,
                        function(filenames){

                          file <- XML::xmlParse(filenames,encoding="UTF-8")

                          file <- XML::xmlToList(file)

                          patient <- file[[2]]

                          rm(file)


                          #check NA function
                          checkNA <- function(x){
                            #if extract fail, return NA
                            out_put <- tryCatch(x,error=function(c){NA})
                            #if extract null, return NA
                            out_put <- ifelse(is.null(out_put),NA,out_put)
                            return(out_put)
                          }
                          ##patient basic ID infor
                          patient_infor <- data.frame(bcr_patient_barcode = checkNA(patient[["bcr_patient_barcode"]][["text"]]),
                                                      vital_status = checkNA(patient[["vital_status"]][["text"]]),
                                                      days_to_death = checkNA(patient[["days_to_death"]][["text"]]),
                                                      days_to_last_followup = checkNA(patient[["days_to_last_followup"]][["text"]]),
                                                      days_to_new_tumor_event=checkNA(patient[["new_tumor_events"]][["new_tumor_event"]][["days_to_new_tumor_event_after_initial_treatment"]][["text"]])
                          )





                          #follow up data
                          #IF have follow up data
                          if(!is.null(patient[["follow_ups"]])){
                            #extract every follow up tumor event YSE or NOT,return a vector TRUE for tumor event happened
                            tumor_event_YESorNOT <- sapply(patient[["follow_ups"]],function(x){
                              x[["new_tumor_events"]][["new_tumor_event_after_initial_treatment"]][["text"]]=="YES"
                            })

                            #IF tumor_event have any yes
                            if(any(tumor_event_YESorNOT)){
                              #For all the follow up select all tumor_event time and get the first tumor event time min()
                              follow_up_min_tumor_event_time  <- lapply(patient[["follow_ups"]][tumor_event_YESorNOT],function(x){
                                min(x[["new_tumor_events"]][["new_tumor_event"]][["days_to_new_tumor_event_after_initial_treatment"]][["text"]])
                              })

                              #select the min tumor_event time in all follow up
                              followupInfor_days_to_new_tumor_event_after_initial_treatment <- min(as.numeric(unlist(follow_up_min_tumor_event_time)))

                              # new_tumor_event is yes
                              followupInfor_new_tumor_event_after_initial_treatment <- "YES"
                            } else {
                              #IF without tumor_event
                              followupInfor_new_tumor_event_after_initial_treatment <- checkNA(patient[["follow_ups"]][["follow_up"]][["new_tumor_events"]][["new_tumor_event_after_initial_treatment"]][["text"]])
                              followupInfor_days_to_new_tumor_event_after_initial_treatment <- NA
                            }

                            #select the latest follow up data by select the largest length of follow_ups infor: patient[["follow_ups"]][[length(patient[["follow_ups"]])]]
                            #followupInfor_barcode
                            #latest_followupInfor_days_to_death
                            #latest_followupInfor_days_to_last_followup
                            #latest_followupInfor_vital_status
                            followup_infor <- data.frame(followupInfor_barcode = checkNA(patient[["follow_ups"]][[length(patient[["follow_ups"]])]][["bcr_followup_barcode"]][["text"]]),
                                                         followupInfor_vital_status = checkNA(patient[["follow_ups"]][[length(patient[["follow_ups"]])]][["vital_status"]][["text"]]),
                                                         followupInfor_days_to_last_followup = checkNA(patient[["follow_ups"]][[length(patient[["follow_ups"]])]][["days_to_last_followup"]][["text"]]),
                                                         followupInfor_days_to_death = checkNA(patient[["follow_ups"]][[length(patient[["follow_ups"]])]][["days_to_death"]][["text"]]),
                                                         followupInfor_days_to_new_tumor_event_after_initial_treatment = followupInfor_days_to_new_tumor_event_after_initial_treatment,
                                                         followupInfor_new_tumor_event_after_initial_treatment = followupInfor_new_tumor_event_after_initial_treatment
                            )

                          } else {
                            #IF do not have follow up data: is.null(patient[["follow_ups"]])
                            followup_infor <- data.frame(followupInfor_barcode = NA,
                                                         followupInfor_vital_status = NA,
                                                         followupInfor_days_to_last_followup = NA,
                                                         followupInfor_days_to_death = NA,
                                                         followupInfor_days_to_new_tumor_event_after_initial_treatment = NA,
                                                         followupInfor_new_tumor_event_after_initial_treatment = NA
                            )

                          }

                          vital_infor <- cbind(patient_infor, followup_infor)


                          return(vital_infor)
                        })

  #convert list to data frame
  #vital_infor is hasn't simplify vital information
  vital_infor <- do.call(rbind,vital_infor)

  #Simplify=F
  if(simplify==F){
    return(vital_infor)
  }

  #Simplify the output for survival analysis
  vital_data <- apply(vital_infor,1,function(x){

    #binding
    if( !is.na(x["followupInfor_days_to_new_tumor_event_after_initial_treatment"])  ){
      x["days_to_new_tumor_event"] <- x["followupInfor_days_to_new_tumor_event_after_initial_treatment"]
    }


    if( !is.na(x["days_to_death"]) ){
      x["days_to_last_followup"] <- x["days_to_death"]
    }

    if( !is.na(x["followupInfor_days_to_death"]) ){
      x["followupInfor_days_to_last_followup"] <- x["followupInfor_days_to_death"]
    }


    x <- x[c(1,2,4,5,6,7,8,10,11)]

    if( !is.na(x["followupInfor_days_to_last_followup"]) ){
      x["vital_status"] <- x["followupInfor_vital_status"]
      x["days_to_last_followup"] <- x["followupInfor_days_to_last_followup"]
    }

    #if the patient without new tumor event, it's pregression free time will be last follow up time
    if(x["followupInfor_new_tumor_event_after_initial_treatment"]=="NO"){
      x["followupInfor_days_to_new_tumor_event_after_initial_treatment"] <- x["days_to_last_followup"]
    }

    x <- x[c(1:2,3,9,4)]

    names(x)[3:5] <- c("day_to_last_follow_OS","new_tumor_event_after_initial_treatment","days_Tumor_free_PFS")

    return(x)

  })


  vital_data <- as.data.frame(t(vital_data))

  #Simplify=TRUE
  return(vital_data)
}
