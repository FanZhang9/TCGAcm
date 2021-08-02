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

#'  Parse TCGA clinical vital data from XML for COX analysis
#'
#' @description `vitalParse()` integate TCGA clinical data from XML including the old vital data and updated latest follow up vital data. And the minimum days to new tumor event after initial treatment.
#'
#' For each value in `list` (in fact, it's a vector) `vitalParse()` extract the vital data from clinical XML files downloaded from the TCGA database
#' including vital_status and days to last follow up (OS), the latest follow_up vital data, which is updataed vital data and it can be easily omited. Meanwhile the days to new tumor event (PFS) is the earliest days to tumor event, after_initial_treatment.
#' So we can use this OS, PFS to perform COX analysis. return data including bcr_patient_barcode, tumor_tissue_site, vital_status, days_to_death, days_to_last_followup, new_tumor_event_after_initial_treatment, days_to_new_tumor_event_after_initial_treatment, days_tumor_free_till_last_followup
#'
#' @param filenames A vector contain the names of TCGA clinical XML files directorys and names (tcga_files_path/xml_file_names), which can easily retrieve by `list.files(path = ".", pattern = "*.xml",recursive = TRUE)`. The working directorys is the key, please check the names (which contain tcga_files_path/XML_file_names) is correct.
#'
#' @param dir file directory contained the XML files, default is current directory "."
#'
#'
#' @author Fan Zhang <fanzhang95@outlook.com>
#'
#'
#'
vitalParse <- function(filenames, dir="."){


  if(!all(grepl(".xml",filenames))){
    stop("list must be a vector of TCGA XML file names")}

  if(!file.exists(file.path(dir,filenames[1]))){
    stop(paste(file.path(dir,filenames[1]),"cann't not find the xml files"))
  }

  #extract the TCGA xml files one by one

  vital_infor <- lapply(filenames,
                        function(filenames){

                          xfile <- xml2::read_xml(file.path(dir,filenames))


                          #xml2::xml_text(xml2::xml_find_all(xfile,"//nte:new_tumor_event_after_initial_treatment"))

                          #xml2::xml_text(xml2::xml_find_all(xfile,"//nte:days_to_new_tumor_event_after_initial_treatment"))



                          #xml2::xml_text(xml2::xml_find_all(xfile,"//nte:days_to_new_tumor_event_after_initial_treatment"))

                          #xml2::xml_text(xml2::xml_find_all(xfile,"//cesc_nte:new_tumor_events"))

                          #xml2::xml_text(xml2::xml_find_all(xfile,"//cesc:follow_ups/*/clin_shared:bcr_followup_barcode"))


                          vital <- tibble::tibble(vital_status = xml2::xml_text(xml2::xml_find_all(xfile,"//clin_shared:vital_status")),
                                                  days_to_death = xml2::xml_text(xml2::xml_find_all(xfile,"//clin_shared:days_to_death")),
                                                  days_to_last_followup = xml2::xml_text(xml2::xml_find_all(xfile,"//clin_shared:days_to_last_followup")))

                          vital <- dplyr::arrange(vital,dplyr::desc(days_to_last_followup),decreasing = TRUE)
                          vital <- vital[1,]

                          new_tumor_event_after_initial_treatment <- xml2::xml_text(xml2::xml_find_all(xfile,"//cesc:follow_ups/*/nte:new_tumor_event_after_initial_treatment | //cesc_nte:new_tumor_events/nte:new_tumor_event_after_initial_treatment"))
                          if(length(new_tumor_event_after_initial_treatment)==0)
                            new_tumor_event_after_initial_treatment <- xml2::xml_text(xml2::xml_find_all(xfile,"//nte:new_tumor_event_after_initial_treatment"))

                          days_to_new_tumor_event_after_initial_treatment <- xml2::xml_text(xml2::xml_find_all(xfile,"//cesc:follow_ups/*/nte:days_to_new_tumor_event_after_initial_treatment | //cesc_nte:new_tumor_events/nte:new_tumor_event_after_initial_treatment"))
                          if(length(days_to_new_tumor_event_after_initial_treatment)==0)
                            days_to_new_tumor_event_after_initial_treatment <- xml2::xml_text(xml2::xml_find_all(xfile,"//nte:days_to_new_tumor_event_after_initial_treatment"))

                          if(length(days_to_new_tumor_event_after_initial_treatment)==0)
                            days_to_new_tumor_event_after_initial_treatment <- rep("",length(new_tumor_event_after_initial_treatment))

                          if(length(days_to_new_tumor_event_after_initial_treatment)!=length(new_tumor_event_after_initial_treatment)){
                            new_tumor_event_after_initial_treatment_backup <- new_tumor_event_after_initial_treatment

                            new_tumor_event_after_initial_treatment[which(new_tumor_event_after_initial_treatment=="YES")] <- days_to_new_tumor_event_after_initial_treatment

                            days_to_new_tumor_event_after_initial_treatment <- gsub("NO","",new_tumor_event_after_initial_treatment)
                            new_tumor_event_after_initial_treatment <- new_tumor_event_after_initial_treatment_backup
                          }



                          tumor_event <- tibble::tibble(new_tumor_event_after_initial_treatment = new_tumor_event_after_initial_treatment,
                                                        days_to_new_tumor_event_after_initial_treatment = days_to_new_tumor_event_after_initial_treatment,
                                                        days_to_death = xml2::xml_text(xml2::xml_find_all(xfile,"//clin_shared:days_to_death")),
                                                        days_to_last_followup = xml2::xml_text(xml2::xml_find_all(xfile,"//clin_shared:days_to_last_followup"))
                          )

                          tumor_event <- dplyr::filter(tumor_event, new_tumor_event_after_initial_treatment!="")

                          if(!nrow(tumor_event)==0){
                            ###without tumor event, select the last follow up time
                            if(all(tumor_event$new_tumor_event_after_initial_treatment=="NO")){
                              tumor_event <- dplyr::filter(tumor_event,new_tumor_event_after_initial_treatment=="NO")
                              days_tumor_free_till_last_followup <- max(tumor_event$days_to_last_followup)
                              new_tumor_event_after_initial_treatment <- "NO"
                              days_to_new_tumor_event_after_initial_treatment <- NA
                            } else if(any(tumor_event$new_tumor_event_after_initial_treatment=="YES")){
                              ##have tumor event, select the first tumor event time
                              tumor_event <- dplyr::filter(tumor_event,new_tumor_event_after_initial_treatment=="YES")
                              days_tumor_free_till_last_followup <- NA
                              new_tumor_event_after_initial_treatment <- "YES"
                              days_to_new_tumor_event_after_initial_treatment <- min(tumor_event$days_to_new_tumor_event_after_initial_treatment)
                            }
                          } else { days_tumor_free_till_last_followup <- NA
                          new_tumor_event_after_initial_treatment <- NA
                          days_to_new_tumor_event_after_initial_treatment <- NA
                          }


                          patient_infor <- tibble::tibble(bcr_patient_barcode = xml2::xml_text(xml2::xml_find_first(xfile,"//shared:bcr_patient_barcode")),
                                                          tumor_tissue_site = xml2::xml_text(xml2::xml_find_first(xfile,"//clin_shared:tumor_tissue_site")),
                                                          vital_status = vital$vital_status,
                                                          days_to_death = vital$days_to_death,
                                                          days_to_last_followup = vital$days_to_last_followup,
                                                          new_tumor_event_after_initial_treatment = new_tumor_event_after_initial_treatment,
                                                          days_to_new_tumor_event_after_initial_treatment = days_to_new_tumor_event_after_initial_treatment,
                                                          days_tumor_free_till_last_followup = days_tumor_free_till_last_followup
                                                          )


  })


  vital_data <- do.call(rbind, vital_infor)



  return(vital_data)
}
