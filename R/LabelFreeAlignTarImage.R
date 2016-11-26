#' LabelFreeAlign_Wizard.
#' 
#' Easy to use Label-Free aligment tool.
#' Just run this function and a wizard will guide you trought the alignment process. 
#'
#' @param storeShifts TRUE if applied shifts must be recorded.
#'
#' @return if storeShifts == TRUE then a list of applied shifts is returned.
#' @export
#'
LabelFreeAlign_Wizard<-function( storeShifts = F)
{
  startWD <- getwd()

  #Ask for the data path
  cat("Select a tar image to align...\n\n")
  path_in_img <- gWidgets2::gfile(text = "Select a tar image to align", type = "open", filter = c("tar" = "tar"), multi = T)
  if(length(path_in_img) == 0){
    setwd(startWD)
    stop("Data importation aborted by user, please source this file again to restart it!\n", call. = F)
  }
  setwd(dirname(path_in_img))

  path_out_img <- paste(tools::file_path_sans_ext(path_in_img),"_Aligned",".tar",sep = "" )

  if(storeShifts)
  {
    Shifts <- list()
  }
  
  for( i in 1:length(path_in_img))
  {
    cat( paste("Starting aligning image:", path_in_img[i], "  ", i ,"of", length(path_in_img), "\n"))

    #Load the input image
    raw<-rMSI::LoadMsiData(path_in_img[i])
    
    #Algin the image
    mShifts <- FullDataSetAligment(raw, subDataMemPercent = 1)
    if(storeShifts)
    {
      Shifts[[length(Shifts) + 1]] <- mShifts
    }

    #Saving Data to compressed format
    rMSI::SaveMsiData(data_file = path_out_img[i], imgData = raw)
    
    #Remove the ramdisk
    rMSI::DeleteRamdisk(raw)

    cat( paste("Image:", path_in_img[i], "Aligned to:", path_out_img[i], "\n\n\n"))
  }

  setwd(startWD)
  
  if(storeShifts)
  {
    return (Shifts)
  }
}
