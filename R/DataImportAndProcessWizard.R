#########################################################################
#     rMSIproc - R package for MSI data processing
#     Copyright (C) 2014 Pere Rafols Soler
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

ImportWizardGui <- function()
{
  oldWarning<-options()$warn
  options(warn = -1)
  
  ## Get the environment for this
  ## instance of the function.
  this <- environment()
  
  ##Class data members 
  sInputDataTypes <- c("XMASS", "rMSI(.tar)", "imzML")
  procParamList <- NULL
  workdir <- path.expand("~/")
  
  ##Class methods
  DataInputTypeChanged <- function(...)
  {
    sDataType <-  gWidgets2::svalue(this$ratio_inDataType)
    this$browseMSIFile$ClearPath()
    this$browseXML$ClearPath()
    this$browseMz$ClearPath()
    
    #XMASS format
    if( sDataType == this$sInputDataTypes[1])
    {
      this$browseMSIFile$SetDirSelection(T)
      this$browseMSIFile$SetFileFilter("xmass")
      this$browseMSIFile$SetLabel("XMASS path:")
      this$browseXML$SetEnabled(T)
      this$browseMz$SetEnabled(T)
      gWidgets2::visible(this$box_pixelres) <- T
      return()
    }
    
    #rMSI TAR format
    if( sDataType == this$sInputDataTypes[2] )
    {
      this$browseMSIFile$SetDirSelection(F)
      this$browseMSIFile$SetFileFilter("tar")
      this$browseMSIFile$SetLabel("rMSI image:")
      this$browseXML$SetEnabled(F)
      this$browseMz$SetEnabled(F)
      gWidgets2::visible(this$box_pixelres) <- F
      return()
    }
    
    #imzML format
    if(sDataType == this$sInputDataTypes[3])
    {
      this$browseMSIFile$SetDirSelection(F)
      this$browseMSIFile$SetFileFilter("imzML")
      this$browseMSIFile$SetLabel("imzML image:")
      this$browseXML$SetEnabled(F)
      this$browseMz$SetEnabled(F)
      gWidgets2::visible(this$box_pixelres) <- F
      return()
    }
    
    stop("Error in ImportWizardGui::DataInputTypeChanged invalid datatype\n")
  }
  
  SetWorkingDir <- function( sDir )
  {
    this$workdir <- sDir
  }
  
  GetWorkingDir <- function()
  {
    return(this$workdir)
  }
  
  RunClicked <- function(h, ...)
  {
    #List tree
    this$procParamList <- list()
    this$procParamList$data <- list()
    this$procParamList$alignment <- list()
    this$procParamList$peakpicking <- list()
    this$procParamList$binning <- list()
    
    #Data source and output list
    this$procParamList$data$source <- list()
    sDataType <-  gWidgets2::svalue(this$ratio_inDataType)
    if( sDataType == this$sInputDataTypes[1] )
    {
      #XMASS
      this$procParamList$data$source$type <- "xmass"
    }
    if( sDataType == this$sInputDataTypes[2] )
    {
      #TAR
      this$procParamList$data$source$type <- "tar"
    }
    if( sDataType == this$sInputDataTypes[3] )
    {
      #imzML
      this$procParamList$data$source$type <- "imzml"
    }
    
    this$procParamList$data$source$datapath <- this$browseMSIFile$GetPath()
    if( length(this$procParamList$data$source$datapath) == 0)
    {
      #Error no file selected
      gWidgets2::gmessage("Error: No file selected.", icon = "error")
      return()
    }
    if( this$procParamList$data$source$type == "xmass" )
    {
      if( !dir.exists(this$procParamList$data$source$datapath) )
      {
        #Error dir does not exists
        gWidgets2::gmessage("Error: The selected XMASS directory does not exists.", icon = "error")
        return()
      }
    }
    else
    {
      for(i in 1:length(this$procParamList$data$source$datapath))
      {
        if( !file.exists(this$procParamList$data$source$datapath[i]) )
        {
          #Error file does not exists
          gWidgets2::gmessage(paste("Error: file", this$procParamList$data$source$datapath[i], "\nDoes not exists."), icon = "error")
          return()
        }
      }
    }

    this$procParamList$data$source$xmlpath <- this$browseXML$GetPath()
    if( this$procParamList$data$source$type == "xmass" )
    {
      if( length(this$procParamList$data$source$xmlpath) == 0)
      {
        #Error no file selected
        gWidgets2::gmessage("Error: No XML file selected.", icon = "error")
        return()
      }
      for(i in 1:length(this$procParamList$data$source$xmlpath))
      {
        if( !file.exists(this$procParamList$data$source$xmlpath[i]) )
        {
          #Error file does not exists
          gWidgets2::gmessage(paste("Error: file", this$procParamList$data$source$xmlpath[i], "\nDoes not exists."), icon = "error")
          return()
        }
      }
    }
    this$procParamList$data$source$spectrumpath <- this$browseMz$GetPath()
    if( this$procParamList$data$source$type == "xmass")
    {
      if( is.null(this$procParamList$data$source$spectrumpath) )
      {
        gWidgets2::gmessage("Error: No spectrum file selected.", icon = "error")
        return()
      }
      #Error file does not exists
      if(!file.exists(this$procParamList$data$source$spectrumpath))
      {
        gWidgets2::gmessage(paste("Error: file", this$procParamList$data$source$spectrumpath, "\nDoes not exists."), icon = "error")
        return()
      }
      
      this$procParamList$data$pixelsize <- gWidgets2::svalue(this$spin_pixelSize)
      if(this$procParamList$data$pixelsize == 0)
      {
        gWidgets2::gmessage(paste("Error: Pixel resolution must be set for XMASS data."), icon = "error")
        return()
      }
    }
    this$procParamList$data$outpath <- this$browseOut$GetPath()
    if( is.null(this$procParamList$data$outpath) )
    {
      gWidgets2::gmessage("Error: No output directory selected.", icon = "error")
      return()
    }
    if( !dir.exists(this$procParamList$data$outpath))
    {
      #Error file does not exists
      gWidgets2::gmessage("Error: output directory does not exists.", icon = "error")
      return()
    }
    
    
    #Alignment list
    this$procParamList$alignment$iterations <- as.integer(gWidgets2::svalue(this$spin_AlignIterations))
    this$procParamList$alignment$maxshift <- as.double(gWidgets2::svalue(this$spin_AlignMaxDev))
    
    #Peak picking list
    this$procParamList$peakpicking$snr <- as.double(gWidgets2::svalue(this$spin_SNR))
    this$procParamList$peakpicking$winsize <- as.integer(gWidgets2::svalue(this$spin_peakWin))
    this$procParamList$peakpicking$oversample <- as.integer(gWidgets2::svalue(this$spin_peakOversample))
    this$procParamList$peakpicking$sgkernsize <- as.integer(gWidgets2::svalue(this$ratio_SGkernSize))
    
    #Binning list
    this$procParamList$binning$tolerance <- as.double(gWidgets2::svalue(this$spin_binTolerance))
    this$procParamList$binning$filter <- as.double(gWidgets2::svalue(this$spin_binFilter)/100)
    
    gWidgets2::dispose(h$obj)
  }
  
  mainWin<-gWidgets2::gwindow(title = "MSI data import and process wizard", visible = F)
  box_main <- gWidgets2::ggroup(horizontal = F, container = mainWin, expand = T, fill = T)
  
  #Data Input box
  frm_dataInput <-  gWidgets2::gframe( "Data Source", container =  box_main, expand = T, fill = T)
  box_dataInput <- gWidgets2::ggroup( horizontal = F, container = frm_dataInput, expand = T, fill = T)
  ratio_inDataType <- gWidgets2::gradio(items = sInputDataTypes, horizontal = T, selected = 2, container = box_dataInput, handler = this$DataInputTypeChanged )
  
  browseMSIFile <- FileBrowseWidget( box_dataInput, setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  browseXML <- FileBrowseWidget( box_dataInput, sLabel = "XML File:", dirSel = F, fFilter = "xml", setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  browseMz <- FileBrowseWidget( box_dataInput, sLabel = "Spectrum:", dirSel = F, fFilter = "txt", multiSel = F, setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  
  #The pixel size widget
  box_pixelres <- gWidgets2::ggroup(horizontal = T, container = box_dataInput)
  gWidgets2::glabel("Pixel resolution [um]:", container = box_pixelres )
  spin_pixelSize<-gWidgets2::gspinbutton(from = 0, to = 5000, by = 1, value = 0, container = box_pixelres, digits = 1)
  
  #Set initial state properly
  DataInputTypeChanged()
  
  #Pre-processing box
  frm_preProcessing <- gWidgets2::gframe( "Pre-Processing parameters", container = box_main, expand = T, fill = T )
  box_proc <- gWidgets2::ggroup(horizontal = T, container = frm_preProcessing, expand = T, fill = T, spacing = 20)
  drawLabelSpin <- function( parent_widget, sText, minVal, maxVal, defaultVal, decPlaces = 0 )
  {
    box_spin <- gWidgets2::ggroup(horizontal = T, container = parent_widget)
    gWidgets2::glabel(sText, container = box_spin )
    gWidgets2::addSpring(box_spin)
    return (gWidgets2::gspinbutton(from = minVal, to = maxVal, by = 1, value = defaultVal, container = box_spin, digits = decPlaces)  )
  }
  
  #Alignment params
  frm_alignment <- gWidgets2::gframe("Alignment", container = box_proc, spacing = 10)
  box_alignment <- gWidgets2::ggroup(horizontal = F, container = frm_alignment)
  spin_AlignIterations <- drawLabelSpin(box_alignment, "Iterations:", 0, 5, 3)
  spin_AlignMaxDev <- drawLabelSpin(box_alignment, "Max Shift [ppm]:", 10, 1000, 200)
  
  #Peak picking params
  frm_peakpick <- gWidgets2::gframe("Peak-Picking", container = box_proc, spacing = 10)
  box_peakpick <- gWidgets2::ggroup(horizontal = F, container = frm_peakpick)
  spin_SNR <- drawLabelSpin(box_peakpick, "SNR Threshold:", 1, 50, 5)
  spin_peakWin <- drawLabelSpin(box_peakpick, "Detector window size:", 5, 50, 10)
  spin_peakOversample <- drawLabelSpin(box_peakpick, "Peak shape over-sampling:", 1, 50, 10)
  frm_SGKernSize <- gWidgets2::gframe("Savitzky-Golay kernel size:", container = box_peakpick)
  ratio_SGkernSize<- gWidgets2::gradio(c(0,5,7,9,11), horizontal = T, container = frm_SGKernSize)
  
  #Binning params
  frm_peakbinning <- gWidgets2::gframe("Peak-Binning", container = box_proc, spacing = 10)
  box_peakbinning <- gWidgets2::ggroup(horizontal = F, container = frm_peakbinning)
  spin_binTolerance <- drawLabelSpin(box_peakbinning, "Tolerance:", 0.01, 1, 0.05, decPlaces = 2)
  spin_binFilter <- drawLabelSpin(box_peakbinning, "Filter [%] :", 1, 100, 5)
  
  #Data output box
  frm_dataOutput <- gWidgets2::gframe( "Data Output", container = box_main, expand = T, fill = T)
  browseOut <- FileBrowseWidget( frm_dataOutput, sLabel = "Output directory:", dirSel = T, fFilter = "txt", setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  
  ##Run button
  box_run <- gWidgets2::ggroup(horizontal = T, container = box_main)
  gWidgets2::addSpring(box_run)
  btn_run <- gWidgets2::gbutton("Process Data", handler = this$RunClicked, container = box_run)
  
  gWidgets2::visible(mainWin) <- T
  
  ## Set the name for the class
  class(this) <- append(class(this),"DataProcessWindow")
  gc()
  
  #Do not return until this window is disposed...
  while(gWidgets2::isExtant(this$mainWin ))
  {
    Sys.sleep(0.1)
  }
  
  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)
  
  
  #Return a structured list with all input data
  return( procParamList )
}
