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
    this$browseXMASS$ClearPath()
    this$browseMz$ClearPath()
    
    #XMASS format
    if( sDataType == this$sInputDataTypes[1])
    {
      this$browseMSIFile$SetFileFilter("xml")
      this$browseMSIFile$SetLabel("XML path:")
      this$browseXMASS$SetEnabled(T)
      this$browseMz$SetEnabled(T)
      gWidgets2::visible(this$box_pixelres) <- T
      return()
    }
    
    #rMSI TAR format
    if( sDataType == this$sInputDataTypes[2] )
    {
      this$browseMSIFile$SetFileFilter("tar")
      this$browseMSIFile$SetLabel("rMSI image:")
      this$browseXMASS$SetEnabled(F)
      this$browseMz$SetEnabled(F)
      gWidgets2::visible(this$box_pixelres) <- F
      return()
    }
    
    #imzML format
    if(sDataType == this$sInputDataTypes[3])
    {
      this$browseMSIFile$SetFileFilter("imzML")
      this$browseMSIFile$SetLabel("imzML image:")
      this$browseXMASS$SetEnabled(F)
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
    
    this$procParamList$smoothing <- list()
    this$procParamList$alignment <- list()
    this$procParamList$calibration <- list()
    this$procParamList$spectraNormalization <- list()
    this$procParamList$peakpicking <- list()
    
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
    
    if(this$procParamList$data$source$type == "xmass")
    {
      this$procParamList$data$source$datapath <- this$browseXMASS$GetPath()
    }
    else
    {
      this$procParamList$data$source$datapath <- this$browseMSIFile$GetPath()
    }
    if( length(this$procParamList$data$source$datapath) == 0)
    {
      #Error no file selected
      gWidgets2::gmessage("Error: No file selected.", icon = "error")
      this$procParamList <- NULL
      return()
    }
    if( this$procParamList$data$source$type == "xmass" )
    {
      if( !dir.exists(this$procParamList$data$source$datapath) )
      {
        #Error dir does not exists
        gWidgets2::gmessage("Error: The selected XMASS directory does not exists.", icon = "error")
        this$procParamList <- NULL
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
          this$procParamList <- NULL
          return()
        }
      }
    }

    if( this$procParamList$data$source$type == "xmass" )
    {
      this$procParamList$data$source$xmlpath <- this$browseMSIFile$GetPath()
      if( length(this$procParamList$data$source$xmlpath) == 0)
      {
        #Error no file selected
        gWidgets2::gmessage("Error: No XML file selected.", icon = "error")
        this$procParamList <- NULL
        return()
      }
      for(i in 1:length(this$procParamList$data$source$xmlpath))
      {
        if( !file.exists(this$procParamList$data$source$xmlpath[i]) )
        {
          #Error file does not exists
          gWidgets2::gmessage(paste("Error: file", this$procParamList$data$source$xmlpath[i], "\nDoes not exists."), icon = "error")
          this$procParamList <- NULL
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
        this$procParamList <- NULL
        return()
      }
      #Error file does not exists
      if(!file.exists(this$procParamList$data$source$spectrumpath))
      {
        gWidgets2::gmessage(paste("Error: file", this$procParamList$data$source$spectrumpath, "\nDoes not exists."), icon = "error")
        this$procParamList <- NULL
        return()
      }
      
      this$procParamList$data$pixelsize <- gWidgets2::svalue(this$spin_pixelSize)
      if(this$procParamList$data$pixelsize == 0)
      {
        gWidgets2::gmessage(paste("Error: Pixel resolution must be set for XMASS data."), icon = "error")
        this$procParamList <- NULL
        return()
      }
    }
    this$procParamList$data$outpath <- this$browseOut$GetPath()
    if( is.null(this$procParamList$data$outpath) )
    {
      gWidgets2::gmessage("Error: No output directory selected.", icon = "error")
      this$procParamList <- NULL
      return()
    }
    if( !dir.exists(this$procParamList$data$outpath))
    {
      #Error file does not exists
      gWidgets2::gmessage("Error: output directory does not exists.", icon = "error")
      this$procParamList <- NULL
      return()
    }

    #Smoothing list
    this$procParamList$smoothing$enabled <- gWidgets2::svalue(this$check_smoothing)
    if(this$procParamList$smoothing$enabled)
    {
      this$procParamList$smoothing$sgkernsize <- as.integer(gWidgets2::svalue(this$ratio_SGkernSize))
    }
    
    #Alignment list
    this$procParamList$alignment$enabled  <- gWidgets2::svalue(this$check_alignment) 
    if( this$procParamList$alignment$enabled )
    {
      this$procParamList$alignment$iterations <- as.integer(gWidgets2::svalue(this$spin_AlignIterations))
      this$procParamList$alignment$maxshift <- as.double(gWidgets2::svalue(this$spin_AlignMaxDev))
    }
    
    #Calibration list
    this$procParamList$calibration$enabled <- gWidgets2::svalue(this$check_calibration)
    if( this$procParamList$calibration$enabled )
    {
      if(gWidgets2::svalue(this$check_peakpicking))
      {
        this$procParamList$calibration$winsize <- as.integer(gWidgets2::svalue(this$spin_peakWin))
      }
      else
      {
        this$procParamList$calibration$winsize <- 20
      }
    }
    
    this$procParamList$spectraNormalization$enabled <- (gWidgets2::svalue(this$check_TICnorm) | gWidgets2::svalue(this$check_MAXnorm) | gWidgets2::svalue(this$check_TicACQnorm))
    if(this$procParamList$spectraNormalization$enabled)
    {
      this$procParamList$spectraNormalization$TIC <- gWidgets2::svalue(this$check_TICnorm)
      this$procParamList$spectraNormalization$MAX <- gWidgets2::svalue(this$check_MAXnorm)
      this$procParamList$spectraNormalization$AcqTIC <- gWidgets2::svalue(this$check_TicACQnorm)
    }
    
    #Peak picking list
    this$procParamList$peakpicking$enabled <- gWidgets2::svalue(this$check_peakpicking)
    if(this$procParamList$peakpicking$enabled)
    {
      this$procParamList$peakpicking$snr <- as.double(gWidgets2::svalue(this$spin_SNR))
      this$procParamList$peakpicking$winsize <- as.integer(gWidgets2::svalue(this$spin_peakWin))
      this$procParamList$peakpicking$oversample <- as.integer(gWidgets2::svalue(this$spin_peakOversample))
      this$procParamList$peakpicking$bintolerance <- as.double(gWidgets2::svalue(this$spin_binTolerance))
      this$procParamList$peakpicking$binfilter <- as.double(gWidgets2::svalue(this$spin_binFilter)/100)
    }
    
    #Number of threads
    this$procParamList$nthreads <- as.integer( gWidgets2::svalue( this$spin_nThreads ))
    
    gWidgets2::dispose(h$obj)
  }
  
  #Checkboxes enabled status
  ChkBoxSmoothingChanged <- function(...)
  {
    gWidgets2::enabled( this$lblSg) <- gWidgets2::svalue(check_smoothing)
    gWidgets2::enabled( this$ratio_SGkernSize) <- gWidgets2::svalue(check_smoothing)
  }
  
  ChkBoxAlignmentChanged <- function(...)
  {
    gWidgets2::enabled( this$spin_AlignIterations) <- gWidgets2::svalue(check_alignment)
    gWidgets2::enabled( this$spin_AlignMaxDev) <- gWidgets2::svalue(check_alignment)
  }
  
  ChkBoxCalibrationChanged <- function(...)
  {
    #Nothing here currently, keeping it for future features...
  }
  
  ChkBoxPeakPickingChanged <- function(...)
  {
    gWidgets2::enabled( this$spin_SNR) <- gWidgets2::svalue(check_peakpicking)
    gWidgets2::enabled( this$spin_peakWin) <- gWidgets2::svalue(check_peakpicking)
    gWidgets2::enabled( this$spin_peakOversample) <- gWidgets2::svalue(check_peakpicking)
    gWidgets2::enabled( this$spin_binTolerance) <- gWidgets2::svalue(check_peakpicking)
    gWidgets2::enabled( this$spin_binFilter) <- gWidgets2::svalue(check_peakpicking)
  }
  
  mainWin<-gWidgets2::gwindow(title = "MSI data import and process wizard", visible = F)
  box_mainH <- gWidgets2::ggroup(horizontal = T, container = mainWin, expand = T, fill = T)
  box_mainV <- gWidgets2::ggroup(horizontal = F, container = box_mainH, expand = T, fill = T)
  
  #Data Input box
  frm_dataInput <-  gWidgets2::gframe( "Data Source", container =  box_mainV, expand = T, fill = T)
  box_dataInput <- gWidgets2::ggroup( horizontal = F, container = frm_dataInput, expand = T, fill = T)
  ratio_inDataType <- gWidgets2::gradio(items = sInputDataTypes, horizontal = T, selected = 2, container = box_dataInput, handler = this$DataInputTypeChanged )
  
  browseMSIFile <- FileBrowseWidget( box_dataInput, setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  browseXMASS <- FileBrowseWidget( box_dataInput, sLabel = "XMASS Dir:", dirSel = T, multiSel = F, setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  browseMz <- FileBrowseWidget( box_dataInput, sLabel = "Spectrum:", dirSel = F, fFilter = "txt", multiSel = F, setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  
  #The pixel size widget
  box_pixelres <- gWidgets2::ggroup(horizontal = T, container = box_dataInput)
  gWidgets2::glabel("Pixel resolution [um]:", container = box_pixelres )
  spin_pixelSize<-gWidgets2::gspinbutton(from = 0, to = 5000, by = 1, value = 0, container = box_pixelres, digits = 1)
  
  #Set initial state properly
  DataInputTypeChanged()
  
  #Pre-processing box
  frm_preProcessing <- gWidgets2::gframe( "Pre-Processing parameters", container = box_mainH, expand = T, fill = T )
  box_procH <- gWidgets2::ggroup(horizontal = T, container = frm_preProcessing, expand = T, fill = T, spacing = 20)
  box_proc1 <- gWidgets2::ggroup(horizontal = F, container = box_procH, expand = T, fill = T, spacing = 20)
  box_proc2 <- gWidgets2::ggroup(horizontal = F, container = box_procH, expand = T, fill = T, spacing = 20)
  drawLabelSpin <- function( parent_widget, sText, minVal, maxVal, defaultVal, decPlaces = 0 )
  {
    box_spin <- gWidgets2::ggroup(horizontal = T, container = parent_widget)
    gWidgets2::glabel(sText, container = box_spin )
    gWidgets2::addSpring(box_spin)
    return (gWidgets2::gspinbutton(from = minVal, to = maxVal, by = 1, value = defaultVal, container = box_spin, digits = decPlaces)  )
  }
  
  #Smoothing params
  frm_SGKernSize <- gWidgets2::gframe("Savitzky-Golay Smoothing:", container = box_proc1)
  box_smoothing <- gWidgets2::ggroup(horizontal = F, container = frm_SGKernSize)
  check_smoothing <- gWidgets2::gcheckbox("Enable smoothing", checked = T, container = box_smoothing, handler = this$ChkBoxSmoothingChanged )
  lblSg <- gWidgets2::glabel("Savitzky-Golay kernel size:", container = box_smoothing)
  ratio_SGkernSize<- gWidgets2::gradio(c(5,7,9,11,13,15), selected = 2, horizontal = T, container = box_smoothing)
    
  #Alignment params
  frm_alignment <- gWidgets2::gframe("Alignment", container = box_proc1, spacing = 10)
  box_alignment <- gWidgets2::ggroup(horizontal = F, container = frm_alignment)
  check_alignment <- gWidgets2::gcheckbox("Enable alignment", checked = T, container = box_alignment, handler = this$ChkBoxAlignmentChanged)
  spin_AlignIterations <- drawLabelSpin(box_alignment, "Iterations:", 1, 5, 3)
  spin_AlignMaxDev <- drawLabelSpin(box_alignment, "Max Shift [ppm]:", 10, 1000, 200)
  
  #Mass Calibration params
  frm_calibration <- gWidgets2::gframe("Mass Calibration", container = box_proc1, spacing = 10)
  box_calibration <- gWidgets2::ggroup(horizontal = F, container = frm_calibration)
  check_calibration <- gWidgets2::gcheckbox("Enable calibration", checked = T, container = box_calibration, handler = this$ChkBoxCalibrationChanged)
  
  #Spectra Normalization params
  frm_spectraNorm <- gWidgets2::gframe("Spectra intensity normalization", container = box_proc1, spacing = 10)
  box_spectraNorm <- gWidgets2::ggroup(horizontal = F, container = frm_spectraNorm)
  check_TICnorm <- gWidgets2::gcheckbox("TIC", checked = T, container = box_spectraNorm )
  check_MAXnorm <- gWidgets2::gcheckbox("MAX", checked = T, container = box_spectraNorm )
  check_TicACQnorm <- gWidgets2::gcheckbox("TICAcq", checked = T, container = box_spectraNorm )
  
  #Peak picking and binning params
  frm_peakpick <- gWidgets2::gframe("Peak-Picking", container = box_proc2, spacing = 10)
  box_peakpick <- gWidgets2::ggroup(horizontal = F, container = frm_peakpick)
  check_peakpicking <- gWidgets2::gcheckbox("Enable peak-picking", checked = T, container = box_peakpick, handler = this$ChkBoxPeakPickingChanged)
  spin_SNR <- drawLabelSpin(box_peakpick, "SNR Threshold:", 1, 100, 5)
  spin_peakWin <- drawLabelSpin(box_peakpick, "Detector window size:", 5, 200, 20)
  spin_peakOversample <- drawLabelSpin(box_peakpick, "Peak shape over-sampling:", 1, 50, 10)
  spin_binTolerance <- drawLabelSpin(box_peakpick, "Peak-Bin Tolerance [ppm]:", 1, 1000, 100, decPlaces = 0)
  spin_binFilter <- drawLabelSpin(box_peakpick, "Peak Filter [%]:", 1, 100, 5)
  
  #Number of processing threads
  frm_procThreads <- gWidgets2::gframe("Processing Threads", container = box_proc2, spacing = 10)
  spin_nThreads <- drawLabelSpin(frm_procThreads, "Max Threads:", 1, parallel::detectCores(), parallel::detectCores(), decPlaces = 0)
  
  #Data output box
  gWidgets2::addSpring(box_mainV)
  frm_dataOutput <- gWidgets2::gframe( "Data Output", container = box_mainV, expand = T, fill = T)
  browseOut <- FileBrowseWidget( frm_dataOutput, sLabel = "Output directory:", dirSel = T, fFilter = "txt", setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  
  ##Run button
  btn_run <- gWidgets2::gbutton("Process Data", handler = this$RunClicked, container = box_proc2)
  
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
