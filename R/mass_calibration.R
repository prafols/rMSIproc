#########################################################################
#     rMSI - R package for MSI data handling and visualization
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


#' CalibrationWindow.
#'
#' @param mass The mass vector of spectrum to calibrate.
#' @param intensity The intensity vector of spectrum to calibrate.
#' @param use_zoo if the zoo package interpolation must be used
#' @param method a string with the method used for interpolation, valid methods are spline, linear or nocal if only target masses must be returned.
#'
#' @return a the calibrated mass axis.
#' @export
#'
CalibrationWindow<-function( mass, intensity, peak_win_size = 20, win_title = "", method = "loess")
{
  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)
  
  ## Get the environment for this
  ## instance of the function.
  this <- environment()
  
  ##Class data members
  spectraWidget <- NULL
  dMass <- mass
  dIntensity <- intensity
  refMz <- NULL #a vector of user selected reference masses
  targetMz <- NULL #a vector of user selected target masses to calibrate
  PeakWindow <- peak_win_size
  rm(mass)
  rm(intensity)
  
  Tbl_ColNames <- list(name = "Name", ref = "Ref. m/z", sel = "Sel. m/z", err = "Error [m/z]", ppm = "Error [ppm]", active = "Active")
  dMassCalibrated <- NULL
  CalMethod <- method
  
  #Create an mzTable
  CreateMzTable <- function( ref_names, ref_mz )
  {
    sel_mz <- rep(NaN, length(ref_mz))
    mzTable <- data.frame( ref_names, ref_mz, sel_mz, (ref_mz - sel_mz), (1e6*(ref_mz - sel_mz)/ref_mz), rep(F, length(ref_mz)))
    colnames(mzTable) <- as.vector(unlist(this$Tbl_ColNames))
    return(mzTable)
  }
  
  #Enable calibration button only if are enough mz points to calibrate
  SetCalibrateButtonActiveState <- function()
  {
    if(length(which( this$Table_Ctl[ , this$Tbl_ColNames$active] ) ) >= 2)
    {
      gWidgets2::enabled(Btn_Calibrate) <-T
    }
    else
    {
      gWidgets2::enabled(Btn_Calibrate) <-F
    }
  }
  
  #A click on mass spectra widgets drives here. 
  SpectrumClicked <- function( channel, mass, tol )
  {
    if(tol < 0.5)
    {
      tol <- 0.5
    }
    idLo <- which.min( abs( (mass - tol) - this$dMass  ) )
    idHi <- which.min( abs( (mass + tol) - this$dMass  ) )
    subMass <- dMass[idLo:idHi]
    subIntensity <- this$dIntensity[idLo:idHi]
    pks <- DetectPeaks(subMass, subIntensity, SNR = 1, WinSize = this$PeakWindow)
    if( length(pks$intensity) == 0 )
    {
      return()
    }
    
    nearestPeak <- which.min(abs(mass - pks$mass))
    PeakMz <- pks$mass[nearestPeak]
    
    iRow <- this$Table_Ctl$get_selected()
    this$Table_Ctl[iRow, this$Tbl_ColNames$sel] <- PeakMz
    this$Table_Ctl[iRow, this$Tbl_ColNames$err] <- this$Table_Ctl[iRow, this$Tbl_ColNames$ref] - PeakMz
    this$Table_Ctl[iRow, this$Tbl_ColNames$ppm] <- 1e6*(this$Table_Ctl[iRow, this$Tbl_ColNames$ref] - PeakMz)/this$Table_Ctl[iRow, this$Tbl_ColNames$ref]
    this$Table_Ctl[iRow, this$Tbl_ColNames$active] <- T
    gWidgets2::svalue(this$Btn_Active) <-  this$Table_Ctl[iRow, this$Tbl_ColNames$active]
    gWidgets2::enabled(Btn_Active) <-T
    this$SetCalibrateButtonActiveState()
    
    this$spectraWidget$SetSelectedMassTol(1, PeakMz, 0)
  }
  
  #Load an ASCII with reference masses
  LoadRefMzAscii <- function (...)
  {
    fname <- gWidgets2::gfile("Select reference m/z file", type = "open", multi = F)
    if( length(fname) > 0)
    {
      dataM <- read.table(fname)
      ref_names <- as.vector(dataM[,1])
      ref_mz <- as.vector(dataM[,2])
      this$Table_Ctl$set_items( CreateMzTable( ref_names, ref_mz ))
      gWidgets2::size(this$Table_Ctl) <- list( width = -1, height = -1,  column.widths = rep(110, length(this$Tbl_ColNames)))
      this$spectraWidget$SetRefMass(ref_mz)
      gWidgets2::enabled(Table_Ctl) <- T
      gWidgets2::enabled(Btn_Active) <- T
      this$SetCalibrateButtonActiveState()
      this$spectraWidget$SetActiveTool("Red")
    }
  }
  
  #Mz Table row selected
  RowSelected <- function ( ... )
  {
    #Zoom in spectrum
    iRow <- this$Table_Ctl$get_selected()
    SelMz <- this$Table_Ctl[ iRow, this$Tbl_ColNames$ref]
    if( SelMz > max(this$dMass) || SelMz < min(this$dMass))
    {
      cat(paste("Selected mass (", SelMz, ") is out of range: ", min(this$dMass), " - " , max(this$dMass),"\n", sep = ""))
      return()
    }
    this$spectraWidget$ZoomMzRange( SelMz - 0.007*SelMz, SelMz + 0.007*SelMz )
    RGtk2::gtkButtonSetLabel( gWidgets2::getToolkitWidget(this$Btn_Active) , paste("m/z",sprintf("%.2f",SelMz), "active") )
    gWidgets2::svalue(this$Btn_Active) <-  this$Table_Ctl[iRow, this$Tbl_ColNames$active]
    
    if(is.nan( this$Table_Ctl[ iRow, this$Tbl_ColNames$sel] ))
    {
      gWidgets2::enabled(Btn_Active) <-F
    }
    else
    {
      gWidgets2::enabled(Btn_Active) <-T
    }
  }
  
  #Button to enable/disable mz
  BtnActiveChanged <- function( ... )
  {
    iRow <- this$Table_Ctl$get_selected()
    if( length(iRow) > 0)
    {
      this$Table_Ctl[iRow, this$Tbl_ColNames$active] <- gWidgets2::svalue(this$Btn_Active)
      this$SetCalibrateButtonActiveState()
    }
  }
  
  #Calibration Button clicked
  BtnCalibrate <- function( ... )
  {
    valid_rows <- which( this$Table_Ctl[, this$Tbl_ColNames$active] )
    this$refMz <- this$Table_Ctl[valid_rows, this$Tbl_ColNames$ref]
    this$targetMz <- this$Table_Ctl[valid_rows, this$Tbl_ColNames$sel]
    if(this$CalMethod == "nocal")
    {
      gWidgets2::dispose(this$this$window)
    }
    else
    {
      this$dMassCalibrated <- calMzAxis(this$dMass, this$refMz, this$targetMz, this$CalMethod)
      gWidgets2::enabled(Chk_ShowCal) <-T
      
      this$spectraWidget$AddSpectra(  this$dMassCalibrated, dIntensity, col = "darkgreen", name = "cal")
      gWidgets2::svalue(this$Chk_ShowRaw) <- F
      gWidgets2::svalue(this$Chk_ShowCal) <- T
      gWidgets2::enabled(Btn_Confirm) <- T
      this$spectraWidget$ZoomResetClicked()
      this$spectraWidget$SetActiveTool("Zoom")
    }
  }
  
  #Clicked on checkbox
  ChkShowRawSpc <- function( ... )
  {
    this$spectraWidget$SetSpectrumEnabled("raw", gWidgets2::svalue(this$Chk_ShowRaw))
  }
  
  #Clicked on checkbox
  ChkShowCalSpc <- function( ... )
  {
    this$spectraWidget$SetSpectrumEnabled("cal", gWidgets2::svalue(this$Chk_ShowCal))
  }
  
  #Validate calibration and quit
  ValidateCalAndQuit <- function( ... )
  {
    gWidgets2::dispose(this$this$window)
  }
  
  #GUI builder
  this$window <- gWidgets2::gwindow ( paste("Spectrum Calibration -",win_title ), visible = F )
    Grp_Top <- gWidgets2::gpanedgroup(horizontal = F, container = window)
  
  Grp_Bot <- gWidgets2::ggroup(horizontal = T, container = Grp_Top)
  Table_Ctl <- gWidgets2::gtable(CreateMzTable( c(NaN), c(NaN) ), multiple = F, container = Grp_Bot, chosen.col = 2) #First empty MZ table
  
  Table_Ctl$set_editable(F)
  gWidgets2::size(this$Table_Ctl) <- list( width = -1, height = -1,  column.widths = rep(110, length(Tbl_ColNames)))
  gWidgets2::enabled(Table_Ctl) <-F
  
  Grp_Btn <- gWidgets2::ggroup(horizontal = F, container = Grp_Bot)
  Btn_LoadRef <- gWidgets2::gbutton("Load Ref m/z", container = Grp_Btn, handler = this$LoadRefMzAscii)
  Btn_Active <- gWidgets2::gcheckbox("m/z active", checked = F, use.togglebutton = T, handler = this$BtnActiveChanged, container = Grp_Btn)
  Btn_Calibrate <- gWidgets2::gbutton("Calibrate", container = Grp_Btn, handler = this$BtnCalibrate)
  gWidgets2::enabled(Btn_Active) <-F
  gWidgets2::enabled(Btn_Calibrate) <-F
  
  Frm_PlotCtl <- gWidgets2::gframe("Plot control", container = Grp_Btn)
  Grp_PlotCtl <- gWidgets2::ggroup(horizontal = F, container = Frm_PlotCtl)
  Chk_ShowRaw <- gWidgets2::gcheckbox("RAW", checked = T, handler = this$ChkShowRawSpc, container = Grp_PlotCtl)
  Chk_ShowCal <- gWidgets2::gcheckbox("CAL", checked = F, handler = this$ChkShowCalSpc, container = Grp_PlotCtl)
  rMSI:::.setCheckBoxText(Chk_ShowRaw, "RAW spectrum", background = NULL, foreground = "darkblue", font_size = NULL, font_weight = "heavy")
  rMSI:::.setCheckBoxText(Chk_ShowCal, "CAL spectrum", background = NULL, foreground = "darkgreen", font_size = NULL, font_weight = "heavy")
  
  Btn_Confirm <- gWidgets2::gbutton("Validate CAL & quit", container = Grp_Btn, handler = this$ValidateCalAndQuit)
  gWidgets2::enabled(Btn_Confirm) <-F
  
  gWidgets2::enabled(Chk_ShowCal) <-F
  
  spectraFrame<-gWidgets2::gframe("", container = Grp_Top,  fill = T, expand = T, spacing = 5 )
  spectraWidget<-rMSI:::.SpectraPlotWidget(parent_widget = spectraFrame, top_window_widget = window, clicFuntion = this$SpectrumClicked, showOpenFileButton = F,  display_sel_red = T, display_sel_spins = F)
  
  gWidgets2::size( window )<- c(1024, 740)
  gWidgets2::size( spectraWidget )<- c(-1, 380)
  gWidgets2::size( Table_Ctl )<- c(-1, 200)
  
  visible(window)<-TRUE
  
  spectraWidget$AddSpectra(  dMass, dIntensity, col = "blue", name = "raw")
  spectraWidget$ZoomResetClicked()
  window$widget$present()
  Table_Ctl$add_handler_selection_changed(this$RowSelected)
  
  ## Set the name for the class
  class(this) <- append(class(this),"CalWindow")
  gc()
  
  #Do not return until this window is disposed...
  while(gWidgets2::isExtant(this$window ))
  {
    Sys.sleep(0.1)
  }
  
  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)
  
  if(this$CalMethod == "nocal")
  {
    #Return a list of masses and targets that can b used for calibration
    return(list( ref = this$refMz, target = this$targetMz))
  }
  else
  {
    #Return the mass axis... it is NULL if no calibration was done
    return(this$dMassCalibrated)
  }
}

#' Apply a new mass axis to a complete image.
#' A GUI to calibrate the mean spectrum will be shown.
#'
#' @param img A image in rMSI data format .
#' @param output_fname full path to store the output image.
#' @param newMzAxis if a new mz axis specified it is used for the calibration window
#' @param useZoo if zoo interpolation must be used (caution this may introduce m/z error if large compensations are made)
#'
#' @export
#'
CalibrateImage<-function(img, output_fname, newMzAxis = NULL)
{
  #Copy the img objet
  calImg <- img
  
  #Replace mass axis using GUI calibration
  if( class( calImg$mean) == "MassSpectrum")
  {
    #Old mean MALDIquant handling
    mIntensity <- calImg$mean@intensity
  }
  else
  {
    mIntensity <- calImg$mean
  }
  
  if(!is.null(newMzAxis))
  {
    if(length(newMzAxis) != length(calImg$mass))
    {
      stop("Eror: Mass axis length is different than intensity legnth!")
    }
    
    calImg$mass <- newMzAxis
    
    if( class( calImg$mean) == "MassSpectrum")
    {
      calImg$mean@mass <- newMzAxis
    }
  }
  
  new_mass <- CalibrationWindow( calImg$mass, mIntensity )
  
  if(is.null(new_mass))
  {
    cat("Calibration process aborted by user\n")
    return()
  }
  
  #Ask for confirmation
  resp <- ""
  while(resp != "y" && resp != "n" && resp != "Y" && resp != "N")
  {
    resp <- readline(prompt = "Proced with calibration of the whole MS image? (this may take some time) [y, n]:")
    if(resp != "y" && resp != "n" && resp != "Y" && resp != "N")
    {
      cat("Invalid response, valid responses are: y, n, Y and N. Try again.\n")
    }
  }
  
  if( resp == "n" || resp =="N")
  {
    cat("Calibration process aborted by user\n")
    return()
  }
  
  calImg$mass <- new_mass
  calImg$mean <- mIntensity #I do not care anymore if it is old or new data, just store new data
  
  #Store
  rMSI::SaveMsiData(calImg, output_fname)
}


#' calMzAxis.
#'
#' @param avgSpc_mz  The mass axis to calibrate.
#' @param ref_mz a vector of reference masses (for exaple the theorical gold peaks).
#' @param target_mz manually slected masses to be fittet to ref_masses (must be the same length than ref_mz).
#' @param method a string with the method used for interpolation, valid methods are loess and linear
#'
#' @return a list containing the calibrated mass axis and the interpolated mass error repect the original mass axis.
#' @export
#'
calMzAxis <- function(avgSpc_mz, ref_mz, target_mz, method = "loess" )
{
  a <- data.frame( targetMass = target_mz, refMass =  ref_mz)
  if(method == "linear")
  {
    fitmodel <- lm(refMass~targetMass, a)
  }
  else if (method == "loess")
  {
    fitmodel <- loess(refMass~targetMass, a, control = loess.control(surface = "direct"))
  }
  else
  {
    stop(paste("The specified method:", method, "is not a valid method\n"))
  }
  
  return( predict(fitmodel, newdata = data.frame( targetMass = testImg$mass)) )
}
