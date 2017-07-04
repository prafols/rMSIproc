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

FileBrowseWidget <- function( parent_widget, sLabel = "File:", dirSel = F, multiSel = T, fFilter = "txt", setdir_fun = NULL, getdir_fun = NULL )
{
  oldWarning<-options()$warn
  options(warn = -1)
  
  ## Get the environment for this
  ## instance of the function.
  this <- environment()
  
  ## Class data members 
  directorySelection <- dirSel
  fileFilter <- fFilter
  sPath <- c()
  MultiFileSelect <- multiSel
  GetWorkDir <- getdir_fun
  SetWorkDir <- setdir_fun
  
  ## Class methods
  BrowseDialog <- function(evt, ...)
  {
    if(this$directorySelection)
    {
      sTitle <- "Select a directory"
      sType <- "selectdir"
      lFilter <- list("All dirs" = list(patterns = c("*")))
    }
    else
    {
      sTitle <- "Select a file"
      sType <- "open"
      lFilter <- list(list(patterns = c(paste("*.", this$fileFilter, sep =""), paste("*.", toupper(this$fileFilter), sep =""))))
      names(lFilter) <- paste(as.character(gWidgets2::svalue(this$lbl)), "*.", this$fileFilter)
    }
    
    init_dir <- path.expand("~/")
    if(!is.null(this$GetWorkDir))
    {
      init_dir <- this$GetWorkDir()
    }
    mPath <- gWidgets2::gfile( text = sTitle, type = sType, multi = (!this$directorySelection) & this$MultiFileSelect, filter = lFilter, initial.dir = init_dir  ) 
    if( length( mPath) > 0)
    {
      this$sPath <- mPath
      if(this$MultiFileSelect && !this$directorySelection)
      {
        mPathTxt <- ""
        for( i in 1:length(this$sPath))
        {
          mPathTxt <- paste(mPathTxt, basename(this$sPath[i]),sep ="")
          if( i < length(this$sPath))
          {
            mPathTxt <- paste(mPathTxt, "\n", sep ="")
          }
        }
        gWidgets2::svalue(this$entry_files) <- mPathTxt
        RGtk2::gtkEntrySetText(gWidgets2::getToolkitWidget(this$entry_dir), dirname(this$sPath[1]))
      }
      else
      {
        RGtk2::gtkEntrySetText(gWidgets2::getToolkitWidget(this$entry_dir), this$sPath)
      }
    }
    if(!is.null(this$SetWorkDir))
    {
      if(dir.exists(gWidgets2::svalue(this$entry_dir)))
      {
        this$SetWorkDir( gWidgets2::svalue(this$entry_dir)) #It is already a dir so not using the dirname funciton
      }
      else
      {
        this$SetWorkDir( dirname(gWidgets2::svalue(this$entry_dir)) )
      }
    }
  }
  
  ##Public Methods
  SetDirSelection <- function( isDirSelection )
  {
    this$directorySelection <- isDirSelection
    gWidgets2::visible(this$entry_files) <- !isDirSelection
  }
  
  SetFileFilter <- function ( fFilter )
  {
    this$fileFilter <- fFilter
  }
  
  SetLabel <- function (sLabel)
  {
    gWidgets2::svalue(this$lbl) <- sLabel
  }
  
  GetPath <- function()
  {
    return(this$sPath)
  }
  
  SetEnabled <- function( enabled )
  {
    gWidgets2::enabled( this$frm ) <- enabled  
  }
  
  ClearPath <- function()
  {
    this$sPath <- c()
    RGtk2::gtkEntrySetText(gWidgets2::getToolkitWidget(this$entry_dir), "")
    gWidgets2::svalue(this$entry_files) <- ""
  }
    
  ## Class constructor
  boxTop <- gWidgets2::ggroup(horizontal = F, container = parent_widget)
  frm <- gWidgets2::gframe( container = boxTop)
  box <- gWidgets2::ggroup(horizontal = T, container = frm, expand = T, fill = T)
  box_lbl <- gWidgets2::ggroup(horizontal = F, container = box)
  lbl <- gWidgets2::glabel(sLabel, container = box_lbl)
  gWidgets2::size(lbl) <- c(110, -1)
  box_entry <- gWidgets2::ggroup(horizontal = F, container = box)
  entry_dir<-gWidgets2::gedit(width = 40, container = box_entry)
  gWidgets2::editable(entry_dir) <- F
  if( multiSel && !dirSel )
  {
    entry_files<-gWidgets2::gtext(width = 250, height = 150, container = box_entry)
    gWidgets2::editable(entry_files) <- F
  }
  box_btn <- gWidgets2::ggroup(horizontal = F, container = box)
  btn <- gWidgets2::gbutton("Browse", handler = this$BrowseDialog, container = box_btn )
  
  ## Set the name for the class
  class(this) <- append(class(this),"FileBrowseWidget")
  gc()
  
  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)
  
  return(this)
}
