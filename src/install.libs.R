#Custom install script to allow coping the fftw dll in windows system

files <- Sys.glob(paste0("*", SHLIB_EXT))
dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)
if(file.exists("symbols.rds"))
{
  file.copy("symbols.rds", dest, overwrite = TRUE)
}

if(WINDOWS)
{
  if(R_ARCH == "/x64")
  {
    cat("64bits windows OS detected. Copying fftw3 dll...\n")
    fftfiles <-file.path(R_PACKAGE_SOURCE, "src/fftw3_win/lib64/libfftw3-3.dll")  
  }
  else
  {
    cat("32bits windows OS detected. Copying fftw3 dll...\n")
    fftfiles <-file.path(R_PACKAGE_SOURCE, "src/fftw3_win/lib32/libfftw3-3.dll")  
  }
  file.copy(fftfiles, dest, overwrite = TRUE)
}

