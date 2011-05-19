TargetSearchGUI <- function() {
  envirGUI <- new.env()
  TS_GUI <- tktoplevel()
  tkwm.deiconify(TS_GUI)
  tkfocus(TS_GUI)
# Gets the GUI icon path and TargetSearch Version
  ts.ico     <- file.path(.find.package("TargetSearch"), "ico/Icon_TSGUI.ico")
  ts.version <- installed.packages()["TargetSearch", "Version"]
# Set the GUI's icon (Windows only)
  if (.Platform$OS.type == "windows")
      try(tkwm.iconbitmap(TS_GUI, ts.ico), silent=TRUE)
      
# Specify correct Path
  tktitle(TS_GUI) <- paste("TargetSearch", ts.version)
  tclServiceMode(FALSE)
# Section for Working Directory
  f1 <- tkframe(TS_GUI, borderwidth=2, relief="groove")
    f1z1 <- tkframe(f1, borderwidth=2)
    f1z2 <- tkframe(f1, borderwidth=2)
      labWD <- tklabel(parent=f1z1, text="Working Directory")
      valWD <- tclVar(getwd())
      ebxWD <- tkentry(parent=f1z2, textvariable=valWD)
      fncGetDir <- function() {
        newDir <- gsub("[{}]", "", tclvalue(tkchooseDirectory()))
        if (newDir != "") { tclvalue(valWD) <- newDir }
      }
      butWD <- tkbutton(parent=f1z2, text="Browse", width=7, command=fncGetDir)
      tkpack(labWD, anchor="w")
      tkpack(ebxWD, side="left", padx=1, expand="yes", fill="x") 
      tkpack(butWD, side="right", padx=1)
    tkpack(f1z1, f1z2, fill="x")
# Section for NetCDF Files
  f2 <- tkframe(TS_GUI, borderwidth=2, relief="groove")
    f2z1 <- tkframe(f2, borderwidth=2)
    f2z2 <- tkframe(f2, borderwidth=2)
      valFileNum <- tclVar("0")
      valFiles <- tclVar("")
      labData <- tklabel(parent=f2z1, text="File Import")
      tkconfigure(labData, "-text", paste("File Import (", tclvalue(valFileNum), " Files)", sep=""))
      rbValue <- tclVar("NetCDF Data")
      
      fncGrayParApex <- function(ApexVal) {
        if(ApexVal == "Apex Data") {
            tkconfigure(ebxBaseline, "-state", "disabled")
            tkconfigure(cbxBaseline, "-state", "disabled")
            tkconfigure(ebxPDthr, "-state", "disabled")
            tkconfigure(ebxPDmin, "-state", "disabled")
            tkconfigure(ebxPDmax, "-state", "disabled")
            tkconfigure(ebxPDSmoothing, "-state", "disabled")
            tkconfigure(butRILoad, "-state", "disabled")
            tkconfigure(butRICreate, "-state", "disabled")
            tkconfigure(butRIEdit, "-state", "disabled")
            tkconfigure(butRISave, "-state", "disabled")
        } else {
            tkconfigure(cbxBaseline, "-state", "normal")
            fncCBBaseline()
            tkconfigure(ebxPDthr, "-state", "normal")
            tkconfigure(ebxPDmin, "-state", "normal")
            tkconfigure(ebxPDmax, "-state", "normal")
            tkconfigure(ebxPDSmoothing, "-state", "normal")
            tkconfigure(butRILoad, "-state", "normal")
            tkconfigure(butRICreate, "-state", "normal")
            fncRIGrayEditSave()
        }
      }
      
      fncGetFiles <- function() {
        tmp.filetype <- ifelse(tclvalue(rbValue)=="Apex Data",
                               "{{Peak Apex Files} {.txt}} {{All files} *}",
                               "{{NetCDF Files} {.cdf}} {{All files} *}") 
        tclvalue(valFiles) <- tclvalue(
          tkgetOpenFile(initialdir = tclvalue(valWD), 
                        filetypes = tmp.filetype, 
                        multiple = TRUE))
        tclvalue(valFileNum) <- length(fncSplitTclStrg(tclvalue(valFiles)))
        tkconfigure(labData, "-text", paste("File Import (", tclvalue(valFileNum), " Files)", sep=""))
        
        fncGrayParApex(as.character(tclvalue(rbValue)))
        invisible()
      }
      rbtCDF <- tkradiobutton(parent=f2z2, variable=rbValue, command=fncGetFiles, value="NetCDF Data")
      rbtPeak <- tkradiobutton(parent=f2z2, variable=rbValue, command=fncGetFiles, value="Apex Data")
      fncSplitTclStrg <- function(x) {
        y <- vector(); i <- 1
        while (nchar(x)!=0) {
          if (substr(x, 1, 1)=="{") {
            j <- min(gregexpr("[}]", x)[[1]])
            y[i] <- substr(x, 2, j-1)
            x <- substr(x, j, nchar(x))
          } else {
            j <- min(gregexpr(" ", x)[[1]])
            if (j == -1) {
              y[i] <- x
              x <- ""
            } else {
              y[i] <- substr(x, 1, j-1)
              x <- substr(x, j, nchar(x))
            }
          }
          while(substr(x, 1, 1)=="}" | substr(x, 1, 1)==" ") {x <- substr(x, 2, nchar(x))}
          i <- i+1
        }
        return(y)
      }
      fncShowFiles <- function() {
        tmp.message <- fncSplitTclStrg(tclvalue(valFiles))
        tmp.number <- length(tmp.message)
        tmp.message <- paste(tmp.message, collapse="\n") 
        if (tmp.number == 0) { 
          tkmessageBox(title = "You selected 0 files.", message = "Please press 'Select' first.", icon = "info", type = "yes")
        }
        if (tmp.number > 50) {
          test <- tkmessageBox(title = "Question", message = "You selected >50 files. Show all?", icon = "question", type = "yesno")
          if (as.character(test)=="yes") {
            tkmessageBox(title = paste(tmp.number, "Files to be processed"), message = tmp.message, icon = "info", type = "ok")
          } else {
            #
          }
        } else {
          tkmessageBox(title = paste(tmp.number, "Files to be processed"), message = tmp.message, icon = "info", type = "ok")
        }
      }
      butShow <- tkbutton(parent=f2z2, text="Show", width=7, command=fncShowFiles)
      butData <- tkbutton(parent=f2z2, text="Select", width=7, command=fncGetFiles)
      tkpack(labData, anchor="w")
      tkpack(tklabel(f2z2, text="NetCDF Data"), rbtCDF, tklabel(f2z2, text="Apex Data"), rbtPeak, side="left", padx=1)  
      tkpack(butShow, butData, side="right", padx=1)
    tkpack(f2z1, f2z2, fill="x")
# Section for Parameters
  f3 <- tkframe(TS_GUI, borderwidth=2, relief="groove")
    # set number of parameter lines and main label
      #for (i in 1:11) { eval(substitute(x <- tkframe(f3, borderwidth=2), list(x=paste("f3z", i, sep="")))) }

      f3zBL <- tkframe(f3, borderwidth=2)  # Baseline
      f3zRI <- tkframe(f3, borderwidth=2)  # RI correction
      f3zPD1 <- tkframe(f3, borderwidth=2) # Peak Detection
      f3zPD2 <- tkframe(f3, borderwidth=2)
      f3zPD3 <- tkframe(f3, borderwidth=2)
      f3zLib <- tkframe(f3, borderwidth=2) # Library
      f3zLib2 <- tkframe(f3, borderwidth=2)
      f3zLib3 <- tkframe(f3, borderwidth=2)
      f3zLib4 <- tkframe(f3, borderwidth=2)
      f3z9 <- tkframe(f3, borderwidth=2)
      f3z10 <- tkframe(f3, borderwidth=2)
      f3z11 <- tkframe(f3, borderwidth=2)
      labPar <- tklabel(parent=f3, text="Processing Parameters")
    # Baseline
      labBaseline <- tklabel(parent=f3, text="    Baseline Correction")
      fncCBBaseline <- function() {
        if (tclvalue(valCBBaseline)==0) {
          tkconfigure(ebxBaseline, "-state", "disabled")
        } else {
          tkconfigure(ebxBaseline, "-state", "normal")
       }
        invisible()
      }
      valCBBaseline <- tclVar(0)
      cbxBaseline <- tkcheckbutton(parent=f3zBL, command=fncCBBaseline, variable=valCBBaseline)
      labBaselineOnOff <- tklabel(parent=f3zBL, text="on/off")
      valBaseline <- tclVar(0.5)
      ebxBaseline <- tkentry(parent=f3zBL, width=4, textvariable=valBaseline)
      tkconfigure(ebxBaseline, "-state", "disabled")
      labBaselineEBX <- tklabel(parent=f3zBL, text="threshold (0..1)")
    # Retention Index
      labRI <- tklabel(parent=f3, text="    Retention Index Corr.")
      valRIData <- tclVar("")
      fncRILoad <- function() {
        tmp <- tclvalue(tkgetOpenFile(initialdir = tclvalue(valWD), filetypes = "{{RI} {.txt}} {{All files} *}"))
        if (tmp != "") {
          tclvalue(valRIData) <- tmp 
          RI_SearchFrames <- read.delim(tclvalue(valRIData), sep="\t", header=T, check.names=F)
          TSPar <- get("TSPar", envir=envirGUI)
          TSPar$RI_SearchFrames <- RI_SearchFrames
          assign("TSPar", value = TSPar, envir = envirGUI)
          tkconfigure(butRIEdit, "-state", "normal")
          tkconfigure(butRISave, "-state", "normal")
        }
      }
      fncRICreate <- function() {
        tmp_RI <- edit(matrix("", nrow=1, ncol=3, dimnames=list(NULL, c("LowLimit", "HighLimit", "Standard"))))
        tmp_RI <- tmp_RI[apply(tmp_RI, 1, function(x) {all(x=="") == FALSE}), 1:3, drop=FALSE]
        if (prod(dim(tmp_RI))>0) {
          TSPar <- get("TSPar", envir = envirGUI)
          TSPar$RI_SearchFrames <- tmp_RI
          assign("TSPar", value = TSPar, envir = envirGUI)
          tkconfigure(butRIEdit, "-state", "normal")
          tkconfigure(butRISave, "-state", "normal")
        } else {
          tkconfigure(butRIEdit, "-state", "disabled")
          tkconfigure(butRISave, "-state", "disabled")
        }
      }
      fncRIGrayEditSave <- function() {
        TSPar <- get("TSPar", envir = envirGUI)
        if(all(TSPar$RI_SearchFrames == "")) {
          tkconfigure(butRIEdit, "-state", "disabled")
          tkconfigure(butRISave, "-state", "disabled")
        } else {
          tkconfigure(butRIEdit, "-state", "normal")
          tkconfigure(butRISave, "-state", "normal")
        }
      }
      fncRIEdit <- function() {
        TSPar <- get("TSPar", envir = envirGUI)
        TSPar$RI_SearchFrames <- edit(TSPar$RI_SearchFrames)
        assign("TSPar", value = TSPar, envir = envirGUI)
      }
      fncRISave <- function() {
        tmp <- tclvalue(tkgetSaveFile(initialdir = tclvalue(valWD), filetypes = "{{RI} {.txt}} {{All files} *}"))
        TSPar <- get("TSPar", envir=envirGUI)
        if(tmp!="") { 
          tclvalue(valRIData) <- tmp 
          write.table(TSPar$RI_SearchFrames, file=tclvalue(valRIData), sep="\t", col.names=TRUE, row.names=TRUE, quote = FALSE)
        }
      }
      butRILoad <- tkbutton(parent=f3zRI, text="Load", command=fncRILoad, width=7)
      butRICreate <- tkbutton(parent=f3zRI, text="Create", command=fncRICreate, width=7)
      butRIEdit <- tkbutton(parent=f3zRI, text="Edit/View", command=fncRIEdit, width=7)
      butRISave <- tkbutton(parent=f3zRI, text="Save", command=fncRISave, width=7)
    # Peak Detection
      labPD <- tklabel(parent=f3, text="    Peak Detection")
      valPDthr <- tclVar(10)
      ebxPDthr <- tkentry(parent=f3zPD1, width=4, textvariable=valPDthr)
      labPDthr <- tklabel(parent=f3zPD1, text="threshold [Intensity Counts]")
      labPD2 <- tklabel(parent=f3, text="")
      valPDmin <- tclVar(85); valPDmax <- tclVar(500)
      ebxPDmin <- tkentry(parent=f3zPD2, width=4, textvariable=valPDmin)
      ebxPDmax <- tkentry(parent=f3zPD2, width=4, textvariable=valPDmax)
      labPDmass <- tklabel(parent=f3zPD2, text="mass range (min/max))")
      labPD3 <- tklabel(parent=f3, text="")
      valPDSmoothing <- tclVar(7)
      ebxPDSmoothing <- tkentry(parent=f3zPD3, width=4, textvariable=valPDSmoothing)
      labPDSmoothing <- tklabel(parent=f3zPD3, text="smoothing (Data Points)")

    # End of peak detection
    
    # Library
      labLib <- tklabel(parent=f3, text="    Library")
      valLibData <- tclVar("")
      fncLibLoad <- function() {
        tmp <- tclvalue(tkgetOpenFile(initialdir = tclvalue(valWD), filetypes = "{{Lib} {.txt}} {{All files} *}"))
        if (tmp != "") {
          tclvalue(valLibData) <- tmp 
          Lib_Data <- read.delim(tclvalue(valLibData), sep="\t", header=T, quote="", comment="", as.is=T, check.names=F)
          TSPar <- get("TSPar", envir=envirGUI)
          TSPar$Library_Data <- Lib_Data
          assign("TSPar", value = TSPar, envir = envirGUI)
          tkconfigure(butLibEdit, "-state", "normal")
          tkconfigure(butLibSave, "-state", "normal")
        }
      }
      fncLibCreate <- function() {
        tmp_Lib <- edit(matrix("", nrow=1, ncol=3, dimnames=list(NULL, c("Name", "RI", "SEL_MASS"))))
        tmp_Lib <- tmp_Lib[apply(tmp_Lib, 1, function(x) {all(x=="") == FALSE}), 1:3, drop=FALSE]
        if (prod(dim(tmp_Lib))>0) {
          TSPar <- get("TSPar", envir = envirGUI)
          TSPar$Library_Data <- tmp_Lib
          assign("TSPar", value = TSPar, envir = envirGUI)
          tkconfigure(butLibEdit, "-state", "normal")
          tkconfigure(butLibSave, "-state", "normal")
        } else {
          tkconfigure(butLibEdit, "-state", "disabled")
          tkconfigure(butLibSave, "-state", "disabled")
        }
      }
      fncLibGrayEditSave <- function() {
        TSPar <- get("TSPar", envir = envirGUI)
        if (all(TSPar$"Library_Data" == "")) {
          tkconfigure(butLibEdit, "-state", "disabled")
          tkconfigure(butLibSave, "-state", "disabled")
        } else {
          tkconfigure(butLibEdit, "-state", "normal")
          tkconfigure(butLibSave, "-state", "normal")
        }
      }
      fncLibEdit <- function() {
        TSPar <- get("TSPar", envir = envirGUI)
        TSPar$Library_Data <- edit(TSPar$Library_Data)
        assign("TSPar", value = TSPar, envir = envirGUI)
      }
      fncLibSave <- function() {
        tmp <- tclvalue(tkgetSaveFile(initialdir = tclvalue(valWD), filetypes = "{{Lib} {.txt}} {{All files} *}"))
        TSPar <- get("TSPar", envir=envirGUI)
        if(tmp!="") { 
          tclvalue(valLibData) <- tmp 
          write.table(TSPar$Library_Data, file=tclvalue(valLibData), sep="\t", col.names=TRUE, row.names=TRUE, quote = FALSE)
        }
      }
      butLibLoad <- tkbutton(parent=f3zLib, text="Load", command=fncLibLoad, width=7)
      butLibCreate <- tkbutton(parent=f3zLib, text="Create", command=fncLibCreate, width=7)
      butLibEdit <- tkbutton(parent=f3zLib, text="Edit/View", command=fncLibEdit, width=7)
      butLibSave <- tkbutton(parent=f3zLib, text="Save", command=fncLibSave, width=7)
      
    # Library Search Options
      labLib2 <- tklabel(parent=f3, text="")
      valPDW1 <- tclVar(2000); valPDW2 <- tclVar(1000); valPDW3 <- tclVar(500)
      ebxPDW1 <- tkentry(parent=f3zLib2, width=4, textvariable=valPDW1)
      ebxPDW2 <- tkentry(parent=f3zLib2, width=4, textvariable=valPDW2)
      ebxPDW3 <- tkentry(parent=f3zLib2, width=4, textvariable=valPDW3)
      labPDW <- tklabel(parent=f3zLib2, text="Search Windows")

    # top Massses section and exclude Masses
      labLib3 <- tklabel(parent=f3, text="")
      valPDtopMass <- tclVar(10)
      ebxPDtopMass <- tkentry(parent=f3zLib3, width=4, textvariable=valPDtopMass)
      labPDtopMass <- tklabel(parent=f3zLib3, text="# of top masses")
      labLib4 <- tklabel(parent=f3, text="")
      valPDexcludeMass <- tclVar("147:149")
      ebxPDexcludeMass <- tkentry(parent=f3zLib4, width=20, textvariable=valPDexcludeMass)
      labPDexcludeMass <- tklabel(parent=f3zLib4, text="excluded masses")


    # Normalization
      labNorm <- tklabel(parent=f3, text="    Normalization")
      rbNorm <- tclVar("dayNorm")
      labNorm1 <- tklabel(f3z9, text="None")
      labNorm2 <- tklabel(f3z9, text="DayNorm")
      labNorm3 <- tklabel(f3z9, text="MedNorm")
      rbtNorm1 <- tkradiobutton(f3z9, variable=rbNorm, value="none")
      rbtNorm2 <- tkradiobutton(f3z9, variable=rbNorm, value="dayNorm")
      rbtNorm3 <- tkradiobutton(f3z9, variable=rbNorm, value="medianNorm")
    # Profiles
      labProf <- tklabel(parent=f3, text="    Final Profiles")
      valProfTS <- tclVar("500"); valProfthr <- tclVar("0.95"); valProfSNr <- tclVar("6")
      ebxProfTS <- tkentry(parent=f3z10, width=4, textvariable=valProfTS)
      labProfTS <- tklabel(parent=f3z10, text="timesplit")
      ebxProfthr <- tkentry(parent=f3z10, width=4, textvariable=valProfthr)
      labProfthr <- tklabel(parent=f3z10, text="correl. thr.")
      labProfz11 <- tklabel(parent=f3, text="")
      ebxProfSNr <- tkentry(parent=f3z11, width=4, textvariable=valProfSNr)
      labProfSNr <- tklabel(parent=f3z11, text="min number of correl. samples")
      # Pack Widgets
      tkpack(cbxBaseline, side="left", padx=0, anchor="w") 
      tkpack(labBaselineOnOff, ebxBaseline, labBaselineEBX, side="left", padx=1, anchor="w")
      tkpack(butRILoad, butRICreate, butRIEdit, butRISave, side="left", padx=1, anchor="w")
      tkpack(ebxPDW1, ebxPDW2, ebxPDW3, labPDW, side="left", anchor="w", padx=1)
      tkpack(ebxPDthr, labPDthr, side="left", anchor="w", padx=1)
      tkpack(ebxPDmin, ebxPDmax, labPDmass, side="left", anchor="w", padx=1)
      tkpack(ebxPDSmoothing, labPDSmoothing, side="left", anchor="w", padx=1)
      tkpack(ebxPDtopMass, labPDtopMass, side="left", anchor="w", padx=1)      
      tkpack(ebxPDexcludeMass, labPDexcludeMass, side="left", anchor="w", padx=1)            
      tkpack(butLibLoad, butLibCreate, butLibEdit, butLibSave, side="left", padx=1, anchor="w")
      tkpack(labNorm1, rbtNorm1, labNorm2, rbtNorm2, labNorm3, rbtNorm3, side="left", padx=1, anchor="w")
      tkpack(ebxProfTS, labProfTS, ebxProfthr, labProfthr, side="left", padx=1, anchor="w") 
      tkpack(ebxProfSNr, labProfSNr, side="left", padx=1, anchor="w")
      tkgrid(labPar, sticky="w")
      tkgrid(labBaseline, f3zBL, sticky="w")
      tkgrid(labPD, f3zPD1, sticky="w")
      tkgrid(labPD2, f3zPD2, sticky="w")
      tkgrid(labPD3, f3zPD3, sticky="w")
      tkgrid(labRI, f3zRI, sticky="w")
      tkgrid(labLib, f3zLib, sticky="w")
      tkgrid(labLib2, f3zLib2, sticky="w")
      tkgrid(labLib3, f3zLib3, sticky="w")
      tkgrid(labLib4, f3zLib4, sticky="w")

      tkgrid(labNorm, f3z9, sticky="w")
      tkgrid(labProf, f3z10, sticky="w")
      tkgrid(labProfz11, f3z11, sticky="w")
# Section for Working Directory
  f4 <- tkframe(TS_GUI, borderwidth=2, relief="groove")
    f4z1 <- tkframe(f4, borderwidth=2)
    f4z2 <- tkframe(f4, borderwidth=2)
    labF41 <- tklabel(f4z1, text="- Parameters -", width=20)
    labF42 <- tklabel(f4z1, text="- Program -", width=20) 
      valMainData <- tclVar("")
      fncPutParameters <- function(TSPar) {
        # write Parameters from TSPar (R List) into TargetSearch GUI
        tclvalue(valWD) <- TSPar$"WorkingDirectory"
        tclvalue(valFileNum) <- TSPar$"Files"$"FileNum"
        tclvalue(rbValue) <- TSPar$"Files"$"FileType"
        fncGrayParApex(TSPar$"Files"$"FileType")
        tclvalue(valFiles) <- TSPar$"Files"$"FilePath"
          tkconfigure(labData, "-text", paste("File Import (", tclvalue(valFileNum), " Files)", sep=""))
        tclvalue(valCBBaseline) <- TSPar$"Baseline"$"On"
          fncCBBaseline()
        tclvalue(valBaseline) <- TSPar$"Baseline"$"Threshold"
        tclvalue(valRIData) <- TSPar$"RI"
        tclvalue(valPDW1) <- TSPar$"PeakDetection"$"SearchWin"[1]
        tclvalue(valPDW2) <- TSPar$"PeakDetection"$"SearchWin"[2]
        tclvalue(valPDW3) <- TSPar$"PeakDetection"$"SearchWin"[3]
        tclvalue(valPDthr) <- TSPar$"PeakDetection"$"Threshold"
        tclvalue(valPDmin) <- TSPar$"PeakDetection"$"MassRange"[1]
        tclvalue(valPDmax) <- TSPar$"PeakDetection"$"MassRange"[2]
        tclvalue(valPDSmoothing) <- TSPar$"PeakDetection"$"Smoothing"
        tclvalue(valPDtopMass)     <- TSPar$"PeakDetection"$"topMass"
        tclvalue(valPDexcludeMass) <- TSPar$"PeakDetection"$"excludeMass"
        tclvalue(valLibData) <- TSPar$"Library"
        tclvalue(rbNorm) <- TSPar$"Normalization"
        tclvalue(valProfTS) <- TSPar$"Profiles"$"TimeSplit"
        tclvalue(valProfthr) <- TSPar$"Profiles"$"Threshold"
        tclvalue(valProfSNr) <- TSPar$"Profiles"$"MinSamNr"
      }
      fncInitialParameters <- function() {
        # write all Parameters into a R list
        TSPar <- list(
          "WorkingDirectory" = tclvalue(valWD),
          "Files" = list(
            "FileNum" = as.numeric(tclvalue(valFileNum)),
            "FileType" = tclvalue(rbValue),
            "FilePath" = fncSplitTclStrg(tclvalue(valFiles))
          ),
          "Baseline" = list(
            "On" = as.logical(as.numeric(tclvalue(valCBBaseline))),
            "Threshold" = as.numeric(tclvalue(valBaseline))
          ),
          "RI" = tclvalue(valRIData),
          "RI_SearchFrames" = "",
          "PeakDetection" = list(
            "SearchWin" = as.numeric(c(tclvalue(valPDW1), tclvalue(valPDW2), tclvalue(valPDW3))),
            "Threshold" = as.numeric(tclvalue(valPDthr)),
            "MassRange" = as.numeric(c(tclvalue(valPDmin), tclvalue(valPDmax))),
            "Smoothing" = as.numeric(tclvalue(valPDSmoothing)),
            "topMass"   = as.numeric(tclvalue(valPDtopMass)),
            "excludeMass" = tclvalue(valPDexcludeMass)
          ),
          "Library" = tclvalue(valLibData),
          "Library_Data" = "",
          "Normalization" = tclvalue(rbNorm),
          "Profiles" = list(
            "TimeSplit" = as.numeric(tclvalue(valProfTS)),
            "Threshold" = as.numeric(tclvalue(valProfthr)),
            "MinSamNr" = as.numeric(tclvalue(valProfSNr))
          )
        )
        return(TSPar)
      }
      fncMainLoad <- function() {
        tclvalue(valMainData) <- tclvalue(tkgetOpenFile(title = "Please select a TargetSearch Parameter File (*.RData)",
          initialdir = tclvalue(valWD), filetypes = "{{Par} {.RData}} {{All files} *}"))
        load(tclvalue(valMainData))
        if (grep("TSPar", ls(envir=envirGUI))==FALSE) {
          TSPar <- fncInitialParameters()
        } else {
          assign("TSPar", TSPar, envir=envirGUI)
        }
        fncRIGrayEditSave()
        fncLibGrayEditSave()
        fncPutParameters(TSPar)
      }
      fncMainRun <- function() {
      #
      # TODO: Check if all Parameters are meaningfull
      #
        UpdateTextGUI <- function(l="- Parameters -", r="- Program -") {
          if(is.na(l) == FALSE) tkconfigure(labF41, "-text", l)
          if(is.na(r) == FALSE) tkconfigure(labF42, "-text", r) 
          .Tcl("update idletasks")
        }
        UpdateTextGUI("Running... (Please wait!)", "Reading Parameters")
        cat("\nTargetSearch", ts.version, "is processing your data...\nPlease wait until finished.\n\n")
        
      # A function to read a masses range to an R-vector
      # Example:  fcnMzRangeToVector("126,147:149,160") = 126, 147, 148, 149, 160.
      fcnMzRangeToVector <- function(x) {
        x <- gsub(" ", "",x[1]) # remove spaces
        if(grepl("^(\\d+[,:]?)+$", x, perl=TRUE)) {
            x <- sub("[,:]$","",x)
            x <- strsplit(unlist(strsplit(x, ",")), ":")
            x <- lapply(x, as.numeric)
            unlist(lapply(x, function(z) if(length(z)==1) z else seq(z[1],z[2])))
        } else {
            return(NULL)
        }
      }

        # Start 'RUN'
        CDFdir <- tclvalue(valWD)
        RIdir  <- tclvalue(valWD)
        setwd(CDFdir)
        if (as.character(tclvalue(rbValue)) == "Apex Data") {
          RI_FILE <- basename( fncSplitTclStrg(tclvalue(valFiles)) )
          RIdir   <- dirname( fncSplitTclStrg(tclvalue(valFiles)) )[1]
          CDF_FILE <- gsub("RI_", "", gsub("txt", "cdf", RI_FILE))
        }
        if (as.character(tclvalue(rbValue)) == "NetCDF Data") {
          CDF_FILE <- basename( fncSplitTclStrg(tclvalue(valFiles)) )
          CDFdir   <- dirname( fncSplitTclStrg(tclvalue(valFiles)) )[1]
          RI_FILE <- gsub("cdf", "txt", paste("RI_", CDF_FILE, sep=""))
        }
        MEASUREMENT_DAY <- sub('^([0-9]+).*$','\\1', CDF_FILE)
        UpdateTextGUI(NA, "Import Samples")
        samples   <- new("tsSample", CDFfiles = CDF_FILE, RIfiles = RI_FILE, days = MEASUREMENT_DAY, CDFpath = CDFdir, RIpath = RIdir)
        if (all(get("TSPar", envir=envirGUI)$Library_Data == "")) {
          Lib <- ImportLibrary(libfile = tclvalue(valLibData),
                               RI_dev = as.numeric(c(tclvalue(valPDW1), tclvalue(valPDW2), tclvalue(valPDW3))))
        } else {
          Lib <- as.data.frame(get("TSPar", envir=envirGUI)$Library_Data, stringsAsFactors=FALSE)
##          print(tclvalue(valPDtopMass))
##          print(tclvalue(valPDexcludeMass))
          # insert the options from the GUI
          tM  <- as.numeric(tclvalue(valPDtopMass))
          eM  <- fcnMzRangeToVector(tclvalue(valPDexcludeMass))
          Lib <- ImportLibrary.tab(libdata=Lib, RI_dev=as.numeric(c(tclvalue(valPDW1), tclvalue(valPDW2), tclvalue(valPDW3))),
                                   TopMasses=tM, ExcludeMasses=eM)
        }
      # RI correction can be ommited in case that RI_..txt Files are already existent
        if (as.character(tclvalue(rbValue)) == "NetCDF Data") {
          UpdateTextGUI(NA, "Detect RI")
          if (tclvalue(valRIData) == "") {
            tkmessageBox(title = "RI-Marker File is missing...", 
                         message = "Please Load/Create a file containing appropriate search windows for RI-Marker first.", 
                         icon = "error", type = "ok")
                         UpdateTextGUI()
                         stop("Please Start Again.")
          }
          rimLimits <- ImportFameSettings(tmp.file = tclvalue(valRIData), mass = NA)
          RImatrix <- RIcorrect(samples = samples, rimLimits = rimLimits,
                                massRange = as.numeric(c(tclvalue(valPDmin), tclvalue(valPDmax))),
                                Window = as.numeric(tclvalue(valPDSmoothing)), 
                                IntThreshold = as.numeric(tclvalue(valPDthr)),
                                baseline = as.logical(as.numeric(tclvalue(valCBBaseline))),
                                baseline.opts = list(threshold = as.numeric(tclvalue(valBaseline))),
                                showProgressBar = TRUE)
          outliers <- FAMEoutliers(samples, RImatrix)
        }
      # search for median RI and updates
        UpdateTextGUI(NA, "Find Peaks")
        Lib    <- medianRILib(samples, Lib)
        cor_RI <- sampleRI(samples = samples, Lib = Lib,
                           r_thres = as.numeric(tclvalue(valProfthr)), 
                           method = tclvalue(rbNorm), 
                           minPairObs = as.numeric(tclvalue(valProfSNr)))
        UpdateTextGUI(NA, "Refine Peak Search")
        peakData <- peakFind(samples = samples, Lib = Lib, cor_RI = cor_RI)
        UpdateTextGUI(NA, "Finishing...")
        MetabProfile <- Profile(samples = samples, Lib = Lib, peakData = peakData, 
                                r_thres = as.numeric(tclvalue(valProfthr)), 
                                method = tclvalue(rbNorm), 
                                minPairObs = as.numeric(tclvalue(valProfSNr)))
        finalProfile <- ProfileCleanUp(Profile = MetabProfile, 
                                       timeSplit = as.numeric(tclvalue(valProfTS)), 
                                       r_thres = as.numeric(tclvalue(valProfthr)),
                                       minPairObs = as.numeric(tclvalue(valProfSNr)))
      # save your results in tabbed text format files
        Write.Results(Lib, finalProfile)
        
      # save workspace
        if (as.character(tclvalue(rbValue)) == "NetCDF Data") {
            save(samples, rimLimits, RImatrix, outliers, Lib, cor_RI, peakData, MetabProfile, finalProfile,
                file=paste("TargetSearch-", Sys.Date(), ".RData", sep = ""))
        } else {
            save(samples, Lib, cor_RI, peakData, MetabProfile, finalProfile,
                file=paste("TargetSearch-", Sys.Date(), ".RData", sep = ""))
        }
        UpdateTextGUI()
        cat(paste("\nTargetSearch ", ts.version, " finished processing your data...\nOutput was stored to ", tclvalue(valWD),"\n\n", sep=""))
      }
      fncMainQuit <- function() {
        tkdestroy(TS_GUI)
      }
      fncMainSave <- function() {
        TSPar <- fncInitialParameters()
        TSPar$RI_SearchFrames <- get("TSPar", envir=envirGUI)$RI_SearchFrames
        TSPar$Library_Data <- get("TSPar", envir=envirGUI)$Library_Data
        tmp <- tclvalue(tkgetSaveFile(title = "Save as TargetSearch Parameter File (*.RData)",
          initialdir = tclvalue(valWD), filetypes = "{{Par} {.RData}} {{All files} *}"))
        if(tmp!="") { save(TSPar, file=tmp) }
      }
    butMainLoad <- tkbutton(parent=f4z2, text="Load", width=9, command=fncMainLoad)
    butMainSave <- tkbutton(parent=f4z2, text="Save", width=9, command=fncMainSave)
    butMainRun <- tkbutton(parent=f4z2, text="Run", width=9, command=fncMainRun)
    butMainQuit <- tkbutton(parent=f4z2, text="Quit", width=9, command=fncMainQuit)
    tkpack(labF41, labF42, side="left", padx=8, pady=0, expand="yes", fill="x")
    tkpack(butMainLoad, butMainSave, butMainRun, butMainQuit, side="left", padx=4, pady=4, expand="yes", fill="x")
   tkpack(f4z1, f4z2) 
# Finalize It  
  tkpack(f1, f2, f3, f4, fill="x")
  tkconfigure(butRIEdit, "-state", "disabled")
  tkconfigure(butRISave, "-state", "disabled")
  tkconfigure(butLibEdit, "-state", "disabled")
  tkconfigure(butLibSave, "-state", "disabled")
  tclServiceMode(TRUE)
  assign("TSPar", fncInitialParameters(), envir=envirGUI)
invisible()
}  
