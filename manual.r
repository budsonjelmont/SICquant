#!/usr/bin/Rscript
# @aly

# To run, launch R gui and type in source("PATH_TO_SCRIPT_\\manual.r", chdir=T)
# xcms2.r must be in same folder
# TODO: 
# Zoom in and out at specific rt windows
# Adjust peak left and right when it would calculate area, from default peak detection
# Slider?
# Drag and drop?
# Fix arrangement of gui

# source("C:\\Users\\superuser.SALOMON\\Documents\\xmlparse\\manual.r", chdir=T)

# http://mcu.edu.tw/~chenmh/teaching/project/r/reference/RTclTkExamples/FileOpenSave.html
# http://stackoverflow.com/questions/15272916/how-to-wait-for-a-keypress-in-r/15283106
# https://stat.ethz.ch/pipermail/r-help/2008-November/179241.html
# http://www.sciviews.org/_rgui/tcltk/Events.html
# http://www.sciviews.org/_rgui/tcltk/InteractiveTkrPlot.html

suppressPackageStartupMessages({
    library(tcltk)
    library(tkrplot)
    library(xcms)
    })

xcmsScript <<- "xcms2_manual.r"
rtwl <<- c(1440, 720, 360, 180, 90, 45, 15, 5)
zoomlevel <<- 1
########################################################################################
# Input: Nothing
# Output: Directory of running script
# 
# This method looks at the argument, tries to get the path of the current running 
#  script, and return it.
# http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
########################################################################################
currentDir <- function() {
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.basename <- gsub("\\\\", "/", dirname(script.name))
    return(script.basename)
}

plotEICr <- function(object,mzrange = numeric(),rtrange = numeric(),scanrange = numeric()) {
    EIC <- rawEIC(object,mzrange=mzrange, rtrange=rtrange, scanrange=scanrange)
    points <- cbind(object@scantime[EIC$scan], EIC$intensity)
    plot(points, type="l", xlab="Seconds",ylab="Intensity")

    if (!is.na(as.numeric(tclvalue(peakrtmin))) && !is.na(as.numeric(tclvalue(peakrtmax)))) {
        lines(abline(v=as.numeric(tclvalue(peakrtmin)), col='blue'))
        lines(abline(v=as.numeric(tclvalue(peakrtmax)), col='blue'))
    }

    mzpoints <- intersect(which(xr@msnPrecursorMz >= mzrange[1]), which(xr@msnPrecursorMz <= mzrange[2]))
    if (length(mzpoints > 0)) {
        rt <- xr@scantime[xr@msnPrecursorScan[mzpoints]]
        rtpoints <- intersect(which(rt >= rtrange[1]), which(rt <= rtrange[2]))
        if (length(rtpoints > 0)) {
            lines(abline(v=rt[rtpoints], col='red', lty=2))
        }
    }
}

doplot <- function(...) {
    if (tclvalue(mass) != "" && tclvalue(ppm) != "" && tclvalue(rt) != "" && tclvalue(rtw) != "" && !is.null(xr)) {
        tryCatch({
            mzrange <- pmMz(as.numeric(tclvalue(mass)), as.numeric(tclvalue(ppm)))
            rtrange <- pmRt(as.numeric(tclvalue(rt)), as.numeric(tclvalue(rtw)))
            xrrtrange <- range(xr@scantime)
            if (rtrange[1] < xrrtrange[1]) rtrange[1] <- xrrtrange[1]
            if (rtrange[2] > xrrtrange[2]) rtrange[2] <- xrrtrange[2]

            plotEICr(xr, mzrange=mzrange, rtrange=rtrange)
            },error=function(ex) {
                plot(1, type="n", axes=F, xlab="", ylab="")
            })

        if (!is.null(manualLeft)) abline(v=manualLeft)
        if (!is.null(manualRight)) abline(v=manualRight)
    } else {
        plot(1, type="n", axes=F, xlab="", ylab="")
    }

    parPlotSize <<- par("plt")
    usrCoords <<- par("usr")
}

resethandler <- function() {
    tclvalue(mass) <- ""
    tclvalue(ppm) <- 10

    tclvalue(peakrt) <- ""
    tclvalue(peakrtmin) <- ""
    tclvalue(peakrtmax) <- ""
    tclvalue(peakint) <- ""
    tclvalue(peaksn) <- ""

    zoom <- 1
    assign("zoomlevel", zoom, envir=.GlobalEnv)

    tclvalue(rtw) <- rtwl[zoom]

    if (!is.null(xr)) tclvalue(rt) <- round(mean(range(xr@scantime)))
    tkrreplot(img)

    manualLeft <<- NULL
    manualRight <<- NULL
}

peakhandler <- function() {
    if (tclvalue(mass) != "" && tclvalue(ppm) != "" && tclvalue(rt) != "" && tclvalue(rtw) != "" && !is.null(xr)) {
        manualLeft <<- NULL
        manualRight <<- NULL
        p <- peakCalc(xr, as.numeric(tclvalue(mass)), as.numeric(tclvalue(rt)), as.numeric(tclvalue(ppm)), NULL, FALSE, FALSE, as.numeric(tclvalue(rtw)))
        if (!is.na(p[["peaks"]][[1]][["peakArea"]])) {
            tclvalue(peakrt) <- p[["peaks"]][[1]][["peakApexRT"]]
            tclvalue(peakrtmin) <- p[["peaks"]][[1]][["peakLeftRT"]]
            tclvalue(peakrtmax) <- p[["peaks"]][[1]][["peakRightRT"]]
            tclvalue(peakint) <- p[["peaks"]][[1]][["peakArea"]]
            # tclvalue(peakint) <- p[["peaks"]][[1]][["peakArea"]]
            tclvalue(peaksn) <- p[["peaks"]][[1]][["sn"]]
        } else {
            tclvalue(peakrt) <- ""
            tclvalue(peakrtmin) <- ""
            tclvalue(peakrtmax) <- ""
            tclvalue(peakint) <- ""
            tclvalue(peaksn) <- ""
        }
    }
    tkrreplot(img)
    update()
}

rawfilehandler <- function() {
    mzXMLPath <- tclvalue(tkgetOpenFile(filetypes="{{mzXML Files} {.mzXML}} {{All files} *}"))
    if (mzXMLPath == "") {
        return
    }
    # xrtmp <- tryCatch(xcmsRaw(name, profstep=0), error=function(ex) NULL)
    xr <<- tryCatch(xcmsRaw(mzXMLPath, profstep=0, includeMSn=TRUE), error=function(ex) NULL)
    # assign("xr", xrtmp, envir=.GlobalEnv)
    # ms2t <- tryCatch(msn2xcmsRaw(xcmsRaw(name, profstep=0, includeMSn=TRUE))@scantime, error=function(ex) NULL)
    if (!is.null(xr)) {
        # ms2spec <- tryCatch(msn2xcmsRaw(xcmsRaw(name, profstep=0, includeMSn=TRUE)), error=function(ex) NULL)
        # ms2t <- findPeaks(ms2spec, method="MS1")
        # assign("ms2", ms2t, envir=.GlobalEnv)
        # assign("ms2", ms2spec, envir=.GlobalEnv)
        resethandler()
        tkwm.title(base,mzXMLPath)
        tclvalue(rawpathvar) <- mzXMLPath
        # print(range(xr@scantime))
        # print(range(xr@env$mz))
    }
}

zoomin <- function() {
    if (tclvalue(mass) != "" && tclvalue(ppm) != "" && tclvalue(rt) != "" && !is.null(xr)) {
        if (zoomlevel < length(rtwl)) {
            zoom <- zoomlevel + 1
            assign("zoomlevel", zoom, envir=.GlobalEnv)
            rtrange <- pmRt(as.numeric(tclvalue(rt)), rtwl[zoom])
            xrrtrange <- range(xr@scantime)
            if (rtrange[1] < xrrtrange[1] + rtwl[zoom]) {
                tclvalue(rt) <- xrrtrange[1] + rtwl[zoom]
            } else if (rtrange[2] > xrrtrange[2] - rtwl[zoom]) {
                tclvalue(rt) <- xrrtrange[2] - rtwl[zoom]
            }
            tclvalue(rtw) <- rtwl[zoom]
            tkrreplot(img)
        }
    }
}

zoomout <- function() {
    if (tclvalue(mass) != "" && tclvalue(ppm) != "" && tclvalue(rt) != "" && !is.null(xr)) {
        if (zoomlevel > 1) {
            zoom <- zoomlevel - 1
            assign("zoomlevel", zoom, envir=.GlobalEnv)
            rtrange <- pmRt(as.numeric(tclvalue(rt)), rtwl[zoom])
            xrrtrange <- range(xr@scantime)
            if (rtrange[1] < xrrtrange[1] + rtwl[zoom]) {
                tclvalue(rt) <- xrrtrange[1] + rtwl[zoom]
            } else if (rtrange[2] > xrrtrange[2] - rtwl[zoom]) {
                tclvalue(rt) <- xrrtrange[2] - rtwl[zoom]
            }
            tclvalue(rtw) <- rtwl[zoom]
            tkrreplot(img)
        }
    }
}

leftshift <- function() {
    if (tclvalue(mass) != "" && tclvalue(ppm) != "" && tclvalue(rt) != "" && !is.null(xr)) {
        rtrange <- range(xr@scantime)
        if (as.numeric(tclvalue(rt)) - 0.25 * rtwl[zoomlevel] <= rtrange[1] + 0.25 * rtwl[zoomlevel]) {
            tclvalue(rt) <- rtrange[1] + 0.25 * rtwl[zoomlevel]
        } else {
            tclvalue(rt) <- as.numeric(tclvalue(rt)) - 0.25 * rtwl[zoomlevel]
        }
        tkrreplot(img)
    }
}

rightshift <- function() {
    if (tclvalue(mass) != "" && tclvalue(ppm) != "" && tclvalue(rt) != "" && !is.null(xr)) {
        rtrange <- range(xr@scantime)
        if (as.numeric(tclvalue(rt)) + 0.25 * rtwl[zoomlevel] >= rtrange[2] - 0.25 * rtwl[zoomlevel]) {
            tclvalue(rt) <- rtrange[2] - 0.25 * rtwl[zoomlevel]
        } else {
            tclvalue(rt) <- as.numeric(tclvalue(rt)) + 0.25 * rtwl[zoomlevel]
        }
        tkrreplot(img)
    }
}

# loadInput <- function() {
#     args <- commandArgs(trailingOnly = TRUE)
#     fminput <- file(args[1], open="r", encoding="UTF-16LE")
#     lines <- readLines(fminput, warn=FALSE)[1]
#     close(fminput)

#     tmpArr <- strsplit(lines, "\\|")[[1]]
#     rawFilePath <- tmpArr[6]
#     fileNameNoExt <- sub("^([^.]*).*", "\\1", basename(rawFilePath))
#     mzXMLPath <- paste(dirname(rawFilePath), paste(fileNameNoExt, ".mzXML", sep=""), sep="/")

#     if (file.exists(mzXMLPath)) {
#         xrtmp <- tryCatch(xcmsRaw(mzXMLPath, profstep=0), error=function(ex) NULL)
#         assign("xr", xrtmp, envir=.GlobalEnv)
#         if (!is.null(xr)) {
#             ms2spec <- tryCatch(msn2xcmsRaw(xcmsRaw(mzXMLPath, profstep=0, includeMSn=TRUE)), error=function(ex) NULL)
#             # ms2t <- findPeaks(ms2spec, method="MS1")
#             # assign("ms2", ms2t, envir=.GlobalEnv)
#             assign("ms2", ms2spec, envir=.GlobalEnv)
#             resethandler()
#             tkwm.title(base,mzXMLPath)
#             tclvalue(rawpathvar) <- mzXMLPath

#             if (tmpArr[2] == "timecourse") {
#                 callFromTimeCourse <<- TRUE
#             } else {
#                 callFromTimeCourse <<- FALSE
#             }
#             tclvalue(mass) <- tmpArr[5]
#             if (substr(tmpArr[4], 1, 2) == "RT") {
#                 tclvalue(rt) <- as.numeric(substr(tmpArr[4], 3, nchar(tmpArr)))*60
#             } else if (substr(tmpArr[4], 1, 2) == "SC") {
#                 tclvalue(rt) <- acqNumToRt(xr, as.numeric(substr(tmpArr[4], 3, nchar(tmpArr))))*60
#             }
#             thisrec <<- tmpArr[7]
#             # if (grepl("XXX", tmpArr[7])) {

#             # }
#             zoom <- 5
#             assign("zoomlevel", zoom, envir=.GlobalEnv)

#             tclvalue(rtw) <- rtwl[zoom]
#             peakhandler()
#         }
#     } else {
#         noFile <<- TRUE
#     }
# }

loadInput <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    fminput <- file(args[1], open="r", encoding="UTF-16LE")
    lines <- readLines(fminput, warn=FALSE)[1]
    close(fminput)

    tmpArr <- strsplit(lines, "\\|")[[1]]
    rawFilePath <- tmpArr[6]
    fileNameNoExt <- sub("^([^.]*).*", "\\1", basename(rawFilePath))
    mzXMLPath <- paste(dirname(rawFilePath), paste(fileNameNoExt, ".mzXML", sep=""), sep="/")

    if (file.exists(mzXMLPath)) {
        xr <<- tryCatch(xcmsRaw(mzXMLPath, profstep=0, includeMSn=TRUE), error=function(ex) NULL)
        # assign("xr", xrtmp, envir=.GlobalEnv)
        if (!is.null(xr)) {
            # ms2spec <- tryCatch(msn2xcmsRaw(xcmsRaw(mzXMLPath, profstep=0, includeMSn=TRUE)), error=function(ex) NULL)
            # ms2 <<- findPeaks(xr, method="MS1")
            # ms2t <- findPeaks(ms2spec, method="MS1")
            # assign("ms2", ms2t, envir=.GlobalEnv)
            # assign("ms2", ms2spec, envir=.GlobalEnv)
            resethandler()
            tkwm.title(base,mzXMLPath)
            tclvalue(rawpathvar) <- mzXMLPath

            if (tmpArr[2] == "timecourse") {
                callFromTimeCourse <<- TRUE
            } else {
                callFromTimeCourse <<- FALSE
            }
            tclvalue(mass) <- tmpArr[5]
            if (substr(tmpArr[4], 1, 2) == "RT") {
                tclvalue(rt) <- as.numeric(substr(tmpArr[4], 3, nchar(tmpArr[4])))
            } else if (substr(tmpArr[4], 1, 2) == "SC") {
                tclvalue(rt) <- acqNumToRt(xr, as.numeric(substr(tmpArr[4], 3, nchar(tmpArr[4]))))
            }
            thisrec <<- tmpArr[7]
            # if (grepl("XXX", tmpArr[7])) {

            # }
            zoom <- 5
            assign("zoomlevel", zoom, envir=.GlobalEnv)

            tclvalue(rtw) <- rtwl[zoom]
            peakhandler()
        }
    } else {
        noFile <<- TRUE
    }
}

update <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    outputdir <- dirname(args[1])
    outputpath <- paste(outputdir, "peakcalcr.txt", sep="/")
    outputf <- file(outputpath)
    if (!is.na(outputdir)) {
        if (!is.na(as.numeric(tclvalue(peakint)))) {
            output <- paste(thisrec, callFromTimeCourse, abs(as.numeric(tclvalue(peakint))), as.numeric(tclvalue(peaksn)), (as.numeric(tclvalue(peakrtmax))-as.numeric(tclvalue(peakrtmin)))/60, as.numeric(tclvalue(peakrt))/60, sep=",")
        } else {
            output <- paste(thisrec, callFromTimeCourse, "ND", 1, "ND", "ND", sep=",")
        }
        write(output, file=outputf, append=FALSE)
        close(outputf)
    }
}

onLeftClick <- function(x,y) {
    tclvalue(peakrt) <- ""
    tclvalue(peakrtmin) <- ""
    tclvalue(peakrtmax) <- ""
    tclvalue(peakint) <- ""
    tclvalue(peaksn) <- ""

    width <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))

    xMin <- parPlotSize[1] * width
    xRange <- (parPlotSize[2] - parPlotSize[1]) * width
    rangeX <- (usrCoords[2] - usrCoords[1])/xRange

    yMin <- parPlotSize[3] * height
    yRange <- (parPlotSize[4] - parPlotSize[3]) * height
    rangeY <- (usrCoords[4] - usrCoords[3])/yRange

    xClick <- as.numeric(x)
    yClick <- height - as.numeric(y)
    xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX
    yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY
    manualLeft <<- xPlotCoord
    manualRight <<- NULL
    tkrreplot(img)
}

onLeftDrag <- function(x,y) {
    width <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))

    xMin <- parPlotSize[1] * width
    xRange <- (parPlotSize[2] - parPlotSize[1]) * width
    rangeX <- (usrCoords[2] - usrCoords[1])/xRange

    yMin <- parPlotSize[3] * height
    yRange <- (parPlotSize[4] - parPlotSize[3]) * height
    rangeY <- (usrCoords[4] - usrCoords[3])/yRange

    xClick <- as.numeric(x)
    yClick <- height - as.numeric(y)
    xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX
    yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY
    manualRight <<- xPlotCoord
    tkrreplot(img)
}

onLeftRelease <- function(x,y) {
    width <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))

    xMin <- parPlotSize[1] * width
    xRange <- (parPlotSize[2] - parPlotSize[1]) * width
    rangeX <- (usrCoords[2] - usrCoords[1])/xRange

    yMin <- parPlotSize[3] * height
    yRange <- (parPlotSize[4] - parPlotSize[3]) * height
    rangeY <- (usrCoords[4] - usrCoords[3])/yRange

    xClick <- as.numeric(x)
    yClick <- height - as.numeric(y)
    xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX
    yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY
    manualRight <<- xPlotCoord
    print(sort(c(manualLeft,manualRight)))
    tkrreplot(img)

    mzrange <- pmMz(as.numeric(tclvalue(mass)), as.numeric(tclvalue(ppm)))
    rtrange <- sort(c(manualLeft, manualRight))
    scrange <- sort(rtRangeToScanRange(xr, rtrange))
    if (!any(is.na(scrange))) {
        roi <- list(list(mzmin=mzrange[1], mzmax=mzrange[2], scmin=scrange[1], scmax=scrange[2]))
        peaks <- centWaveGB(xr, peakwidth=PEAKWIDTHG, ppm=as.numeric(tclvalue(ppm)), prefilter=PREFILTER, snthresh=SNTHRESH, 
            integrate=INTEGRATE, mzdiff=MZDIFF, ROI.list=roi, verbose.columns=TRUE, fitgauss=TRUE)

        print("ROI")
        print(peaks[])

        if (nrow(peaks) == 0) {
            target <- list(mzrange=mzrange, rtrange=rtrange, peakwidth=PEAKWIDTHG)
            targetedPeaks <- centWaveGB(xr, peakwidth=PEAKWIDTHG, ppm=massWindow, prefilter=PREFILTER, snthresh=SNTHRESH, 
                integrate=INTEGRATE, mzdiff=MZDIFF, targets=list(target), verbose.columns=TRUE, fitgauss=TRUE)

            print("targeted")
            print(targetedPeaks[])
            if (nrow(targetedPeaks) == 0) {
                baseline <- 1
                sdnoise <- 1
            } else {
                maxPeak <- targetedPeaks[which.max(targetedPeaks[,"intb"]),,drop=FALSE]
                baseline <- maxPeak[,"baseline"]
                sdnoise <- maxPeak[,"sdnoise"]
            }
        } else {
            maxPeak <- peaks[which.max(peaks[,"intb"]),,drop=FALSE]
            baseline <- maxPeak[,"baseline"]
            sdnoise <- maxPeak[,"sdnoise"]
        }
    }
    peakEIC <- rawEIC(xr, mzrange=mzrange, rtrange=rtrange)
    pseudoPeak <- pseudoPeakArea(xr, mzrange=mzrange, scanrange=scrange, maxIntensity=max(peakEIC$intensity), baseline=baseline, sdnoise=sdnoise)

    tclvalue(peakrt) <- xr@scantime[peakEIC$scan[which.max(peakEIC$intensity)]]
    tclvalue(peakrtmin) <- rtrange[1]
    tclvalue(peakrtmax) <- rtrange[2]
    tclvalue(peakint) <- pseudoPeak$intb
    tclvalue(peaksn) <- pseudoPeak$sn
    update()
}

scriptDir <- currentDir()
if (length(scriptDir) == 0) scriptDir <- "."
xcmsScriptPath <- paste(sep="/", scriptDir, xcmsScript)
source(xcmsScriptPath)

xr <<- NULL
# ms2 <<- NULL
thisrec <<- ""

parPlotSize <- c()
usrCoords <- c()
manualLeft <<- NULL
manualRight <<- NULL

# width <<- NULL
# height <<- NULL
# xMin <<- NULL
# xRange <<- NULL
# rangeX <<- NULL
# yMin <<- NULL
# yRange <<- NULL
# rangeY <<- NULL

mass <- tclVar()
ppm <- tclVar(10)
rt <- tclVar()
rtw <- tclVar(rtwl[1])
rawpathvar <- tclVar()
mpeakleft <- tclVar()
mpeakright <- tclVar()

peakrt <- tclVar()
peakrtmin <- tclVar()
peakrtmax <- tclVar()
peakmz <- tclVar()
peakint <- tclVar()
peaksn <- tclVar()

args <- commandArgs(trailingOnly = TRUE)

noFile <<- FALSE
callFromTimeCourse <<- FALSE

invisible({
    base <- tktoplevel()
	tkwm.resizable(base,FALSE,FALSE)
    menu <- tkmenu(base)
    tkconfigure(base,menu=menu)
    filemenu <- tkmenu(menu,tearoff=FALSE)
    tkadd(filemenu,"command",label="Quit",command=function() tkdestroy(base))
    tkadd(menu,"cascade",label="File",menu=filemenu)
    tkfocus(base)

    tkwm.title(base,'Manual')

    mainfrm <- tkframe(base)

    rawpathL <- tkframe(mainfrm)
    tkpack(tklabel(rawpathL,text='mzXML File Path',anchor="w"),side='left')
    if (length(args) == 0) {
        tkpack(tkbutton(rawpathL,text='Find',command=rawfilehandler),side='right')
    }
    tkpack(tkentry(rawpathL,width=75,textvariable=rawpathvar),side='right')
    tkpack(mainfrm,rawpathL)

    plotf <- tkframe(mainfrm)
    img <- tkrplot(plotf,doplot,hscale=1,vscale=1)

    tkpack(plotf,img)
    plotbutfrm <- tkframe(plotf,borderwidth=0,relief="groove")
    peakinbut <- tkbutton(plotbutfrm,text='+',width=2,command=zoomin)
    tkpack(plotbutfrm,peakinbut,side="left")
    peakoutbut <- tkbutton(plotbutfrm,text='-',width=2,command=zoomout)
    tkpack(plotbutfrm,peakoutbut,side="left")
    peakleftbut <- tkbutton(plotbutfrm,text='<',width=2,command=leftshift)
    tkpack(plotbutfrm,peakleftbut,side="left")
    peakrightbut <- tkbutton(plotbutfrm,text='>',width=2,command=rightshift)
    tkpack(plotbutfrm,peakrightbut,side="left")
    tkpack(plotf,plotbutfrm,side="bottom")

    inputfieldfrm <- tkframe(mainfrm)

    line1L <- tkframe(inputfieldfrm)
    tkpack(tklabel(line1L,text='',width=20),side='left')
    tkpack(inputfieldfrm,line1L)
    inputL <- tkframe(inputfieldfrm)
    tkpack(tklabel(inputL,text='Input',width=20,anchor="w"),side='left')
    tkpack(inputfieldfrm,inputL)
    massL <- tkframe(inputfieldfrm)
    tkpack(tklabel(massL,text='Mass',width=28,anchor="w"),side='left')
    tkpack(tkentry(massL,width=10,textvariable=mass),side='right')
    tkpack(inputfieldfrm,massL)
    ppmL <- tkframe(inputfieldfrm)
    tkpack(tklabel(ppmL,text='PPM',width=28,anchor="w"),side='left')
    tkpack(tkentry(ppmL,width=10,textvariable=ppm),side='right')
    tkpack(inputfieldfrm,ppmL)

    line2L <- tkframe(inputfieldfrm)
    tkpack(tklabel(line2L,text='',width=20),side='left')
    tkpack(inputfieldfrm,line2L)
    peakL <- tkframe(inputfieldfrm)
    tkpack(tklabel(peakL,text='Peak',width=20,anchor="w"),side='left')
    tkpack(inputfieldfrm,peakL)
    peakrtL <- tkframe(inputfieldfrm)
    tkpack(tklabel(peakrtL,text='Retention Time',width=20,anchor="w"),side='left')
    tkpack(tkentry(peakrtL,width=20,textvariable=peakrt,state="readonly"),side='right')
    tkpack(inputfieldfrm,peakrtL)
    peakrtminL <- tkframe(inputfieldfrm)
    tkpack(tklabel(peakrtminL,text='Retention Time Min',width=20,anchor="w"),side='left')
    tkpack(tkentry(peakrtminL,width=20,textvariable=peakrtmin,state="readonly"),side='right')
    tkpack(inputfieldfrm,peakrtminL)
    peakrtmaxL <- tkframe(inputfieldfrm)
    tkpack(tklabel(peakrtmaxL,text='Retention Time Max',width=20,anchor="w"),side='left')
    tkpack(tkentry(peakrtmaxL,width=20,textvariable=peakrtmax,state="readonly"),side='right')
    tkpack(inputfieldfrm,peakrtmaxL)
    peakintL <- tkframe(inputfieldfrm)
    tkpack(tklabel(peakintL,text='Intensity',width=20,anchor="w"),side='left')
    tkpack(tkentry(peakintL,width=20,textvariable=peakint,state="readonly"),side='right')
    tkpack(inputfieldfrm,peakintL)
    peaksnL <- tkframe(inputfieldfrm)
    tkpack(tklabel(peaksnL,text='Signal to Noise',width=20,anchor="w"),side='left')
    tkpack(tkentry(peaksnL,width=20,textvariable=peaksn,state="readonly"),side='right')
    tkpack(inputfieldfrm,peaksnL)
    line3L <- tkframe(inputfieldfrm)
    tkpack(tklabel(line3L,text='',width=20),side='left')
    tkpack(inputfieldfrm,line3L)

    butfrm <- tkframe(inputfieldfrm,borderwidth=0,relief="groove")
    plotbut <- tkbutton(butfrm,text='Replot',width=7,command=function()tkrreplot(img))
    peakbut <- tkbutton(butfrm,text='Peak',width=7,command=peakhandler)
    resetbut <- tkbutton(butfrm,text='Reset',width=7,command=resethandler)
    quitbut <- tkbutton(butfrm,text='Quit',width=7,command=function()tkdestroy(base))

    tkpack(butfrm,plotbut,side="left")
    tkpack(butfrm,peakbut,side="left")
    tkpack(butfrm,resetbut,side="left")
    tkpack(butfrm,quitbut,side="left")

    if (length(args) == 1) {
        butfrm2 <- tkframe(inputfieldfrm,borderwidth=0,relief="groove")
        updatebut <- tkbutton(butfrm2,text='Update',width=7,command=update)
        tkpack(butfrm2,updatebut,side="left")
        reloadbut <- tkbutton(butfrm2,text='Reload',width=7,command=loadInput)
        tkpack(butfrm2,reloadbut,side="left")
        tkpack(inputfieldfrm,butfrm2,side="bottom")
    }

    tkpack(inputfieldfrm,butfrm,side="bottom")

    tkpack(mainfrm,inputfieldfrm,side="left")
    tkpack(mainfrm,plotf,side="right")

    if (length(args) == 1) {
        loadInput()
    } else if (length(args) > 1) {
        tk_messageBox(type="ok", message="Too many arguments.")
    }

    if (noFile) {
        tk_messageBox(type="ok", message="mzXML file does not exist.")
    }

    tkbind(img, "<Button-1>",onLeftClick)
    tkbind(img, "<B1-Motion>",onLeftDrag)
    tkbind(img, "<ButtonRelease-1>",onLeftRelease)
    tkconfigure(img,cursor="hand2")

    tkwait.window(base)
})