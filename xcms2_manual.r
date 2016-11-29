#!/usr/bin/env Rscript
# @aly

# Packages required for all functions in this script to run
# Remember to install and load them in scripts that source this script
# suppressPackageStartupMessages({
#     library(e1071)
#     library(xcms)
#     library(nortest)
#     })

########################################################################################
# Parameters for centwave/getPeaks.
# Put peakwidths (PEAKWIDTHL) in ascending order. For quantitation, the order where the 
#  results are returned is in decending order (reverse).
# PEAKWIDTHG is used for getPeaks/pseudoPeakCalc baseline adjusted peak area.
########################################################################################
RTWINDOW <<- 30

PREFILTER <<- c(1,100)
SNTHRESH <<- -Inf
INTEGRATE <<- 2
MZDIFF <<- -10000.0

PEAKWIDTHG <<- c(10,200)


########################################################################################
# Input: Mass, ppm
# Output: Mass range
# 
# This method calculates a mass range given a mass and a mass window (ppm).
# Numbers used are based on visual basic AutoFill
########################################################################################
pmMz <- function (x, ppm=10) {
    c(x - x * ppm * 0.000001, x + x * ppm * 0.000001)
}

########################################################################################
# Input: Retention time, retention time window (one side)
# Output: Retention time range
# 
# This method calculates a retention time window given a retention time and a retention
#  time window.
########################################################################################
pmRt <- function (x, rtwidth=10) {
    c(x - rtwidth, x + rtwidth)
}

########################################################################################
# Input: xcmsRaw object, retention time range
# Output: Scan number range
# 
# This method takes a retention time range, and converts it to a scan number range, by 
#  finding the closest retention time of the scan that is less than the minimum, and the
#  the closest retention time of the scan that is greater than the maximum, of the 
#  provided retention time range.
# This method is based on xcms source code. 
########################################################################################
rtRangeToScanRange <- function (xcmsRaw, rtrange) {
    rtrange <- range(rtrange)
    scanidx <- (xcmsRaw@scantime >= rtrange[1]) & (xcmsRaw@scantime <= rtrange[2])
    c(match(TRUE, scanidx), length(scanidx) - match(TRUE, rev(scanidx)))
}

########################################################################################
# Input: xcmsRaw object, acquisition number
# Output: Retention time closest to acquisition number
# 
# This method finds the retention time closest to the acquisition number provided, by 
#  first finding the closest scan number (scan index) and the finds the associated 
#  retention time.
########################################################################################
acqNumToRt <- function (xcmsRaw, aqNum) {
    xcmsRaw@scantime[as.integer(acqNumToScNum(xcmsRaw, aqNum))]
}

########################################################################################
# Input: xcmsRaw object, acquisition number
# Output: Scan number closest to acquisition number
# 
# This method finds the scan number (scan index) closest to the provided acquisition 
#  number by going through the list of acquisition numbers with associated scan numbers
#  and returning that scan number.
########################################################################################
acqNumToScNum <- function (xcmsRaw, aqNum) {
    which.min(abs(xcmsRaw@acquisitionNum - as.integer(aqNum)))
}

########################################################################################
# Input: mzXML path
# Output: xcmsRaw object
# 
# This method converts the mzXML at the path provided to a xcmsRaw object with 
#  "profstep=0".
########################################################################################
convert <- function (mzXMLPath) {
    return(xcmsRaw(mzXMLPath, profstep=0))
}

########################################################################################
# Input: xcmsRaw object, mass range (vector), RT range (vector)
# Output: EIC (raw)
# 
# This method geneerates an EIC with specified mass and retention time windows, and uses
#  the rawEIC method (not getEIC). Adjustments are made to the output from rawEIC to 
#  include retention times instead of scan numbers, by using the indices of scan numbers
#  as the indices of retention times (scantime).
# Final output that is returned has retention times and intensities associated with
#  them.
########################################################################################
getChromatogram <- function(xcmsRaw, mzrange, rtrange) {
    EIC <- rawEIC(xcmsRaw, mzrange, rtrange)
    EIC$rt <- unlist(lapply(EIC$scan, function(scan) xcmsRaw@scantime[scan]))
    EIC$scan <- NULL
    return(EIC)
}

########################################################################################
# Input: xcmsRaw object, mass list (vector), retention time, mass window (ppm), 
#  svm (for quality peak metrics), metric (boolean for calculating quality peak 
#  metrics), peakCalcLoader (boolean for only returning detected peaks, only used for
#  peakCalcLoader), rtw (default to constant at top of file, used for manual peak
#  detector)
# Output: List of peaks detected, standard deviation of peak apex retention times, and 
#  standard deviation of peak shapes
# 
# This follows what the "Run" method does in PEAKCALC.bas in visual basic.
# This method first runs selectMass on the masses provided, along with other parameters.
#  Then it goes through the peaks returned, either detected with centWave, or calculated
#  with getPeaks. Then it tries to use the retention time of a detected peak or the 
#  peak with the highest intensity, and redo the peak detection/calculation for the 
#  other peaks supplied. Once that is done, then it will calculate the quality peak 
#  metrics for those all peaks if required. Then the standard deviations are calcualted
#  with the retention times of the peaks (apex) and peak shape. This part is based on
#  the visual basic PEAKCALC.bas.
########################################################################################
peakCalc <- function(xcmsRaw, masslist, rt, massWindow, svmObject=NULL, metric=FALSE, peakCalcLoader=FALSE, rtw=RTWINDOW) {
    p <- lapply(masslist, function(mass) selectMass(xcmsRaw, mass, rt, massWindow, peakCalcLoader, rtw))
    nPeaks <- length(p)
    peakAreas <- unlist(lapply(p, '[[', "peakArea"))

    if (!all(peakAreas > 0) || any(is.na(peakAreas))) {
        if (any(na.omit(peakAreas) > 0)) {
            maxPos <- which.max(peakAreas)
        } else {
            maxPos <- which.min(peakAreas)
            maxPos <- ifelse(length(maxPos) > 0, maxPos, -1)
        }

        if (maxPos > 0) {
            for (j in (1:nPeaks)[-maxPos]) {
                if (p[[j]][["peakArea"]] < 0 || is.na(p[[j]][["peakArea"]])) {
                    p[[j]] <- selectMass(xcmsRaw, masslist[j], rt, massWindow, peakCalcLoader, rtw)
                }
            }
        }
    }

    if (metric) {
        for (l in 1:nPeaks) {
            if (!is.na(p[[l]][["peakArea"]])) {
                m <- metrics(xcmsRaw, svmObject, rt, p[[l]][["peakApexRT"]], p[[l]][["peakLeftRT"]], p[[l]][["peakRightRT"]], masslist[l], massWindow)
                if (length(m) < 1) {
                    min <- head(xcmsRaw@scantime, n=1)
                    max <- tail(xcmsRaw@scantime, n=1)
                    peakLeft <- p[[l]][["peakApexRT"]] - 30
                    if (peakLeft < min) {
                        peakLeft <- min
                    }
                    peakRight <- p[[l]][["peakApexRT"]] + 30
                    if (peakRight > max) {
                        peakRight <- max
                    }
                    tryCatch(p[[l]][["metric"]] <- metrics(xcmsRaw, svmObject, rt, p[[l]][["peakApexRT"]], peakLeft, peakRight, masslist[l], massWindow), 
                        error=function(ex) p[[l]][["metric"]] <- "no score")
                } else {
                    p[[l]][["metric"]] <- m
                }
            } else {
                p[[l]][["metric"]] <- "no score"
            }
        }
    }

    stDevPeakApex <- standardDeviation(unlist(lapply(p, '[[', "peakApexRT")))

    lr <- rep(NA, nPeaks)
    for (k in 1:nPeaks) {
        if (is.na(p[[k]][["peakApexRT"]]) || p[[k]][["peakLeftRT"]] == p[[k]][["peakRightRT"]]) {
            lr[k] <- 0
        } else {
            lr[k] = (2 * p[[k]][["peakApexRT"]] - p[[k]][["peakLeftRT"]] - p[[k]][["peakRightRT"]]) / 
            (p[[k]][["peakRightRT"]] - p[[k]][["peakLeftRT"]])
        }
    }

    stDevPeakShape <- standardDeviation(lr)

    return(list(peaks=p, stDevPeakApex=stDevPeakApex, stDevPeakShape=stDevPeakShape))
}

########################################################################################
# Input: xcmsRaw object, mass, retention time, mass window (ppm), peakCalcLoader 
#  (boolean for only returning detected peaks, only used for peakCalcLoader), rtw 
#  (default to constant in top of file, used for manual peak detector)
# Output: Peak and associated information.
# 
# This follows what the "SelectMass" method does in PEAKCALC.bas in visual basic, but 
#  with modifications to run xcms and its peak detection, centWave.
# First this creates a mass range and a scan number range. These are used for the ROI
#  that will be used for centWave. Then centWave uses that roi, along with the 
#  parameters at the top of the script (and loops through multiple peak widths). The 
#  first run of centWave is with the whole scan range of the xcmsRaw object, while 
#  the second run of centWave is with the scan range calculated, and the third run is 
#  on double of that scan range. Then three sets of those can run, depending if the 
#  previous return a result or not. The first is the default centWave, the second is
#  with doubling the mass window, and checking if the peak detected fall within the
#  original mass window, and the third is with manually generating a ROI but with 
#  specific parameters (the main difference should be MINCENTROIDS, which in centWave
#  is defaulted to 4 as minimum). Then the outputs are merged, with the peaks detected 
#  with a larger peak width parameter at the top of the list. With that list, then it 
#  is filtered by retention times (peaks at the top of the list, or those that are 
#  detected with a larger peak width parameter), and NA's are removed. 
#  Now it goes through the case where centWave either detects something or not. If it 
#  does, then it will go through the list of peaks and find the one with the largest
#  intensity. If not, then it will run getPeaks with the center being the highest 
#  intensity point, and try to calculate a baseline adjusted intensity and signal to 
#  noise with a part of centWave method extracted from xcms source. 
# 
# http://metabolomics-forum.com/viewtopic.php?f=25&t=154&sid=1af1b9cff6303e147ab8e0c37fcd8626
########################################################################################
selectMass <- function(xcmsRawO, mass, rt, massWindow, peakCalcLoader=FALSE, rtw) {
    emptyMat <- matrix(c(rep(NA, 10)), nrow=1)
    noPeak <- list(peakArea=NA, peakApexRT=NA, peakLeftRT=NA, peakRightRT=NA, peakWidth=NA, sn=NA, metric=NA)
    mzrange <- pmMz(mass, massWindow)
    rtrange <- pmRt(rt, rtw)
    scrange <- rtRangeToScanRange(xcmsRawO, rtrange)
    if (any(is.na(scrange))) return(noPeak)
    roi <- list(list(mzmin=mzrange[1], mzmax=mzrange[2], scmin=scrange[1], scmax=scrange[2]))
    peaks <- centWaveGB(xcmsRawO, peakwidth=PEAKWIDTHG, ppm=massWindow, prefilter=PREFILTER, snthresh=SNTHRESH, 
        integrate=INTEGRATE, mzdiff=MZDIFF, ROI.list=roi, verbose.columns=TRUE, fitgauss=TRUE)

    if (nrow(peaks) == 0) { # if "No Peaks" or "omz == 0" (mass outside of raw file mass range) error
        if (peakCalcLoader) { # do not want pseudo peaks in detecting peaks with peakCalcLoader
            return(noPeak)
        }
        eic <- tryCatch(rawEIC(xcmsRawO, mzrange=mzrange, scanrange=scrange), error = function(ex) NA)
        if (all(eic$intensity == 0) || is.na(eic)) {
            return(noPeak)
        }

        target <- list(mzrange=mzrange, rtrange=rtrange, peakwidth=PEAKWIDTHG)
        targetedPeaks <- centWaveGB(xcmsRawO, peakwidth=PEAKWIDTHG, ppm=massWindow, prefilter=PREFILTER, snthresh=SNTHRESH, 
            integrate=INTEGRATE, mzdiff=MZDIFF, targets=list(target), verbose.columns=TRUE, fitgauss=TRUE)

        if (nrow(targetedPeaks) == 0) {
            return(noPeak)
        } else {
            peak <- targetedPeaks[which.max(targetedPeaks[,"intb"]),,drop=FALSE]
            print(peak[])
            tmpeic <- rawEIC(xcmsRawO, mzrange=c(peak[,"mzmin"], peak[,"mzmax"]), rtrange=c(peak[,"rtmin"], peak[,"rtmax"]))
            maxort <- xcmsRawO@scantime[tmpeic$scan[which.max(tmpeic$intensity)]]
            # return(list(peakArea=-peak[,"intb"], peakApexRT=peak[,"rt"], peakLeftRT=peak[,"rtmin"], peakRightRT=peak[,"rtmax"], peakWidth=peak[,"rtmax"]-peak[,"rtmin"], sn=peak[,"sn"], metric=NA))
            return(list(peakArea=-peak[,"intb"], peakApexRT=maxort, peakLeftRT=peak[,"rtmin"], peakRightRT=peak[,"rtmax"], peakWidth=peak[,"rtmax"]-peak[,"rtmin"], sn=peak[,"sn"], metric=NA))
        }
    } else {
        # pl <- unlist(filteredPeaks[which.min(abs(filteredPeaks[,"V4"] - as.numeric(rt))),]) # Pick peak with rt closest to input rt
        # pl <- unlist(filteredPeaks[which.max(peaks[,"intb"]),]) # Pick peak with highest intensity
        peak <- peaks[which.max(peaks[,"intb"]),,drop=FALSE]
        print(peak[])
        tmpeic <- rawEIC(xcmsRawO, mzrange=c(peak[,"mzmin"], peak[,"mzmax"]), rtrange=c(peak[,"rtmin"], peak[,"rtmax"]))
        maxort <- xcmsRawO@scantime[tmpeic$scan[which.max(tmpeic$intensity)]]
        # return(list(peakArea=peak[,"intb"], peakApexRT=peak[,"rt"], peakLeftRT=peak[,"rtmin"], peakRightRT=peak[,"rtmax"], peakWidth=peak[,"rtmax"]-peak[,"rtmin"], sn=peak[,"sn"], metric=NA))
        return(list(peakArea=peak[,"intb"], peakApexRT=maxort, peakLeftRT=peak[,"rtmin"], peakRightRT=peak[,"rtmax"], peakWidth=peak[,"rtmax"]-peak[,"rtmin"], sn=peak[,"sn"], metric=NA))
    }
}


pseudoPeakArea <- function(object, mzrange, scanrange, maxIntensity, baseline, sdnoise) {
    maxint <- unname(maxIntensity)

    sn <- round((maxint - baseline) / sdnoise)

    tmpeic <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)
    tmpint <- tmpeic$intensity

    irt <- unlist(object@scantime[tmpeic$scan])
    int_b <- tmpint - baseline
    int_b[which(int_b < 0)] <- 0
    area_b <- sum(diff(irt)*sapply(2:length(int_b), function(x){mean(int_b[c(x-1, x)])}))

    return(list(intb=unname(area_b), sn=unname(sn)))
}


# findPeaks <- function(x, thresh=0) {
#     pks <- which(diff(sign(diff(x, na.pad=FALSE)),na.pad=FALSE) < 0) + 2
#     if(!missing(thresh)) {
#      pks[x[pks-1]-x[pks] > thresh]
#     } else pks
# }

# findValleys <- function(x, thresh=0) {
#     pks <- which(diff(sign(diff(x, na.pad=FALSE)),na.pad=FALSE) > 0) + 2
#     if(!missing(thresh)) {
#      pks[x[pks-1]-x[pks] > thresh]
#     } else pks
# }





# Ginkgo Bioworks' XCMS, based on XCMS 1.39.2
# https://github.com/benjiec/xcms
# from xcmsRaw.R
centWaveOnROI <- function(object, scanrange, basenames, verbosenames,
                                               roi, roi_index, peakwidth=c(20,50), snthresh=10,
                                               mzCenterFun="wMean",
                                               integrate=1, fitgauss=FALSE,
                                               verbose.columns=FALSE, targeted=FALSE, sleep=0) {

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(object@scantime))) / 2)

    if (length(z <- which(scalerange==0)))
        scalerange <- scalerange[-z]

    if (length(scalerange) < 1)
        stop("No scales ? Please check peak width!\n")

    if (length(scalerange) > 1)
        scales <- seq(from=scalerange[1], to=scalerange[2], by=2)  else
    scales <- scalerange;

    minPeakWidth <- scales[1];
    noiserange <- c(minPeakWidth*3, max(scales)*3);
    maxGaussOverlap <- 0.5;
    minPtsAboveBaseLine <- max(4,minPeakWidth-2);
    minCentroids <- minPtsAboveBaseLine ;
    scRangeTol <- maxDescOutlier <- floor(minPeakWidth/2);

    scantime <- object@scantime
    Nscantime <- length(scantime)
    peaks <- peakinfo <- NULL

    if (!is.null(roi)) {  # here to avoid re-indenting everything
        N <- roi$scmax - roi$scmin + 1

        mzrange <- c(roi$mzmin,roi$mzmax)
        scrange <- c(roi$scmin,roi$scmax)
        ## scrange + noiserange, used for baseline detection and wavelet analysis
        sr <- c(max(scanrange[1],scrange[1] - max(noiserange)),min(scanrange[2],scrange[2] + max(noiserange)))
        eic <- xcms:::rawEIC(object,mzrange=mzrange,scanrange=sr)
        d <- eic$intensity
        td <- sr[1]:sr[2]
        scan.range <- c(sr[1],sr[2])
        ## original mzROI range
        mzROI.EIC <- xcms:::rawEIC(object,mzrange=mzrange,scanrange=scrange)
        omz <- xcms:::rawMZ(object,mzrange=mzrange,scanrange=scrange)
        if (all(omz == 0))
            return (peaks)
        od  <- mzROI.EIC$intensity
        otd <- mzROI.EIC$scan
        if (all(od == 0))
            return (peaks)

        ##  scrange + scRangeTol, used for gauss fitting and continuous data above 1st baseline detection
        ftd <- max(td[1], scrange[1] - scRangeTol) : min(td[length(td)], scrange[2] + scRangeTol)
        fd <- d[match(ftd,td)]

        ## 1st type of baseline: statistic approach
        if (N >= 10*minPeakWidth)  ## in case of very long mass trace use full scan range for baseline detection
            noised <- xcms:::rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity
        else
            noised <- d
        ## 90% trimmed mean as first baseline guess
        noise <- xcms:::estimateChromNoise(noised, trim=0.05, minPts=3*minPeakWidth)

        # ## any continuous data above 1st baseline ?
        # if (!targeted && !xcms:::continuousPtsAboveThreshold(fd,threshold=noise,num=minPtsAboveBaseLine))
        #     return(peaks)

        ## 2nd baseline estimate using not-peak-range
        lnoise <- xcms:::getLocalNoiseEstimate(d,td,ftd,noiserange,Nscantime, threshold=noise,num=minPtsAboveBaseLine)

        ## Final baseline & Noise estimate
        baseline <- max(1,min(lnoise[1],noise))
        sdnoise <- max(1,lnoise[2])
        sdthr <-  sdnoise * snthresh

        if (!targeted && !xcms:::continuousPtsAboveThreshold(fd,threshold=noise,num=minPtsAboveBaseLine))
            # return(peaks)
            return(list(peaks=peaks,baseline=baseline))

        ## is there any data above S/N * threshold ?
        if (!(any(fd - baseline >= sdthr)))
            return(peaks)

        wCoefs <- xcms:::MSW.cwt(d, scales=scales, wavelet='mexh')
        if (!(!is.null(dim(wCoefs)) && any(wCoefs- baseline >= sdthr)))
            return(peaks)

        if (td[length(td)] == Nscantime) ## workaround, localMax fails otherwise
            wCoefs[nrow(wCoefs),] <- wCoefs[nrow(wCoefs)-1,] * 0.99
        localMax <- xcms:::MSW.getLocalMaximumCWT(wCoefs)
        rL <- xcms:::MSW.getRidge(localMax)
        wpeaks <- sapply(rL,
                         function(x) {
                             w <- min(1:length(x),ncol(wCoefs))
                             any(wCoefs[x,w]- baseline >= sdthr)
                         })
        if (any(wpeaks)) {
            wpeaksidx <- which(wpeaks)
            new_peaks <- list()
            ## check each peak in ridgeList
            for (p in 1:length(wpeaksidx)) {
                opp <- rL[[wpeaksidx[p]]]
                pp <- unique(opp)
                if (length(pp) >= 1) {
                    dv <- td[pp] %in% ftd
                    if (any(dv)) { ## peaks in orig. data range
                        ## Final S/N check
                        if (any(d[pp[dv]]- baseline >= sdthr)) {
                            ## try to decide which scale describes the peak best
                            coef <- numeric(length(opp))
                            irange = rep(ceiling(scales[1]/2),length(opp))
                            for (k in 1:length(opp)) {
                                kpos <- opp[k]
                                r1 <- ifelse(kpos-irange[k] > 1,kpos-irange[k],1)
                                r2 <- ifelse(kpos+irange[k] < length(d),kpos+irange[k],length(d))
                                if (dim(wCoefs)[2] >= k) { coef[k] <- wCoefs[opp[k],k] }
                                else {
                                    cat("Missing CWT coefficient for", scantime[td[opp[k]]], "scale", k, "\n");
                                    coef[k] <- 0
                                }
                            }
                            maxpc <- which.max(coef)
                            if (length(maxpc) > 1) {
                                # use the narrowest shape that fits the wavelet
                                best.scale.nr <- maxpc[1]
                            } else  best.scale.nr <- maxpc

                            best.scale <-  scales[best.scale.nr]
                            best.scale.pos <- opp[best.scale.nr]
                            best.center <- best.scale.pos

                            pprange <- min(pp):max(pp)
                            ## maxint <- max(d[pprange])
                            lwpos <- max(1,best.scale.pos - best.scale)
                            rwpos <- min(best.scale.pos + best.scale,length(td))
                            p1 <- match(td[lwpos],otd)[1]
                            p2 <- match(td[rwpos],otd); p2 <- p2[length(p2)]
                            if (is.na(p1)) p1<-1
                            if (is.na(p2)) p2<-N
                            mz.value <- omz[p1:p2]
                            if (length(mz.value[which(mz.value > 0)]) == 0)
                              next
                            mz.int <- od[p1:p2]
                            maxint <- max(mz.int)

                            ## re-calculate m/z value for peak range
                            mzrange <- range(mz.value[which(mz.value > 0)])
                            mzmean <- do.call(mzCenterFun,list(mz=mz.value, intensity=mz.int))

                            ## Compute dppm only if needed
                            dppm <- NA
                            if (verbose.columns)
                                if (length(mz.value) >= (minCentroids+1))
                                    dppm <- round(min(xcms:::running(abs(diff(mz.value)) /(mzrange[2] *  1e-6),fun=max,width=minCentroids))) else
                            dppm <- round((mzrange[2]-mzrange[1]) /  (mzrange[2] *  1e-6))

                            p = c(mzmean,mzrange,           ## mz
                                             NA,NA,NA,                   ## rt, rtmin, rtmax,
                                             NA,                         ## intensity (sum)
                                             NA,                         ## intensity (-bl)
                                             maxint,                     ## max intensity
                                             round((maxint - baseline) / sdnoise),  ##  S/N Ratio
                                             NA,                         ## Gaussian RMSE
                                             NA,NA,NA,                   ## Gaussian Parameters
                                             roi_index,                  ## ROI Position
                                             dppm,                       ## max. difference between the [minCentroids] peaks in ppm
                                             best.scale,                 ## Scale
                                             td[best.scale.pos], td[lwpos], td[rwpos],  ## Peak positions guessed from the wavelet's (scan nr)
                                             # NA,NA )                    ## Peak limits (scan nr)
                                             NA,NA,NA,baseline,sdnoise)

                            # check and ignore duplicates
                            duplicated <- FALSE
                            if (length(new_peaks) > 0) {
                              for (npi in 1:length(new_peaks)) {
                                b = new_peaks[[npi]]
                                if (all(b[which(!is.na(b))] == p[which(!is.na(p))])) { duplicated <- TRUE; break; }
                              }
                            }
                            if (!duplicated) {
                              new_peaks[[length(new_peaks)+1]] <- p
                              peaks <- rbind(peaks, p, deparse.level=0)
                              ## Peak positions guessed from the wavelet's
                              pinfo <- c(best.scale, best.scale.nr, best.scale.pos, lwpos, rwpos, best.center)
                              peakinfo <- rbind(peakinfo, pinfo, deparse.level=0)
                            }
                        }
                    }
                }
            }  ##for
        } ## if


        ##  postprocessing
        if (!is.null(peaks)) {
            colnames(peaks) <- c(basenames, verbosenames)

            colnames(peakinfo) <- c("scale","scaleNr","scpos","scmin","scmax","scmid")

            for (p in 1:dim(peaks)[1]) {
                ## find minima, assign rt and intensity values
                if (integrate == 1) {
                    lm <- xcms:::descendMin(wCoefs[,peakinfo[p,"scaleNr"]], istart= peakinfo[p,"scpos"])
                    gap <- all(d[lm[1]:lm[2]] == 0) ## looks like we got stuck in a gap right in the middle of the peak
                    if ((lm[1]==lm[2]) || gap )## fall-back
                        lm <- xcms:::descendMinTol(d, startpos=c(peakinfo[p,"scmin"], peakinfo[p,"scmax"]), maxDescOutlier)
                } else if (integrate == 2)
                    lm <- xcms:::descendMinTol(d,startpos=c(peakinfo[p,"scmin"],peakinfo[p,"scmax"]),maxDescOutlier)
                else {
                    lm <- c(peakinfo[p,"scmin"], peakinfo[p,"scmax"])
                    # narrow rt range down by skipping regions less than 10% of
                    # max at each end, considering those as trailing noise
                    maxo <- max(d[lm[1]:lm[2]])
                    pd <- d[lm[1]:lm[2]];
                    mino <- 0.1*maxo
                    lm.l <- xcms:::findEqualGreaterUnsorted(pd,mino)
                    lm.l <- max(1, lm.l-1)
                    lm.r <- xcms:::findEqualGreaterUnsorted(rev(pd),mino)
                    lm.r <- max(1, lm.r-1)
                    lm <- lm + c(lm.l-1, -(lm.r-1))
                }

                ## narrow down peak rt boundaries by skipping zeros
                pd <- d[lm[1]:lm[2]];
                ## instead of skipping zeros, skip regions with weak intensities
                maxo <- max(d[lm[1]:lm[2]])
                mino <- 0.01*maxo
                lm.l <- xcms:::findEqualGreaterUnsorted(pd,mino)
                lm.l <- max(1, lm.l - 1)
                lm.r <- xcms:::findEqualGreaterUnsorted(rev(pd),mino)
                lm.r <- max(1, lm.r - 1)
                lm <- lm + c(lm.l - 1, -(lm.r - 1) )

                peakrange <- td[lm]
                peaks[p,"rtmin"] <- scantime[peakrange[1]]
                peaks[p,"rtmax"] <- scantime[peakrange[2]]

                peaks[p,"maxo"] <- max(d[lm[1]:lm[2]])

                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]])/(peakrange[2] - peakrange[1])
                if (is.na(pwid))
                    pwid <- 1

                # see https://groups.google.com/forum/#!topic/xcms-devel/bdkvPGSWOU0
                # changing integration to account for non-uniform scan densities
                irt <- scantime[peakrange[1]:peakrange[2]]
                int <- d[lm[1]:lm[2]]
                area <- sum(diff(irt)*sapply(2:length(int), function(x){mean(int[c(x-1, x)])}))
                peaks[p, "into"] <- area

                int_b <- d[lm[1]:lm[2]] - baseline
                int_b[which(int_b < 0)] <- 0
                area_b <- sum(diff(irt)*sapply(2:length(int_b), function(x){mean(int_b[c(x-1, x)])}))
                peaks[p, "intb"] <- area_b

                peaks[p,"lmin"] <- lm[1]
                peaks[p,"lmax"] <- lm[2]
                peaks[p,"pwid"] <- pwid

                if (fitgauss) {
                    ## perform gaussian fits, use wavelets for inital parameters
                    md <- max(d[lm[1]:lm[2]]);d1 <- d[lm[1]:lm[2]]/md; ## normalize data for gaussian error calc.
                    pgauss <- xcms:::fitGauss(td[lm[1]:lm[2]],d[lm[1]:lm[2]],pgauss =
                                       list(mu=peaks[p,"scpos"],sigma=peaks[p,"scmax"]-peaks[p,"scmin"],h=peaks[p,"maxo"]))
                    rtime <- peaks[p,"scpos"]
                    if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                        gtime <- td[match(round(pgauss$mu),td)]
                        if (!is.na(gtime)) {
                            rtime <- gtime
                            peaks[p,"mu"] <- pgauss$mu; peaks[p,"sigma"] <- pgauss$sigma; peaks[p,"h"] <- pgauss$h;
                            peaks[p,"egauss"] <- sqrt((1/length(td[lm[1]:lm[2]])) * sum(((d1-xcms:::gauss(td[lm[1]:lm[2]],pgauss$h/md,pgauss$mu,pgauss$sigma))^2)))
                        }
                    }
                    peaks[p,"rt"] <- scantime[rtime]
                    ## avoid fitting side effects
                    if (peaks[p,"rt"] < peaks[p,"rtmin"])
                        peaks[p,"rt"] <- scantime[peaks[p,"scpos"]]
                } else peaks[p,"rt"] <- scantime[td[peakinfo[p,"scmid"]]]
            }
            peaks <- joinOverlappingPeaksGB(td,d,otd,omz,od,scantime,scan.range,peaks,maxGaussOverlap,mzCenterFun=mzCenterFun)
        }

        if ((sleep >0) && (!is.null(peaks))) {
            tdp <- scantime[td]; trange <- range(tdp)
            egauss <- paste(round(peaks[,"egauss"],3),collapse=", ")
            cdppm <- paste(peaks[,"dppm"],collapse=", ")
            csn <- paste(peaks[,"sn"],collapse=", ")
            par(bg = "white")
            l <- layout(matrix(c(1,2,3),nr=3,nc=1,byrow=T),heights=c(.5,.75,2));
            par(mar= c(2, 4, 4, 2) + 0.1)
            xcms:::plotRaw(object,mzrange=mzrange,rtrange=trange,log=TRUE,title='')
            title(main=paste(f,': ', round(mzrange[1],4),' - ',round(mzrange[2],4),' m/z , dppm=',cdppm,', EGauss=',egauss ,',  S/N =',csn,sep=''))
            par(mar= c(1, 4, 1, 2) + 0.1)
            image(y=scales[1:(dim(wCoefs)[2])],z=wCoefs,col=terrain.colors(256),xaxt='n',ylab='CWT coeff.')
            par(mar= c(4, 4, 1, 2) + 0.1)
            plot(tdp,d,ylab='Intensity',xlab='Scan Time');lines(tdp,d,lty=2)
            lines(scantime[otd],od,lty=2,col='blue') ## original mzbox range
            abline(h=baseline,col='green')
            bwh <- length(sr[1]:sr[2]) - length(baseline)
            if (odd(bwh)) {bwh1 <-  floor(bwh/2); bwh2 <- bwh1+1} else {bwh1<-bwh2<-bwh/2}
            if  (any(!is.na(peaks[,"scpos"])))
            {   ## plot centers and width found through wavelet analysis
                abline(v=scantime[na.omit(peaks[(peaks[,"scpos"] >0),"scpos"])],col='red')
            }
            abline(v=na.omit(c(peaks[,"rtmin"],peaks[,"rtmax"])),col='green',lwd=1)
            if (fitgauss) {
                tdx <- seq(min(td),max(td),length.out=200)
                tdxp <- seq(trange[1],trange[2],length.out=200)
                fitted.peaks <- which(!is.na(peaks[,"mu"]))
                for (p in fitted.peaks)
                {   ## plot gaussian fits
                    yg<-xcms:::gauss(tdx,peaks[p,"h"],peaks[p,"mu"],peaks[p,"sigma"])
                    lines(tdxp,yg,col='blue')
                }
            }
            Sys.sleep(sleep)
        }

    } ## if TRUE
    # return(peaks)
    return(list(peaks=peaks))
}

centWaveGB <- function(object, ppm=25, peakwidth=c(20,50), snthresh=10,
                                                    prefilter=c(3,100), mzCenterFun="wMean",
                                                    integrate=1, mzdiff=-0.001,
                                                    fitgauss=FALSE, scanrange=numeric(), noise=0, ## noise.local=TRUE,
                                                    targets=NULL,
                                                    sleep=0, verbose.columns=FALSE, ROI.list=list()) {
    if (!xcms:::isCentroided(object))
        warning("It looks like this file is in profile mode. centWave can process only centroid mode data !\n")

    mzCenterFun <- paste("mzCenter", mzCenterFun, sep=".")
    if (!exists(mzCenterFun, mode="function"))
        stop("Error: >",mzCenterFun,"< not defined ! \n")

    scanrange.old <- scanrange
    if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else
        scanrange <- range(scanrange)

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    if (!(identical(scanrange.old,scanrange)) && (length(scanrange.old) >0))
        cat("Warning: scanrange was adjusted to ",scanrange,"\n")

    # basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","intb","maxo","sn")
    basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","intb","maxo","sn")
    verbosenames <- c("egauss","mu","sigma","h","f", "dppm", "scale","scpos","scmin","scmax","lmin","lmax","pwid","baseline","sdnoise")
    targeted <- !is.null(targets)
    roi_list <- ROI.list

    if (targeted) {
      # Targeted mode: let's build a ROI list
      rois <- list()

      for (ti in 1:length(targets)) {
        target <- targets[[ti]]
        rtrange <- target$rtrange
        mzrange <- target$mzrange
        t_peakwidth <- target$peakwidth
        if (is.null(t_peakwidth)) t_peakwidth <- peakwidth
        t_snthresh <- target$snthresh
        if (is.null(t_snthresh)) t_snthresh <- snthresh
        t_integrate <- target$integrate
        if (is.null(t_integrate)) t_integrate <- integrate
        t_fitgauss <- target$fitgauss
        if (is.null(t_fitgauss)) t_fitgauss <- fitgauss

        # Convert RT to scan index
        si <- which(object@scantime >= rtrange[1] & object@scantime <= rtrange[2])
        scmin <- si[1]
        scmax <- si[length(si)]
        
        rawmz <- xcms:::rawMZ(object, mzrange=mzrange, scanrange=c(scmin, scmax))
        gz <- which(rawmz > 0)
        if (length(gz) > 0) { 
          # Define a big ROI for the entire region
          mzmin <- min(rawmz[gz])
          mzmax <- max(rawmz[gz])
          roi <- list()
          roi$mzmin <- mzmin
          roi$mzmax <- mzmax
          roi$scmin <- scmin
          roi$scmax <- scmax
          roi$peakwidth <- t_peakwidth
          roi$snthresh <- t_snthresh
          roi$integrate <- t_integrate
          roi$fitgauss <- t_fitgauss
          #cat('ROI', roi$mzmin, roi$mzmax, roi$scmin, roi$scmax, '\n')
          rois[[length(rois)+1]] <- roi
        }

        # Optionally, can also define a ROI for each consecutive scans with
        # data for this mz range. This will produce many more peaks, but more
        # closely mimic the un-targeted ROI finding behavior.
        if (FALSE) {
          cur_sc_start <- -1
          cur_min_mz <- -1
          cur_max_mz <- -1

          for (p in 1:length(rawmz)+1) { # add one more to handle when rawmz[-1] > 0
            if (p <= length(rawmz))
              mz <- rawmz[p]
            else
              mz <- 0
            if (mz > 0) {
              if (cur_sc_start < 0) { cur_sc_start <- p }
              if (cur_min_mz < 0 || mz < cur_min_mz) { cur_min_mz <- mz; }
              if (cur_max_mz < 0 || mz > cur_max_mz) { cur_max_mz <- mz; }
            }
            else {
              if (cur_sc_start >= 0) {
                roi <- list()
                roi$mzmin <- cur_min_mz
                roi$mzmax <- cur_max_mz
                roi$scmin <- si[cur_sc_start]
                roi$scmax <- si[p-1]
                roi$peakwidth <- t_peakwidth
                roi$snthresh <- t_snthresh
                roi$integrate <- t_integrate
                roi$fitgauss <- t_fitgauss
                #cat('ROI', roi$mzmin, roi$mzmax, roi$scmin, roi$scmax, '\n')
                rois[[length(rois)+1]] <- roi
                cur_sc_start <- -1
                cur_min_mz <- -1
                cur_max_mz <- -1
              }
            }
          }
        }
      }

      roi_list <- rois
    }

    ## If no ROIs are supplied then search for them.
    else if (length(roi_list) == 0) {
        cat("\n Detecting mass traces at", ppm, "ppm ... \n"); flush.console();
        roi_list <- xcms:::findmzROI(object, scanrange=scanrange,
                              dev=ppm*1e-6, minCentroids=4, prefilter=prefilter, noise=noise)
        if (length(roi_list) == 0) {
            cat("No ROIs found ! \n")
            if (verbose.columns) {
                nopeaks <- new("xcmsPeaks", matrix(nrow=0, ncol=length(basenames)+length(verbosenames)))
                # nopeaks <- matrix(nrow=0,ncol=length(basenames)+length(verbosenames))
                colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
                nopeaks <- new("xcmsPeaks", matrix(nrow=0, ncol=length(basenames)))
                # nopeaks <- matrix(nrow=0,ncol=length(basenames))
                colnames(nopeaks) <- c(basenames)
            }
            return(invisible(nopeaks))
        }
    }

    lf <- length(roi_list)

    if (!targeted && lf > 0) { # add defaults to each ROI
        for (f in 1:lf) {
            roi <- roi_list[[f]]
            roi$peakwidth <- peakwidth
            roi$snthresh <- snthresh
            roi$integrate <- integrate
            roi$fitgauss <- fitgauss
            roi_list[[f]] <- roi
        }
    }

    peaklist <- list()


    cat('\n Detecting chromatographic peaks ... \n % finished: '); lp <- -1;

    if (lf > 0) {
      for (f in 1:lf) {
        ## Show progress
        perc <- round((f/lf) * 100)
        if ((perc %% 10 == 0) && (perc != lp))
        {
            cat(perc," ",sep="");
            lp <- perc;
        }
        flush.console()

        feat <- roi_list[[f]]
        #cat('ROI', feat$mzmin, feat$mzmax, feat$scmin, feat$scmax, '\n')
        #flush.console()
        peakwidth <- feat$peakwidth
        snthresh <- feat$snthresh
        integrate <- feat$integrate
        fitgauss <- feat$fitgauss

        # peaks <- centWaveOnROI(object, scanrange, basenames, verbosenames,
        #                        feat, f, peakwidth=peakwidth, snthresh=snthresh,
        #                        mzCenterFun=mzCenterFun, integrate=integrate, fitgauss=fitgauss,
        #                        verbose.columns=verbose.columns,
        #                        targeted=targeted, sleep=sleep)

        peakstmp <- centWaveOnROI(object, scanrange, basenames, verbosenames,
                               feat, f, peakwidth=peakwidth, snthresh=snthresh,
                               mzCenterFun=mzCenterFun, integrate=integrate, fitgauss=fitgauss,
                               verbose.columns=verbose.columns,
                               targeted=targeted, sleep=sleep)
        peaks <- peakstmp$peaks
        # baseline <- peakstmp$baseline
        # sdnoise <- peakstmp$sdnoise

        if (!is.null(peaks)) {
            # peaks <- cbind(peaks, baseline=baseline, sdnoise=sdnoise)
            peaklist[[length(peaklist)+1]] <- peaks
        }
      } ## f
    }

    if (length(peaklist) == 0) {
        cat("\nNo peaks found !\n")

        if (verbose.columns) {
            nopeaks <- new("xcmsPeaks", matrix(nrow=0, ncol=length(basenames)+length(verbosenames)))
            # nopeaks <- matrix(nrow=0,ncol=length(basenames)+length(verbosenames))
            colnames(nopeaks) <- c(basenames, verbosenames)
        } else {
            nopeaks <- new("xcmsPeaks", matrix(nrow=0, ncol=length(basenames)))
            # nopeaks <- matrix(nrow=0,ncol=length(basenames))
            colnames(nopeaks) <- c(basenames)
        }

        return(invisible(nopeaks))
    }

    p <- do.call(rbind,peaklist)

    if (!verbose.columns)
        p <- p[,basenames,drop=FALSE]

    uorder <- order(p[,"into"], decreasing=TRUE)
    pm <- as.matrix(p[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE])
    uindex <- xcms:::rectUnique(pm,uorder,mzdiff,ydiff = -0.00001) ## allow adjacent peaks
    pr <- p[uindex,,drop=FALSE]
    cat("\n",dim(pr)[1]," Peaks.\n")

    peaks=invisible(new("xcmsPeaks", pr))
    # peaks=pr
}

# from cwTools.R
joinOverlappingPeaksGB <- function(td,d,otd,omz,od,scantime,scan.range,peaks,maxGaussOverlap=0.5,mzCenterFun) {

    gausspeaksidx <- which(!is.na(peaks[,"mu"]))
    Ngp <- length(gausspeaksidx)
    if (Ngp == 0)
        return(peaks)

    newpeaks <- NULL

    gpeaks <- peaks[gausspeaksidx,,drop=FALSE]
    if (dim(peaks)[1] - Ngp > 0)
        notgausspeaks <- peaks[-gausspeaksidx,,drop=FALSE]

    if (Ngp > 1) {
        comb <- which(upper.tri(matrix(0,Ngp,Ngp)),arr=TRUE)
        overlap <- rep(FALSE,dim(comb)[1])
        for (k in 1:dim(comb)[1]) {
            p1 <- comb[k,1]; p2 <- comb[k,2]
            overlap[k] <- xcms:::gaussCoverage(xlim=scan.range,h1=gpeaks[p1,"h"],mu1=gpeaks[p1,"mu"],s1=gpeaks[p1,"sigma"],
                                        h2=gpeaks[p2,"h"],mu2=gpeaks[p2,"mu"],s2=gpeaks[p2,"sigma"]) >= maxGaussOverlap
        }
    } else overlap <- FALSE

    if (any(overlap) && (Ngp > 1)) {
        jlist <- list()
        if (length(which(overlap)) > 1) {
            gm <- comb[overlap,]
            ## create list of connected components
            cc <- list()
            cc[[1]] <- gm[1,] ## copy first entry
            for (j in 2:dim(gm)[1]) { ## search for connections
                ccl <- unlist(cc)
                nl <- sapply(cc, function(x) length(x))
                ccidx <- rep(1:length(nl),nl)
                idx <- match(gm[j,],ccl)
                if (any(!is.na(idx))) { ## connection found, add to list
                    pos <- ccidx[ idx[which(!is.na(idx))[1]] ]
                    cc[[pos]] <- c(cc[[pos]],gm[j,])
                } else  ## create new list element
                    cc[[length(cc) + 1]] <- gm[j,]

            }
            ccn <- list()
            lcc <- length(cc)
            ins <- rep(FALSE,lcc)
            if (lcc > 1) {
                jcomb <- which(upper.tri(matrix(0,lcc,lcc)),arr=TRUE)
                for (j in 1:dim(jcomb)[1]) {
                    j1 <- jcomb[j,1]; j2 <- jcomb[j,2]
                    if (any(cc[[j1]] %in% cc[[j2]]))
                        ccn[[length(ccn) +1]] <- unique(c(cc[[j1]],cc[[j2]]))
                    else {
                        if (!ins[j1]) {
                            ccn[[length(ccn) +1]] <- unique(cc[[j1]])
                            ins[j1] <- TRUE
                        }
                        if (!ins[j2]) {
                            ccn[[length(ccn) +1]] <- unique(cc[[j2]])
                            ins[j2] <- TRUE
                        }
                    }
                }
            } else ccn <- cc;

            size <- sapply(ccn, function(x) length(x))
            s2idx <- which(size >= 2)

            if (length(s2idx) > 0) {
                for (j in 1:length(s2idx)) {
                    pgroup <- unique(ccn[[ s2idx[j] ]])
                    jlist[[j]] <- pgroup
                }
            } else stop('(length(s2idx) = 0) ?!?')
        } else jlist[[1]] <- comb[overlap,]

        ## join all peaks belonging to one cc
        for (j in 1:length(jlist)) {
            jidx <- jlist[[j]]
            newpeak <- gpeaks[jidx[1],,drop=FALSE]
            newmin <- min(gpeaks[jidx,"lmin"])
            newmax <- max(gpeaks[jidx,"lmax"])
            newpeak[1,"scpos"] <- -1 ## not defined after join
            newpeak[1,"scmin"] <- -1 ##    ..
            newpeak[1,"scmax"] <- -1 ##    ..
            newpeak[1,"scale"] <- -1 ##    ..

            newpeak[1,"maxo"] <- max(gpeaks[jidx,"maxo"])
            newpeak[1,"sn"]   <- max(gpeaks[jidx,"sn"])
            newpeak[1,"lmin"] <- newmin
            newpeak[1,"lmax"] <- newmax
            newpeak[1,"rtmin"] <- scantime[td[newmin]]
            newpeak[1,"rtmax"] <- scantime[td[newmax]]
            newpeak[1,"rt"] <- weighted.mean(gpeaks[jidx,"rt"],w=gpeaks[jidx,"maxo"])

            ## Re-assign m/z values
            p1 <- match(td[newmin],otd)[1]
            p2 <- match(td[newmax],otd); p2 <- p2[length(p2)]
            if (is.na(p1)) p1 <- 1
            if (is.na(p2)) p2 <- length(omz)
            mz.value <- omz[p1:p2]
            mz.int <- od[p1:p2]

            ## re-calculate m/z value for peak range
            mzmean <- do.call(mzCenterFun,list(mz=mz.value,intensity=mz.int))
            mzrange <- range(mz.value[which(mz.value > 0)])
            newpeak[1,"mz"] <- mzmean
            newpeak[1,c("mzmin","mzmax")] <- mzrange

            ## re-fit gaussian
            md <- max(d[newmin:newmax]);d1 <- d[newmin:newmax]/md;
            pgauss <- xcms:::fitGauss(td[newmin:newmax],d[newmin:newmax],pgauss = list(mu=td[newmin] + (td[newmax]-td[newmin])/2,sigma=td[newmax]-td[newmin],h=max(gpeaks[jidx,"h"])))
            if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                newpeak[1,"mu"]    <- pgauss$mu
                newpeak[1,"sigma"] <- pgauss$sigma
                newpeak[1,"h"]     <- pgauss$h
                newpeak[1,"egauss"]<- sqrt((1/length(td[newmin:newmax])) * sum(((d1-xcms:::gauss(td[newmin:newmax],pgauss$h/md,pgauss$mu,pgauss$sigma))^2)))
            } else { ## re-fit after join failed
                newpeak[1,"mu"]       <- NA
                newpeak[1,"sigma"]    <- NA
                newpeak[1,"h"]        <- NA
                newpeak[1,"egauss"]   <- NA
            }

            if (is.null(newpeaks)) newpeaks <- newpeak  else
            newpeaks <- rbind(newpeaks,newpeak)
        }
        ## add the remaining peaks
        jp <- unique(unlist(jlist))
        if (dim(peaks)[1] - length(jp) > 0)
            newpeaks <- rbind(newpeaks,gpeaks[-jp,])

    } else
        newpeaks <- gpeaks

    grt.min <- newpeaks[,"rtmin"]
    grt.max <- newpeaks[,"rtmax"]

    if (nrow(peaks) - Ngp > 0) { ## notgausspeaks
        for (k in 1:nrow(notgausspeaks)) {## here we can only check if they are completely overlapped by other peaks
            if (!any((notgausspeaks[k,"rtmin"] >= grt.min) & (notgausspeaks[k,"rtmax"] <= grt.max)))
                newpeaks <- rbind(newpeaks,notgausspeaks[k,])
        }
    }

    rownames(newpeaks) <- NULL
    newpeaks
}

mzCenter.wMean <- function(mz,intensity) {
    weighted.mean(mz, intensity)
}

mzCenter.mean <- function(mz,intensity) {
    mean(mz)
}

mzCenter.apex <- function(mz,intensity) {
    mz[which.max(intensity)]
}

mzCenter.wMeanApex3 <- function(mz,intensity) {
    iap <- which.max(intensity)
    st <- max(1,iap-1)
    en <- min(iap+1,length(mz))
    weighted.mean(mz[st:en], intensity[st:en])
}

mzCenter.meanApex3 <- function(mz,intensity) {
    iap <- which.max(intensity)
    st <- max(1,iap-1)
    en <- min(iap+1,length(mz))
    mean(mz[st:en])
}





########################################################################################
# Input: List
# Output: Standard deviation (number)
# 
# This method will calculate the standard deviation of a list (of size greater than 2)
#  in the same manner as the visual basic PEAKCALC.bas did. 
########################################################################################
standardDeviation <- function(list) {
    dblAnswer <- 0
    lengthList <- length(list[na.omit(list) > 0])

    if (lengthList > 1) {
        dblSum <- sum(list, na.rm=TRUE)
        dblMean <- dblSum / lengthList
        dblSumSqdDevs <- sum((list[list > 0] - dblMean) ^ 2, na.rm=TRUE)
        dblAnswer <- sqrt(dblSumSqdDevs / (lengthList - 1))
    }
    return(dblAnswer)
}


########################################################################################
# Metrics ##############################################################################

# ms2Call <- function (rTime, rt) {
#     return(abs(rt-mean(rTime)) / 60)
# }

ms2Call <- function (cwrt, rt) {
    return(abs(rt-cwrt) / 60)
}

noiseDetector <- function (rTime, iList, rEdgev, lEdgev) {
    rightDetect <- tail(rTime, n=1)
    leftDetect <- head(rTime, n=1)
    maxiList <- max(iList)*0.2
    noiseScore <- 0
    index <- 1

    while (leftDetect < lEdgev) {
        if (iList[match(leftDetect, rTime)] > maxiList) {
            noiseScore = noiseScore + 1
        }
        index = index + 1
        leftDetect = rTime[index]
    }

    index <- length(rTime)

    while (rightDetect > rEdgev) {
        if (iList[match(rightDetect, rTime)] > maxiList) {
            noiseScore = noiseScore + 1
        }
        index = index - 1
        rightDetect = rTime[index]
    }
    return(noiseScore)
}

normalityTests <- function (iList, rEdge, lEdge) {
    tryCatch({
        return(shapiro.test(iList[lEdge:(rEdge+1)])$p.value)
    }, error = function(ex) {
        return(1)
    })
}

derivativeChanges <- function (rTime, iList, rEdge, lEdge) {
    curr.slope = ((iList[2])-(iList[1]))/((rTime[2])-(rTime[1]))
    num.changes = 0
    curr.intensity = iList[2]
    negthreshold = .5
    posthreshold=1.5

    tryCatch({
        for (i in (lEdge+2):rEdge-1) {
            new.slope = ((iList[i+1])-(iList[i]))/((rTime[i+1])-(rTime[i]))
            new.intensity = iList[i+1]
            if (new.slope * curr.slope < 0) {
                if (((new.slope < 0) & (new.intensity <= curr.intensity/negthreshold)) |
                            ((new.slope > 0) & (new.intensity >= curr.intensity*posthreshold))) {
                    num.changes = num.changes + 1
                    curr.intensity = new.intensity
                    curr.slope = new.slope
                }
            } else if (new.slope * curr.slope > 0) {
                if (((new.slope > 0) & (new.intensity > curr.intensity)) |
                            (( new.slope < 0) & (new.intensity < curr.intensity))) {
                    curr.intensity = new.intensity
                    curr.slope = new.slope
                }
            } else if (new.slope * curr.slope == 0) {
                if (new.intensity != curr.intensity) {
                    curr.intensity = new.intensity
                    num.changes = num.changes+1
                    curr.slope = new.slope
                }
            }
        }
        return(num.changes)
    }, error = function(ex) {
        return(0)
    })
}

zeroCounts <- function (iList, rEdge, lEdge) {
    num.changes = 0
    initial.intensity = iList[lEdge]

    tryCatch({
        for (i in lEdge:rEdge) {
            if ((initial.intensity == 0) & (num.changes == 0)) {
                if (iList[i+1]!= 0) {
                    num.changes = 1
                }
            } else {
                if ((iList[i] != 0) & (iList[i+1] == 0)) {
                    num.changes = num.changes +1
                }
            }
        }
        return(num.changes)
    }, error = function(ex) {
        return(0)
    })
}

maxIntensity <- function (iList, rEdge, lEdge) {
    return(max(iList[lEdge:(rEdge+1)]))
}

rtWindow <- function (rEdgev, lEdgev) {
    return((rEdgev - lEdgev) / 60)
}

svmObjectm <- function (path) {
    dat <- read.csv(path, colClasses=c('NULL',NA,NA,NA,NA,'NULL',NA,NA,NA,NA))

    BScut <- function (x, cut) {
        x <- ifelse(x <= cut, 1, 0)
    }

    cutoff <- 3
    labels <- unlist(lapply(dat[,2], BScut, cutoff))

    gamma <- 0.15
    cost <- 0.55

    return(svm(Rank ~ ., data = dat, cost = cost, gamma = gamma))
}

metrics <- function (xcmsRaw, svmObject, rt, cwrt, rtmin, rtmax, mz, ppm) {
    EIC <- getChromatogram(xcmsRaw, mzrange=pmMz(mz, ppm), rtrange=c(rtmin-(rtmax-rtmin), rtmax+(rtmax-rtmin)))
    rTime <- EIC$rt
    iList <- EIC$intensity
    rEdge <- which.min(abs(rTime-rtmax))
    lEdge <- which.min(abs(rTime-rtmin))
    if (rEdge == lEdge) {
        if (lEdge > 1) {
            lEdge <- lEdge - 1
        }
        if (rEdge < length(rTime)) {
            rEdge <- rEdge + 1
        }
    }
    rEdgev <- rTime[rEdge]
    lEdgev <- rTime[lEdge]

    return(predict(svmObject, matrix(c(ms2Call(cwrt, rt), 
        noiseDetector(rTime, iList, rEdgev, lEdgev), 
        normalityTests(iList, rEdge, lEdge), 
        derivativeChanges(rTime, iList, rEdge, lEdge), 
        zeroCounts(iList, rEdge, lEdge), 
        maxIntensity(iList, rEdge, lEdge), 
        rtWindow(rEdgev, lEdgev)), nrow=1, ncol=7), decision.values=FALSE))
}