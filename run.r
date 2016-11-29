#!/usr/bin/Rscript
# @aly

########################################################################################
# To make this script run, one needs to run R and call this script. A shortcut can be 
#  found and used/adjusted.
########################################################################################

########################################################################################
# Names of config file (path), old directory where input files will be moved/archived,
#  and autoFill script. Old directory will be created in directory of config file.
########################################################################################
config <<- "E:/XMLinput/config.txt"
olddirName <<- "OLD"
autoFillScript <<- "autoFill.r"


########################################################################################
# Input: Config file path
# Output: Parsed config file list
# 
# This method parses the config file generated by filemaker and splits by ";". The first
#  value is the key (name), and the second value is the value of the key. This 
#  is basically a dictionary in R. File names, constants, etc. are in this config file 
#  and parsed config file list.
########################################################################################
constants <- function(f) {
    parsedConstants <- list()
    cfgFile <- file(f, open="r", encoding="UTF-16LE")
    lines <- readLines(cfgFile, warn=FALSE)
    close(cfgFile)
    for (line in lines) {
        split <- unlist(strsplit(line, ";"))
        parsedConstants[[split[1]]] <- split[2]
    }
    return(parsedConstants)
}

########################################################################################
# Input: String (message)
# Output: Nothing
# 
# This method prompts the user when it is called, providing the option to continue or 
#  quit the current running script, based on [key press] and [enter]. It also will try
#  again if an improper key was entered.
# http://stackoverflow.com/questions/10266963/moving-files-between-folders
########################################################################################
waitForInput <- function(error="") {
    cat(paste("Error occured: ", error, "\n", sep=""))
    cat("To continue, type [y]. To cancel, type [n]. Then press [enter].\n")
    key <- scan("stdin", n=1, what="")
    if (key == "y") {
        cat("Continuing job.\n")
    } else if (key == "n") {
        cat("Canceling job.\n")
        quit()
    } else {
        cat("Invalid key input, try again\n")
        waitForInput()
    }
}

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

########################################################################################
# Input: Input file path, output file path
# Output: Nothing
# 
# This method will "move" a file by renaming it and recursively create any folders it 
#  needs to do so to move the old file path to the new file path.
# http://stackoverflow.com/questions/10266963/moving-files-between-folders
########################################################################################
my.file.rename <- function(from, to) {
    todir <- dirname(to)
    if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
    file.rename(from = from,  to = to)
}

########################################################################################
# Input: Directory path
# Output: Nothing
# 
# This method will remove all files in the directory provided.
########################################################################################
removeOldFiles <- function(dir) {
    oldFiles <- list.files(dir)
    unlink(lapply(oldFiles, function(names) paste(dir, names, sep="/")), recursive = TRUE, force = FALSE)
}

########################################################################################
# Input: Nothing
# Output: Nothing
# 
# This method is the job detector, similarly to how the visual basic AutoFill was 
#  implemented. It first looks and stores the directory of the script running. Then it
#  will create an "old" directory, for storing files once they are processed. There is a
#  infinite while loop that will first check for the config file, if it exists. If it 
#  does, then it will parse it. The timestamp when the config file is stored. Files in 
#  the temporary folder where raw files will be copied over and converted is cleared 
#  out. Then it will loop through numbers 1-9, first for label free, then for SILAC. If
#  it finds a job, then it will process it by running autoFill, only if the would be 
#  output file does not exist, so previous output files would not be overwritten. If it
#  can process a job, it will run the autoFill script. Once that is done, it will copy
#  the file 1 and file 2 xml over to the old folder, under the timestamp assigned 
#  earlier. Then a copy of the output is made and also placed in that folder. This will
#  repeat until there are no numbers, and also for SILAC. When that is done, then it 
#  move other files, such as align.txt, xcal_RT_from_rep.xml, and the config file, over
#  to the old folder and the timestamp. Then variables are cleared off and the wait file
#  is deleted. Now the job detector is back to the state it an take in another config
#  file and job.
########################################################################################
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    parallel <- FALSE
    if (length(args) > 0) {
        if (toupper(args[1]) == "PARALLEL") {
            parallel <- TRUE
        }
    }
    scriptDir <- currentDir()
    parsedConstants <- NA
    waitFile <- NA
    olddir <- paste(dirname(config), "/", olddirName, "/", sep="")
    dir.create(olddir, showWarnings = FALSE)
    timeStamp <- NA

    while (TRUE) {
        while (!file.exists(config)) {
            print("Waiting for configuration file.")
            Sys.sleep(5)
        }
        if (is.na(parsedConstants)) {
            print("Configuration file found. Parsing...")
            Sys.sleep(2)
            parsedConstants <- constants(config)
            waitFile <- paste(parsedConstants[["xmldir"]], parsedConstants[["wait"]], sep="")
            timeStamp <- format(Sys.time(), "%Y%m%d/%I%M%S%p/")
            dir.create(parsedConstants[["xmldir"]], showWarnings = FALSE)
            dir.create(parsedConstants[["tempdir"]], showWarnings = FALSE)
            removeOldFiles(parsedConstants[["tempdir"]])
            dir.create(paste(olddir, timeStamp, sep=""), recursive=TRUE)
        }
        if (file.exists(waitFile)) {
            print("Jobs found.")
            Sys.sleep(2)
            print("Looking for Label Free jobs.")
            Sys.sleep(2)
            for (i in 1:9) {
                repPath1 <- paste(parsedConstants[["xmldir"]], parsedConstants[["rep"]], i, parsedConstants[["f1lf"]], sep="")
                repPath2 <- paste(parsedConstants[["xmldir"]], parsedConstants[["rep"]], i, parsedConstants[["f2lf"]], sep="")
                outPath <- paste(parsedConstants[["xmldir"]], parsedConstants[["comlf"]], i, ".txt", sep="")
                if (file.exists(repPath1) && file.exists(repPath2)) {
                    print(paste("Label Job for replicate ", i, " found.", sep=""))
                    newNameRep1 <- paste(olddir, timeStamp, basename(repPath1), sep="")
                    newNameRep2 <- paste(olddir, timeStamp, basename(repPath2), sep="")
                    Sys.sleep(1)
                    if (file.exists(outPath)) {
                        print("Cannot run job. Previous job results are still in folder.")
                        Sys.sleep(1)
                        print("Looking for next job.")
                    }
                    else {
                        print(paste("Running Label Free job for replicate ", i, ".", sep=""))
                        if (parallel) {
                            cmd <- paste("Rscript", shQuote(paste(scriptDir, autoFillScript, sep="/")), 
                                shQuote(config), i, "LABELFREE", "PARALLEL", sep=" ")
                        } else {
                            cmd <- paste("Rscript", shQuote(paste(scriptDir, autoFillScript, sep="/")), 
                                shQuote(config), i, "LABELFREE", "NOPARALLEL", sep=" ")
                        }
                        print(cmd)
                        system(cmd, show.output.on.console = TRUE)
                        if (!file.exists(outPath)) {
                            waitForInput(paste("Error occured in Autofill for replicate ", i, sep=""))
                        } else {
                            file.copy(outPath, olddir)
                            newOutputPath <- paste(olddir, parsedConstants[["comlf"]], i, ".txt", sep="")
                            my.file.rename(newOutputPath, paste(olddir, timeStamp, parsedConstants[["comlf"]], i, ".txt", sep=""))
                        }
                        my.file.rename(repPath1, newNameRep1)
                        my.file.rename(repPath2, newNameRep2)
                    }
                    if (file.exists(repPath1)) {
                        if (!file.exists(newNameRep1)) {
                            file.copy(repPath1, newNameRep1)
                        }
                    }
                    if (file.exists(repPath2)) {
                        if (!file.exists(newNameRep2)) {
                            file.copy(repPath2, newNameRep2)
                        }
                    }
                    Sys.sleep(2)
                }
                Sys.sleep(2)
            }
            print("Looking for SILAC jobs.")
            Sys.sleep(2)
            for (i in 1:9) {
                repPath1 <- paste(parsedConstants[["xmldir"]], parsedConstants[["rep"]], i, parsedConstants[["f1s"]], sep="")
                repPath2 <- paste(parsedConstants[["xmldir"]], parsedConstants[["rep"]], i, parsedConstants[["f2s"]], sep="")
                outPath <- paste(parsedConstants[["xmldir"]], parsedConstants[["coms"]], i, ".txt", sep="")
                if (file.exists(repPath1) && file.exists(repPath2)) {
                    print(paste("SILAC for replicate ", i, " found.", sep=""))
                    newNameRep1 <- paste(olddir, timeStamp, basename(repPath1), sep="")
                    newNameRep2 <- paste(olddir, timeStamp, basename(repPath2), sep="")
                    Sys.sleep(1)
                    if (file.exists(outPath)) {
                        print("Cannot run job. Previous job results are still in folder.")
                        Sys.sleep(1)
                        print("Looking for next job.")
                    }
                    else {
                        print(paste("Running SILAC job for replicate ", i, ".", sep=""))
                        if (parallel) {
                            cmd <- paste("Rscript", shQuote(paste(scriptDir, autoFillScript, sep="/")), 
                                shQuote(config), i, "SILAC", "PARALLEL", sep=" ")
                        } else {
                            cmd <- paste("Rscript", shQuote(paste(scriptDir, autoFillScript, sep="/")), 
                                shQuote(config), i, "SILAC", "NOPARALLEL", sep=" ")
                        }
                        print(cmd)
                        system(cmd, show.output.on.console = TRUE)
                        if (!file.exists(outPath)) {
                            waitForInput(paste("Error occured in Autofill for replicate ", i, sep=""))
                        } else {
                            file.copy(outPath, olddir)
                            newOutputPath <- paste(olddir, parsedConstants[["coms"]], i, ".txt", sep="")
                            my.file.rename(newOutputPath, paste(olddir, timeStamp, parsedConstants[["coms"]], i, ".txt", sep=""))
                        }
                        my.file.rename(repPath1, newNameRep1)
                        my.file.rename(repPath2, newNameRep2)
                    }
                    if (file.exists(repPath1)) {
                        if (!file.exists(newNameRep1)) {
                            file.copy(repPath1, newNameRep1)
                        }
                    }
                    if (file.exists(repPath2)) {
                        if (!file.exists(newNameRep2)) {
                            file.copy(repPath2, newNameRep2)
                        }
                    }
                    Sys.sleep(2)
                }
                Sys.sleep(2)
            }
            removeOldFiles(parsedConstants[["tempdir"]])
            newConfigName <- paste(olddir, timeStamp, basename(config), sep="")
            my.file.rename(config, newConfigName)
            newAlignName <- paste(olddir, timeStamp, parsedConstants[["align"]], sep="")
            my.file.rename(paste(parsedConstants[["xmldir"]], parsedConstants[["align"]], sep=""), newAlignName)
            newxCalName <- paste(olddir, timeStamp, parsedConstants[["xcal"]], sep="")
            my.file.rename(paste(parsedConstants[["xmldir"]], parsedConstants[["xcal"]], sep=""), newxCalName)
            file.remove(waitFile)
            parsedConstants <- NA
            timeStamp <- NA
            waitFile <- NA
        }
        else {
            print("Waiting for job.")
        }
        Sys.sleep(5)
    }
}


########################################################################################


main()