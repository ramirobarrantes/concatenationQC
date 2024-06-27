#!/usr/bin/env Rscript

options(echo=FALSE,warning=FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

sampleSheet <- args[1]
type <- args[2] # Nanopore, geomx, singular
outputFilename <- args[3]
directories <- args[4]
directories <- strsplit(directories,",")[[1]]

suppressWarnings(library(readxl))
library(tools)
library(ini)
library(stringr)
#checking files

createNanoporeSamplesheets <- function(sampleSheet,outputFilename) {
  data <- read.table(sampleSheet,header=TRUE,sep=",")
  df <- do.call(rbind,lapply(1:96, function(n) {
    barcode <- paste0("barcode",str_pad(n, 2, pad = "0"))
    fastqList <- c()
    for (j in 1:nrow(data)) {
      dir <- paste0(data[j,2],"/",barcode)
      pattern <- paste0(".*",barcode,".*.fastq.gz")
      fastqs <- list.files(dir,pattern=pattern,full.names=TRUE)
      fastqList <- c(fastqList,fastqs)
    }
    data.frame(Barcode=barcode, fastq=paste(fastqList,collapse=","))
  }))
  write.table(df,file=outputFilename,row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
}

createGeomxConcatSamplesheets <- function(sampleSheet,outputFilename) {
  data <- read.csv(sampleSheet)
  for (i in 1:nrow(data)) {
    fastqFiles <- gsub(" ","",strsplit(data[i,"FASTQ_file_location"],split=",")[[1]])
    iniFile <- data[i,"INI_Files"]
    if (!file.exists(iniFile)) {
      stop(paste0("INI file ",iniFile," does not exist"))
    }
    for (j in 1:length(fastqFiles)) {
      fastqFile <- fastqFiles[j]
      if (!file.exists(fastqFile)) {
        stop(paste0("Fastq file ",fastqFile," does not exist"))
      }
    }
  }

  geomxInput <- data.frame()
  for (i in 1:nrow(data)) {
    fastqFiles <- gsub(" ","",strsplit(data[i,"FASTQ_file_location"],split=",")[[1]])
    iniFile <- data[i,"INI_Files"]
    ini <- read.ini(iniFile)
    iniOutput <- paste0(outputDirectory,"/",basename(tools::file_path_sans_ext(iniFile)))
    if (!file.exists(iniOutput)) {
      dir.create(iniOutput)
    }
    
    aois <- names(ini$AOI_List)
    df<-do.call(rbind,lapply(1:length(aois),function(j) {
      aoi <- aois[j]
      R1 <- c()
      R2 <- c()
      for (k in 1:length(fastqFiles)) {
        fastqFile <- fastqFiles[k]
        filesR1 <- list.files(path=fastqFile,pattern=paste0(aoi,".*_R1_.*"),full.names=TRUE)
        R1 <- c(R1,filesR1)
        filesR2 <- list.files(path=fastqFile,pattern=paste0(aoi,".*_R2_.*"),full.names=TRUE)
        R2 <- c(R2,filesR2)
      }
      rbind(data.frame(ID=paste0(aoi,"_S111_L001_R1_001"),fastq=paste(R1,collapse=",")),
            data.frame(ID=paste0(aoi,"_S111_L001_R2_001"),fastq=paste(R2,collapse=",")))
      
    }))
    sampleSheet <- paste0(iniOutput,"/",basename(iniOutput),"_sampleSheet.csv")
    write.table(df,file=sampleSheet,row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
    df <- data.frame(fastqDir=paste0(iniOutput,"/rawData"),outDir=paste0(iniOutput,"/dcc"),iniFile=iniFile)
    geomxInput <- rbind(geomxInput,df)
  }
}

createSingularSamplesheets <- function(sampleSheet,directories,outputFilename) {
  tmp <- read.csv(sampleSheet)
  idx <- which(tmp[,1]=="Sample_ID")
  data <- tmp[(idx+1):nrow(tmp),]
  colnames(data) <- tmp[idx,]
  sampleIds <- unique(data$Sample_ID)
  df <- do.call(rbind,lapply(1:length(sampleIds), function(i) {
    sampleID <- sampleIds[i]
    R1 <- c()
    R2 <- c()
    for (j in 1:length(directories)) {
      dir <- directories[j]
      filesR1 <- list.files(path=dir,pattern=paste0(sampleID,".*_R1_.*.fastq.gz"),full.names=TRUE)
      R1 <- c(R1,filesR1)
      filesR2 <- list.files(path=dir,pattern=paste0(sampleID,".*_R2_.*.fastq.gz"),full.names=TRUE)
      R2 <- c(R2,filesR2)
    }
    rbind(data.frame(id=paste0(sampleID,"_R1"),fastqList=paste0('"',paste(R1,collapse=","),'"')),
          data.frame(id=paste0(sampleID,"_R2"),fastqList=paste0('"',paste(R2,collapse=","),'"')))
  }))
  write.table(df,file=outputFilename,row.names=FALSE,quote=FALSE,col.names=TRUE,sep=",")
}


type <- toupper(type)
if (type=="GEOMX") {
  concatenationSampleSheets <- createGeomxConcatSamplesheets(sampleSheet,outputFilename)
  cmd <- paste0("sbatch runConcatenation.sh ",sampleSheet," ",iniOutput)
}
if (type=="NANOPORE") {
  concatenationSampleSheets <- createNanoporeSamplesheets(sampleSheet,outputFilename)
#  cmd <- paste0("sbatch runConcatenation.sh ",sampleSheet," ",outputDirectory)
}
if (type=="SINGULAR") {
  concatenationSampleSheets <- createSingularSamplesheets(sampleSheet,directories,outputFilename)
}


message("Done.")
#system(cmd)


#cat("Once the concatenation finished, please run the following from the runNGSPipeline.sh location:\n\n ")
#for (i in 1:nrow(geomxInput)) {
#  cat(paste0("./runNGSPipeline.sh ",geomxInput$fastqDir[i]," ",geomxInput$outDir[i]," ",geomxInput$iniFile[i]),"\n")
#}
