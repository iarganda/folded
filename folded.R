
# remotes::install_github("patrickreidy/callierr")
library(jpeg)
library(rjson)
library(raster)
library(rgdal)
library(rasterVis)
library(datasets)
library(tidyr)
library(scales)
library(colorspace)
library(dplyr)
library(colourvalues)
library(circular)
# filogenia
library(phytools)
library(ggplot2)
library(geiger)
library(ape)
library(paleotree)
library(RColorBrewer)
library(data.table)
library(surface)  
'%!in%' <- function(x,y)!('%in%'(x,y))

# Supplementary data:
#
# Oscar Sanisidro, Ignacio Arganda Carreras, and Juan L. Cantalapiedra. 2022. Folded: a toolkit to describe mammalian herbivore dentition from 2D images. Methods in Ecology and Evolution.
#

options(stringsAsFactors = FALSE)
options(scipen=999)
  
# Get Working directory
  get_os <- function(){
    if (.Platform$OS.type == "windows") { 
      "win"
    } else if (Sys.info()["sysname"] == "Darwin") {
      "mac" 
    } else if (.Platform$OS.type == "unix") { 
      "unix"
    } else {
      stop("Unknown OS")
    }
  }
  
# Working directory (must be in the same folder as the R script)
  WD <- "Rhinocerotidae"
  UseTeethInFolder <- TRUE
  PlotInMM <- TRUE # plot ROI's in milimeters? (otherwise in pixels)
  PlotSpecimen <- FALSE # plot the tooth?
  
### Perform analysis ####
  
  wd_backup <- paste(getwd(),'/',WD,sep="")
  setwd(wd_backup)
  if(UseTeethInFolder){  ImageFileVector <- list.files(pattern=".jpg"); ScaleVector <- rep(22,times=length(ImageFileVector)) }

# loop for all teeth in the folder defined in WD
  DataList <- list()
  OPC_list <- list()
  for(m in 1:length(ImageFileVector)){ 
    setwd(wd_backup)
    ImageFile <- ImageFileVector[m]
    print(ImageFile)
    rasterImage <- raster(ImageFile)
    if(PlotSpecimen){plot(rasterImage, axes=F, asp=1, box=F, col=  gray.colors(4, start = 0, end = 1),  legend=FALSE, main=tools::file_path_sans_ext(ImageFile))}
 
  # Cambiar a la carpeta con resultados del diente en particular sacada de FIJI
    setwd( paste(getwd(),"/",tools::file_path_sans_ext(ImageFileVector)[m],"_analysis/",sep="") )
     
  # Sacar Listas de Archivos
    options(warn=-1)
    CsvFileNames <- list.files(pattern=".csv")
    CsvFileNames <- tools::file_path_sans_ext(CsvFileNames)
    CsvFileNames <- CsvFileNames[grepl( tools::file_path_sans_ext(ImageFile), CsvFileNames, fixed = TRUE)]
    CsvFileNamesMatrix <- t(sapply(strsplit(CsvFileNames, "_", fixed=T), "[", i = seq_len(max(sapply(strsplit(CsvFileNames, "_", fixed=T), length)))))
    CsvFileNamesMatrix[CsvFileNamesMatrix==""] <- NA
    CsvFileNamesMatrix <- CsvFileNamesMatrix[,colSums(is.na(CsvFileNamesMatrix))<nrow(CsvFileNamesMatrix)]
    if(ncol(CsvFileNamesMatrix)<4){ CsvFileNamesMatrix <- cbind(matrix(NA,ncol=4-ncol(CsvFileNamesMatrix),nrow=nrow(CsvFileNamesMatrix)),CsvFileNamesMatrix)}
    CsvSummaryNames <- CsvFileNames[which(CsvFileNamesMatrix[,ncol(CsvFileNamesMatrix)]%in%"summary")]
    CsvCoherencyDataNames <- CsvFileNames[which(CsvFileNamesMatrix[,ncol(CsvFileNamesMatrix)]%in%"complexity")]
    CsvPatchesDataNames <- CsvFileNames[which(CsvFileNamesMatrix[,ncol(CsvFileNamesMatrix)]%in%"kmeans")]
    CsvFractalDataNames <- CsvFileNames[grepl('FractalDimension', CsvFileNames, fixed = TRUE)]
    CsvEnamelDataNames <- CsvFileNames[which(CsvFileNamesMatrix[,ncol(CsvFileNamesMatrix)]%in%"enamel")]
    CsvDentineDataNames <- CsvFileNames[which(CsvFileNamesMatrix[,ncol(CsvFileNamesMatrix)]%in%"dentine")]
    CsvToothDataNames <- CsvFileNames[which(CsvFileNamesMatrix[,ncol(CsvFileNamesMatrix)]%in%"remainingTooth")]
    options(warn=0)
    SpeciesName <- strsplit(CsvCoherencyDataNames,'_')[[1]]
    SpeciesName[SpeciesName==""] <- NA
    SpeciesName <- paste(SpeciesName[!is.na(SpeciesName)][1:2],collapse="_")
    
  # Leer Orientacion (MEANorientation, MEDIANorientation, SDorientation, VARIANCEorientation) 
    FijiData <- read.csv(paste(CsvCoherencyDataNames,".csv",sep=""))
    DistanceData <- as.numeric(FijiData[,'Accumulated.distance..mm.'])
    # HAY UN ERROR EN EL SCRIPT; LA LONGITUD LINEAL DEL ESMALTE NO SE ACUMULA, SE REINICIA A CADA ROI
    DistanceData <- c(0, diff(DistanceData))
    DistanceData[DistanceData<0] <- 0
    DistanceData <- cumsum(DistanceData)
    #DistanceData <- as.numeric(FijiData[,'Accumulated.distance..mm.'])
    SpatialData <- FijiData[,c('X','Y')]
    OrientationData <- as.numeric(FijiData[,'Orientation'])
    ThicknessData <- as.numeric(FijiData[,'Thickness..pix.'])
    ThicknessData2 <- as.numeric(FijiData[,'Thickness..mm.'])
    CoherData <- as.numeric(FijiData[,'Folding..1.coherency.'])
    ROIData <- as.numeric(FijiData[,'Branch'])
  # Leer el resto de parametros
    KmeansData <- read.csv(paste(CsvPatchesDataNames,".csv",sep=""))
    FractalData <- read.csv(paste(CsvFractalDataNames,".csv",sep=""))
    ComplexityData <- read.csv(paste(CsvSummaryNames,".csv",sep=""))
    toothData <-  read.csv(paste(CsvToothDataNames,".csv",sep=""))
    enamelData <- read.csv(paste(CsvEnamelDataNames,".csv",sep=""))
    dentineData <- read.csv(paste(CsvDentineDataNames,".csv",sep=""))
     
  # Plotear el diente y los ROI
   VecUnique <- match(unique(do.call(paste, as.data.frame(cbind(FijiData[,'Skeleton'], FijiData[,'Branch'])))),do.call(paste, as.data.frame(cbind(FijiData[,'Skeleton'], FijiData[,'Branch']))))
   VecSkel <- tapply(seq_along(FijiData[,'Branch']), FijiData[,'Branch'], identity)[unique(FijiData[,'Branch'])]
    for(k in 1:length(VecSkel)){
      FijiData2plot <- FijiData[VecSkel[[k]],c('X','Y')]
      plot(FijiData2plot, type="l",xlim=c(min(FijiData[,'X']),max(FijiData[,'X'])),ylim=c(min(FijiData[,'Y']),max(FijiData[,'Y'])), asp=1, yaxt='n', xaxt='n', ann=FALSE, axes=F)
      par(new=TRUE)
    }
   plot(FijiData[VecUnique,'X'],FijiData[VecUnique,'Y'], pch=16, xlim=c(min(FijiData[,'X']),max(FijiData[,'X'])),ylim=c(min(FijiData[,'Y']),max(FijiData[,'Y'])), asp=1, yaxt='n', xaxt='n', ann=FALSE, axes=F)
   
  # Orientation
    {
    # Reescalar (el original del script da un rango de entre -1.60 y +1.60, lo cual no tiene sentido)
      OrientationData2 <- scales::rescale(OrientationData, to = c(-90,90), from = (range(c(-1.6,1.6),na.rm=T, finite=F)))
    # DEFINIR ESCALA!!!!!!!!
      OrientScale <- c( seq(from=0, to=90, by=(90-0)/180), seq(from=-90, to=0, by=(90-0)/180))
      # OrientScale <- seq(from= 0, to= 180, by=180/360); OrientScale <- OrientScale[1:(length(OrientScale)-1)]
      OrientScale <- OrientScale[1:(length(OrientScale)-1)]
      OrientScale[OrientScale==-90] <- NA
      OrientScale <- OrientScale[!is.na(OrientScale)]
      OrientScale <- cbind(OrientScale, 0:359)
    # Plotear Piechart
      OrientScaleHEX <- hex(HSV(cbind(H=OrientScale[,2],S=1,V=1)))
      OrientationDataRGB <- OrientScaleHEX
      OrientScaleHEX <- rev(c(OrientScaleHEX[181:360], rep(NA, times=360), OrientScaleHEX[1:180]))
      OrientScaleLabel <- rep("", times=720)
      OrientScaleLabelWhere <- seq(from=0, to=length(OrientScaleLabel), by=length(OrientScaleLabel)/12)
      OrientScaleLabel[OrientScaleLabelWhere] <- c("+30","+60","+90","","","","","","-90","-60","-30","0")
      par(mfrow=c(1,1),oma = c(0, 0, 0, 0))
      pie(rep(1,times=length(OrientScaleHEX)), col=OrientScaleHEX, border=NA, labels=OrientScaleLabel)
    # Ploteo
      ablinepos2 <- c(1,1+which(diff(as.numeric(as.factor(as.character(ROIData))))!=0))
      plot.new()
      # OrientationCol <- OrientationData2plot
      if(!PlotInMM){xlab.plot <- "length (pixels)"}else{xlab.plot <- "length (mm)"}
      par(mfrow=c(1,1), oma = c(8, 0, 8, 0))
      # barplot(rep(1,times=length(OrientationCol)), col=OrientationCol, space=0, border=NA, yaxt='n', xlab=xlab.plot, ylab="orientation", main="orientation (color barcode)")
      abline(v=ablinepos2, col="black", lwd=3)
      plot(DistanceData, OrientationData2, type="l", col="black", bty="n", pch="o", xlim=c(0,max(DistanceData)), xlab=xlab.plot, ylab="Orientation value", main="orientation (degrees)", lty=1)
      abline(v=ablinepos2, col="black", lwd=3)
      mtext(as.character(ImageFile), outer = TRUE, cex = 1.5)
     
    # Sumatorio Orientation
      # pasar los datos a tipo "circular"
        OrientationData.circular <- circular(as.numeric(OrientationData2), units = "degrees", zero=pi, rotation="counter", modulo = "pi")
        # plot.circular(OrientationData.circular,col=alpha("black",1), axes=T, stack=TRUE, zero=-pi/2, bins=360, shrink=4, asp=1)
        # sacar parametros circulares
        MEANorientation <- round(mean.circular(OrientationData.circular)[[1]],1) # 'mean angle' or 'preferred direction'
        MEDIANorientation <- round(median.circular(OrientationData.circular)[[1]],2)
        SDorientation <- round(sd.circular(OrientationData.circular),2) #See equation (2.33) at pag. 36 in Fisher (1993) for its definition. 
        VARIANCEorientation <- round(var.circular(OrientationData.circular)[[1]],2)
        print(summary(OrientationData.circular))
    }
    
  # Leer Anisotropy (1/Coherency; TOTALcoherency, MAXcoherency, MEANcoherency, SDcoherency, VARcoherency)  
    PlotScaleCoherency <- T
    {
      # Montar CoherencyData
        CoherencyData <- cbind('distance'=DistanceData, 'coherency'=CoherData, 'ROI'=ROIData)
      # Re escalar de 0-255 a 0-1 para plotear colores
        FijiDataHSV <- HSV(cbind(h=rep(1,times=nrow(CoherencyData)), s=rep(0,times=nrow(CoherencyData)), v=scales::rescale(CoherencyData[,2],to=c(0,1))))
        FijiDataHEX <- hex(FijiDataHSV)
      # plotear thickness
        ablinepos1 <- round(as.numeric(CoherencyData[,1])[c(1,1+which(diff(as.numeric(as.factor(CoherencyData[,'ROI'])))!=0))],0)
        ablinepos2 <- c(1,1+which(diff(as.numeric(as.factor(CoherencyData[,'ROI'])))!=0))
        if(!PlotInMM){xlab.plot <- "length (pixels)"}else{xlab.plot <- "length (mm)"}
        par(mfrow=c(2,1),oma = c(0, 0, 2, 0))
        plot(CoherencyData[,1], CoherencyData[,2], type="l", col="black", pch="o", frame.plot = FALSE, xlab=xlab.plot, ylab="anisotropy", lty=1)
        abline(v=ablinepos1, col="black", lwd=3)
        barplot(rep(1,times=nrow(CoherencyData)), col=FijiDataHEX, space=0, border=NA, yaxt='n', xlab=xlab.plot, ylab="anisotropy")
        mtext(as.character(ImageFile), outer = TRUE, cex = 1.5)
        par(mfrow=c(1,1),oma = c(0, 0, 0, 0))
      # Sumatorio thickness
        TOTALcoherency <- round(sum(as.numeric(CoherencyData[,'coherency'])),0)
        MAXcoherency <- round(max(as.numeric(CoherencyData[,'coherency'])),2)
        MEANcoherency <- round(mean(as.numeric(CoherencyData[,'coherency'])),4)
        SDcoherency <- round(sd(as.numeric(CoherencyData[,'coherency'])),2)
        VARcoherency <- round(var(as.numeric(CoherencyData[,'coherency'])),2)
        GAINcoherency <- (cumsum(as.numeric(CoherencyData[,'coherency']))-as.numeric(CoherencyData[,'coherency'])[1])[length(CoherencyData[,'coherency'])]
    }
    
  # Leer Thickness (TOTALthickness, MEANthickness, SDthickness, VARthickness)
    {
      # hacer loess con los valores de thickness (a veces da picos raros o valores de 0 en zonas muy delgadas)
      # panojilla
        ThicknessData2plot <- ThicknessData2
      # Re escalar de 0-255 a 0-1
        FijiDataHSV <- HSV(cbind(h=rep(1,times=length(ThicknessData2plot)), s=rep(0,times=length(ThicknessData2plot)), v=scales::rescale(ThicknessData2plot,to=c(0,1))))
        FijiDataHEX <- hex(FijiDataHSV)
      
      # plotear thickness
        ablinepos1 <- round(as.numeric(DistanceData)[c(1,1+which(diff(as.numeric(as.factor(CoherencyData[,'ROI'])))!=0))],0)
        ablinepos2 <- c(1,1+which(diff(as.numeric(as.factor(ROIData)))!=0))
        if(!PlotInMM){xlab.plot <- "length (pixels)"}else{xlab.plot <- "length (mm)"}
        par(mfrow=c(2,1),oma = c(0, 0, 2, 0))
        plot(DistanceData, ThicknessData2plot, type="l", col="black", pch="o", frame.plot = FALSE, xlab=xlab.plot, ylab="thickness", lty=1)
        # lines(DistanceData, ThicknessData2plot, type = "l", lty = 1, col="gray")
        abline(v=ablinepos1, col="black", lwd=3)
        barplot(rep(1,times=length(ThicknessData)), col=FijiDataHEX, space=0, border=NA, yaxt='n', xlab=xlab.plot, ylab="thickness")
        abline(v=ablinepos2, col="black", lwd=3)
        mtext(as.character(ImageFile), outer = TRUE, cex = 1.5)
        par(mfrow=c(1,1),oma = c(0, 0, 0, 0))
          
      # Grosor predominante 
        HistBreaks <- 120
        MainThickness <- ThicknessData2plot
        par(mfrow=c(1,1))
        MainThickHist <- hist(MainThickness, breaks=HistBreaks, plot=FALSE)
        # plot(MainThickHist$mids, MainThickHist$density, type="l", col="black", pch="o", frame.plot = FALSE, main="density plot: Thickness",xlab="thickness",ylab="Frequency", lty=1)
        print(paste("most frequent thickness:",round(median(ThicknessData2plot),1),sep=" "))
        print(paste("mean thickness:",round(mean(ThicknessData2plot),1),sep=" "))
        # Que orientacion tiene el grosor predominante?
        print(paste("dominant orientation of the thicker parts:",round(mean(ThicknessData2plot),1),'pixels',sep=" "))   
          
      # Sumatorio thickness
        TOTALthickness <- round(sum(as.numeric(ThicknessData2)),2)
        MEANthickness <- round(mean(as.numeric(ThicknessData2)),2)
        SDthickness <- round(sd(as.numeric(ThicknessData2)),2)
        VARthickness <- round(var(as.numeric(ThicknessData2)),2)
    }
    
  # Leer No. patches y fractalidad
    {
    TOTALpatchNo <- ComplexityData['X2D.OPC']
    TOTALfractalDim <- FractalData['D'] 
    OPC_list[[m]] <- KmeansData
    names(OPC_list)[[m]] <- ImageFile
    }
    
  # Leer Area y perimetro
    Atotal <- toothData[3]+enamelData[3]+dentineData[3]
      Abasetooth <- Atotal
    Aenamel <- enamelData[3]
    Adentine <- dentineData[3]
      Aocclusal <- Aenamel + Adentine
      AocclusRel <- round(Aocclusal/Atotal,2)
    PerimTotal <- (ComplexityData['Occlusal.enamel.length..OEL..in.pixels']*ComplexityData['Scale.bar..mm.'])/ComplexityData['Scale.bar..pixels.']
  # Medidas adicionales
    RELcoherencyAtotal <- round(TOTALcoherency/Atotal,3)
    RELcoherencyPerim <- round(TOTALcoherency/PerimTotal,3)
    RELcoherencyAocclusal <- round(TOTALcoherency/Aocclusal,3)
    RELthickness <- round(TOTALthickness/Atotal,2)
    mean_anis_x_thick <- round(mean(CoherData*ThicknessData),2)
    MEANthickness1A <- MEANthickness/PerimTotal # media de LocalThickness / Perimetro. LA ORIGINAL!!!!
    MEANthickness1B <- MEANthickness/Aocclusal # media de LocalThickness /  area total
    MEANthickness2A <- Aenamel/PerimTotal # area de esmalte dividido por perimetro esmalte
    MEANthickness2B <- Aenamel/Aocclusal # area de esmalte dividido por el area total
    RELenamel <- Aenamel/Aocclusal
  # Medidas de la literatura
    IndentIndex <- round((PerimTotal^2)/(4*pi*Atotal),2) # segun Gailer and Kaiser, 2014
    OEI <- ComplexityData['Occlusal.enamel.length..OEL..in.pixels'] # Occlusal Enamel Index (Famoso et al.). Famoso llama a PerimTotal Occlusal Enamel Length (OEL; Famoso et al.): the total length of enamel bands exposed on the occlusal surface as measured through the center of the enamel band
    EI <- ComplexityData['Enamel.Index..EI.'] # Becerra et al.
    D <- ComplexityData['Indentation.index..D.'] # indentation index
  # Meter los datos en bruto en una lista
    DataListSpecimen <- list(OrientationData, CoherData, ThicknessData)
    names(DataListSpecimen) <- c(paste(SpeciesName,"_Orientation",sep=""),paste(SpeciesName,"_Coherency",sep=""),paste(SpeciesName,"_Thickness",sep=""))
    DataList[[m]] <- DataListSpecimen

  # Vector de medidas total 
    AreaPerimData <- c(Atotal, Aenamel, RELenamel, Adentine, Abasetooth, Aocclusal, AocclusRel, PerimTotal, IndentIndex, OEI, EI) 
    names(AreaPerimData) <- c('Atotal (mm2)', 'Aenamel (mm2)', 'AenamelRel (Aenamel/Aocclusal)', 'Adentine (mm2)', 'Abasetooth (mm2)', 'Aocclusal (mm2)', 'AocclusRel (Aocc/Atotal)', 'PerimTotal (mm)','Indent index', 'OEI', 'EI(mm-1)'); print(AreaPerimData)
    AnalysisData <- c(MEANorientation, MEDIANorientation, SDorientation, VARIANCEorientation, TOTALcoherency, MEANcoherency, SDcoherency, VARcoherency, RELcoherencyAtotal, RELcoherencyAocclusal, GAINcoherency, 
                      mean_anis_x_thick, TOTALthickness, MEANthickness, MEANthickness1A, MEANthickness1B, MEANthickness2A, MEANthickness2B, SDthickness, VARthickness, RELthickness, TOTALfractalDim, TOTALpatchNo)
    names(AnalysisData) <- c("Mean orient(º)", "Median orient(º)", "SD orient", "variance orient", "total complexity", "Mean complexity", "SD complexity", "variance complexity", "Relat complexity (C/At)", "Relat complexity (C/Aocc)", "complexity gain", 
                              "Mean complexity*thick", "total thickness", "Mean thickness","Mean thickness(LocTh/Perim)", "Mean thickness(LocTh/At)", "Mean thickness(Aenam/Perim)", "Mean thickness(Aenam/At)", "SD thickness", "variance thickness", "Relat thickness", "Fractal Dimension (D)", "2D OPC"); print(AnalysisData)
     
  # Volver al directorio principal de trabajo
    setwd(wd_backup)    
  # Exportar como csv todos los parametros del diente
    
    write.csv(c(AreaPerimData, AnalysisData), file=paste(tools::file_path_sans_ext(ImageFile),".csv",sep=""))
    
  # Borrar para el siguiente loop
    rm(FijiData); rm(DistanceData); rm(SpatialData); rm(OrientationData); rm(ThicknessData); rm(ThicknessData2);
    rm(CoherData); rm(ROIData); rm(KmeansData); rm(FractalData); rm(ComplexityData); rm(toothData); rm(enamelData); rm(dentineData);
    rm(OEI); rm(EI); rm(Atotal); rm(Aenamel); rm(RELenamel); rm(Adentine); rm(Abasetooth); rm(Aocclusal); rm(AocclusRel); rm(PerimTotal);
    rm(IndentIndex); rm(MEANorientation); rm(MEDIANorientation);  rm(SDorientation); rm(VARIANCEorientation); rm(TOTALcoherency);  rm(MEANcoherency); 
    rm(SDcoherency); rm(VARcoherency);  rm(RELcoherencyAtotal); rm(RELcoherencyAocclusal); rm(GAINcoherency);  rm(mean_anis_x_thick); rm(TOTALthickness);
    rm(MEANthickness1A);  rm(MEANthickness1B); rm(MEANthickness2A); rm(MEANthickness2B);rm(SDthickness);  rm(VARthickness); rm(RELthickness); rm(TOTALfractalDim); rm(TOTALpatchNo)
   }
 
  
  
  # Other functions
  leadtrail <- function (x, rm = c("zeros", "na"), lead = c(TRUE, FALSE), trail = c(TRUE, FALSE)) {
    rm <- match.arg(rm, c("zeros", "na"))
    lead <- lead[1]
    trail <- trail[1]
    if (rm == "zeros") {
      idx <- which(x == 0)
    }
    else {
      idx <- which(is.na(x))
    }
    n <- length(x)
    l <- length(idx)
    if (lead == TRUE & l > 0) {
      if (idx[1] == 1) {
        d.idx <- diff(idx)
        loc <- which(d.idx > 1)[1]
        if (is.na(loc)) {
          loc <- l
        }
        lead.rm <- 1:loc
      }
      else {
        lead.rm <- NULL
      }
    }
    else {
      lead.rm <- NULL
    }
    if (trail == TRUE & l > 0) {
      if (tail(idx, 1) == n) {
        d.idx <- diff(rev(idx))
        loc <- which(d.idx != -1)[1]
        if (is.na(loc)) {
          loc <- l
        }
        trail.rm <- (n - loc + 1):n
      }
      else {
        trail.rm <- NULL
      }
    }
    else {
      trail.rm <- NULL
    }
    keep <- rep(TRUE, n)
    keep[lead.rm] <- FALSE
    keep[trail.rm] <- FALSE
    y <- x[keep]
    return(y)
  }
  bezierCurve <- function(x, y, n=10){
    outx <- NULL
    outy <- NULL
    
    i <- 1
    for (t in seq(0, 1, length.out=n))
    {
      b <- bez(x, y, t)
      outx[i] <- b$x
      outy[i] <- b$y
      
      i <- i+1
    }
    
    return (list(x=outx, y=outy))
  }
  bez <- function(x, y, t){
    outx <- 0
    outy <- 0
    n <- length(x)-1
    for (i in 0:n)
    {
      outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
      outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
    }
    
    return (list(x=outx, y=outy))
  }
  rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
  {
  GetDir2 <- GetDir()
  os <- get_os()
  if(os=="mac"){
    GetDir2b <- rev(gregexpr("/", GetDir2)[[1]])[1]
    GetDir <- substr(GetDir2,1,GetDir2b)
    setwd(GetDir)
  }
  if(os=="Windows" | os=="win"){
    GetDir2b <- rev(gregexpr("\\\\", GetDir2)[[1]])[1]
    cdirCommand <- paste("cd ", GetDir(), sep="")
    GetDir <- substr(GetDir2,1,GetDir2b)
    GetDir <- gsub("\\\\", "/", GetDir)
    setwd(GetDir)
  }
  getwd() 
}
 
 
  wd_backup <- paste(getwd(),'/',WD,sep="")
  setwd(wd_backup)
  AnisThreshold <- 100 # algunos giros producen blancos muy fuertes (intersecciones tb), con esto se asegura que no pasen de 100 (valor recomendado)
  if(UseTeethInFolder){  ImageFileVector <- list.files(pattern=".jpg"); ScaleVector <- rep(22,times=length(ImageFileVector)) }
  
# Loop para varios dientes (si se tiene ya una matriz, ir al ultimo paso)
 
# Leer y combinar especimenes (si ya estan hechos, no hace falta correr los pasos anteriores)
  CombineAllTeeth <- TRUE
  RemoveIndividualCSV <- TRUE
  options(scipen = 1000)
  # if(file.exists('AllTeeth.csv')){read.csv(file="AllTeeth.csv")}else{write.csv(NA, file="AllTeeth.csv")}
  if(CombineAllTeeth){
      ImageFileNamesOrig <- list.files(pattern=".jpg")
      ImageFileNames <- tools::file_path_sans_ext(ImageFileNamesOrig)
      CsvFileNamePath <- paste(getwd(), "/",ImageFileNames,".csv",sep="")
      DataMatrix <- matrix(NA, nrow=length(read.csv(CsvFileNamePath[1])),ncol=2)
      row.names(DataMatrix) <- names(read.csv(CsvFileNamePath[1]))
      if(length(CsvFileNamePath)==1){lengthCsvFileName <- 1}else{lengthCsvFileName <- 1:length(CsvFileNamePath)}
      for(i in lengthCsvFileName){ 
        DataTooth <- t(read.csv(CsvFileNamePath[i]))
        DataMatrix <- cbind(DataMatrix, DataTooth)
        }
      DataMatrix <- DataMatrix[,colSums(is.na(DataMatrix))<nrow(DataMatrix)]
      colnames(DataMatrix) <- ImageFileNames
    } # DataMatrix

# Plotear tabla combinada
  PlotVars <- FALSE
  VarsToPlot <- c(1:nrow(DataMatrix))
  par(mfrow=c(2,2), mai = c(0.8, 0.8, 0.5, 0.3))
  if(PlotVars){for(n in VarsToPlot){
  RowToPlot <- n
  plot(DataMatrix[RowToPlot,],type="l", ylab=row.names(DataMatrix)[RowToPlot], ylim=c(0,max(DataMatrix[RowToPlot,])), xlab='tooth ID', xaxt="n")
    axis(side=1, at=c(1:ncol(DataMatrix)), labels = FALSE)
    text(x=c(1:ncol(DataMatrix)), y=par()$usr[3]-0.02*(par()$usr[4]-par()$usr[3]), labels = colnames(DataMatrix), srt=45, adj=1, xpd = TRUE)
  }}
  par(mfrow=c(1,1))
  colnames(DataMatrix) <- gsub("___","", colnames(DataMatrix))
  colnames(DataMatrix) <- gsub("__","", colnames(DataMatrix))
  DataMatrix <- DataMatrix[2:nrow(DataMatrix),]

# Leer sin correr lo anterior
  # DataMatrix <- read.csv(file="AllTeeth.csv", row.names=1)
  
# Size corrected variables
  SizeCorrect <- TRUE
  # SizeTransformedVar <- c("total.complexity", "Mean.complexity","complexity.gain", "total.complexity","X2D.OPC", "log10_X2D.OPC", "log10_total.complexity", "log10_Mean.complexity", "log10_X2D.OPC") # "Mean.thickness.from.image."
  SizeTransformedVar <- c("total.complexity", "Mean.complexity","complexity.gain", "total.complexity","X2D.OPC", "Mean.thickness") # "Mean.thickness.from.image."
  print(row.names(DataMatrix))
  SizeVar <- "Atotal..mm2."
  if(SizeCorrect){
    for(k in 1:nrow(DataMatrix)){
      if(any(row.names(DataMatrix)[k]==SizeTransformedVar) && all(row.names(DataMatrix)!=paste(row.names(DataMatrix)[k],"/",SizeVar,sep=""))){
        DataMatrix <- rbind(DataMatrix,DataMatrix[k,]/DataMatrix[SizeVar,])
        row.names(DataMatrix)[nrow(DataMatrix)] <- paste(row.names(DataMatrix)[k],"/",SizeVar,sep="")
      }
    }
  }
  
# Pasar a log10 algunas de las variables
  LogTransform <- TRUE
  # LogTransformedVar <- c("Atotal..mm2.", "total.complexity", "Mean.complexity","complexity.gain", "total.complexity","X2D.OPC") # "Mean.thickness.from.image."
  LogTransformedVar <- c("Atotal..mm2.", "total.complexity", "Mean.complexity","complexity.gain", "Mean.thickness", "Mean.thickness.Aenam.At.", "total.complexity","X2D.OPC",
                         "Mean.thickness/Atotal..mm2.", "Mean.complexity/Atotal..mm2.", "X2D.OPC/Atotal..mm2.") # "Mean.thickness.from.image."
  if(LogTransform){
    for(k in 1:nrow(DataMatrix)){
      if(any(row.names(DataMatrix)[k]==LogTransformedVar) && all(row.names(DataMatrix)!=paste("log10_",row.names(DataMatrix)[k],sep=""))){
        DataMatrix <- rbind(DataMatrix,log10(DataMatrix[k,]))
        row.names(DataMatrix)[nrow(DataMatrix)] <- paste("log10_",row.names(DataMatrix)[k],sep="")
      }
     } 
  }

# Unir hipsodoncia y otras variables
  AddExternalVars <- TRUE
  if(WD=="Rhinocerotidae" && AddExternalVars){
    print(row.names(DataMatrix))
    library("readxl")
    Hypsodonty0 <- as.data.frame(read_excel("C:/Users/Oscar/Dropbox/R Scripts Files/Orientation Analysis 2D/Methods in Ecology and Evolution/TableS2.xlsx", sheet = 1))
    Hypsodonty0 <- Hypsodonty0[,c('image file name','hypsodonty')]
    Hypsodonty <- Hypsodonty0[,2]+1
    names(Hypsodonty) <- gsub(" ", "_", Hypsodonty0[,1])
    DataMatrix <- rbind(DataMatrix, hypsodonty=Hypsodonty)
    #plot(DataMatrix["log10_X2D.OPC/Atotal..mm2.",], Hypsodonty,pch=16,xlab="log10_X2D.OPC/Atotal..mm2.")
    #text(DataMatrix["log10_X2D.OPC/Atotal..mm2.",], Hypsodonty, labels = colnames(DataMatrix),cex = 0.6, pos = 4, col = "black")
    #plot(DataMatrix["log10_Mean.complexity/Atotal..mm2.",], Hypsodonty,pch=16,xlab="log10_Mean.complexity/Atotal..mm2.")
    #text(DataMatrix["log10_Mean.complexity/Atotal..mm2.",], Hypsodonty, labels = colnames(DataMatrix),cex = 0.6, pos = 4, col = "black")
  }
  
  print(DataMatrix)
  write.csv(t(DataMatrix), file="AllTeeth.csv")
  if(RemoveIndividualCSV){ file.remove(CsvFileNamePath) }
 
  
#______________TEST________________________________________________________________________________________________________________________________________________________________________________
  
    DataMatrix2 <- as.data.frame(DataMatrix)
    Var_normality_test <- matrix(NA,nrow=nrow(DataMatrix2),ncol=3); row.names(Var_normality_test) <- row.names(DataMatrix2); colnames(Var_normality_test) <- c("Norm_ShapW","Norm_p.value","Norm?")
    DataMatrixPhyloNames <- apply( matrix(unlist(strsplit(colnames(DataMatrix2), "_")), ncol=3,byrow=T)[,1:2],1,paste,collapse="_")
    colnames(DataMatrix2) <- DataMatrixPhyloNames
    par(mfrow=c(2,2), mai = c(0.8, 0.8, 0.5, 0.3))
    plot_diagnosis <- FALSE # Use plot_diagnosis function to create diagnostic plots
    for(n in 1:nrow(as.data.frame(DataMatrix2))){ 
      if(plot_diagnosis){qqPlot(DataMatrix2[n,], pch=20, ylab=row.names(DataMatrix2)[n])}
      shap.results <- shapiro.test(as.numeric(DataMatrix2[n,]))
      print(row.names(DataMatrix2)[n])
      print(shap.results) # p value non-sign: normal distribution of residuals
      if(shap.results$p.value > 0.05){
        print("p value non-sign: normal distribution of residuals")
        Var_normality_test[n,'Norm?'] <- 'normal'
      }
      if(shap.results$p.value <= 0.05){
        print("p value sign: non-normal distribution of residuals")
        Var_normality_test[n,'Norm?'] <- 'not normal'
      }
      Var_normality_test[n,'Norm_ShapW'] <- round(shap.results$statistic,3)
      Var_normality_test[n,'Norm_p.value'] <- round(shap.results$p.value,3)
      write.csv(Var_normality_test, 'normality_tests.csv')
  }
    print(Var_normality_test)
    # Leer y datar el arbol
    FAD_file <- "C:/Users/Oscar/Dropbox/R Scripts Files/Orientation Analysis 2D/Rhinocerotidae/Refs/FAD.csv"
    FAD_names <- c("FAD","LAD") 
    Tree_file <- "C:/Users/Oscar/Dropbox/R Scripts Files/Orientation Analysis 2D/Rhinocerotidae/Refs/Elasmotheriini+Rhinocerotina.nex"
    
# Distribucion normal de las variables por separado
    {
        # leer los FAD
          FADs <- read.csv(FAD_file, sep=";")
          FADs[,1] <- gsub("__","", FADs[,1])
          ages <- FADs
          row.names(ages) <- FADs[,1]
          ages[,1] <- ages[,2]
          colnames(ages) <- FAD_names
        # Leer el arbol 
          Tree <- read.nexus(Tree_file)
        # sacar nombre de las especies y comprobar que todas estan en el CSV y viceversa
          Tree <- drop.tip(Tree, Tree$tip.label[!Tree$tip.label%in%colnames(DataMatrix2)], trim.internal = TRUE)
          TreeTipLabel <- Tree$tip.label
        # volver a seleccionar ages y fad/lad
          ages <- ages[row.names(ages) %in% TreeTipLabel,]
          FADs <- FADs[FADs[,1] %in% TreeTipLabel,]
        # Datar el arbol con todas las especies
          ages <- ages[rownames(ages) %in% Tree$tip.label,]
          Tree <- timePaleoPhy(ladderize(Tree),ages,type="equal",vartime=2,ntrees=1,add.term=T,randres=TRUE, inc.term.adj = FALSE, dateTreatment = "randObs", plot=FALSE)
    }

# Regresiones lineales
    library("car")
    library("ggpubr")
    library("caper")
    library("nlme")
    {
    # OLS / GLS
      print(row.names(DataMatrix))
      SelectedVars <- c(
        "log10_Mean.complexity","Atotal..mm2.",
        "log10_X2D.OPC", "Atotal..mm2.",
        "log10_Mean.thickness","Atotal..mm2.",
        "log10_Mean.complexity", "log10_Mean.thickness",
        "log10_X2D.OPC", "log10_Mean.thickness",
        "log10_X2D.OPC", "log10_Mean.complexity",
        "log10_Mean.complexity/Atotal..mm2.", "log10_Mean.thickness/Atotal..mm2.",
        "log10_X2D.OPC/Atotal..mm2.", "log10_Mean.thickness/Atotal..mm2.",
        "log10_Mean.complexity","hypsodonty",
        "log10_X2D.OPC", "hypsodonty",
        "log10_Mean.thickness","hypsodonty"
        )
      SelectedVars_Matrix <- matrix(SelectedVars, ncol=2, byrow=T)
      corr_data <- as.data.frame(t(DataMatrix2))
      corr_data$Species <- rownames(corr_data)
          
        # Generalized least squares (OLS) regression
          Corr_table <- matrix(NA,nrow=nrow(SelectedVars_Matrix),ncol=5)
          row.names(Corr_table) <- apply(SelectedVars_Matrix,1,paste,collapse=" ~ ")
          colnames(Corr_table) <- c("intercept","slope", "P", "RSE", "R2")
          for(m in 1:nrow(SelectedVars_Matrix)){
            print( SelectedVars_Matrix[m,])
            SelectedVars <-  SelectedVars_Matrix[m,]
            a <- as.numeric(corr_data[,SelectedVars[1]])
            b <- as.numeric(corr_data[,SelectedVars[2]])
            # Ordinary least squares (OLS) regression
            ols_model <- lm(a ~ b, data=corr_data) #  If "REML" the model is fit by maximizing the restricted log-likelihood. If "ML" the log-likelihood is maximized. Defaults to "REML".
            ols_summary <- summary(ols_model)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "P"] <- round(ols_summary$coefficients[2,'Pr(>|t|)'],3)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "RSE"] <- round(ols_summary$sigma,3)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "R2"] <- round(ols_summary$r.squared, 3)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "intercept"] <-  round(ols_model$coefficients['(Intercept)'],3)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "slope"] <- unlist(round(unlist(ols_model$coefficients['b']),3))
          }
          Corr_table_OLS <- cbind(rep("OLS",times=nrow(Corr_table)),Corr_table)
          
        # PGLS
          Corr_table <- matrix(NA,nrow=nrow(SelectedVars_Matrix),ncol=5)
          row.names(Corr_table) <- apply(SelectedVars_Matrix,1,paste,collapse=" ~ ")
          colnames(Corr_table) <- c("intercept","slope", "P", "RSE", "R2")
          for(m in 1:nrow(SelectedVars_Matrix)){
            print(SelectedVars_Matrix[m,])
            SelectedVars <-  SelectedVars_Matrix[m,]
            comp_data <- comparative.data(phy=Tree, data=corr_data, names.col=Species, vcv=TRUE, vcv.dim=3, na.omit=F)
            pgls_model <- caper::pgls(get(SelectedVars[1]) ~ get(SelectedVars[2]), data=comp_data)
            pgls_summary <- summary(pgls_model)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "P"] <- round(pgls_summary$coefficients[,'Pr(>|t|)'][2][[1]],2)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "RSE"] <-  round(pgls_model$sterr[2][[1]],2)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "intercept"] <- round(pgls_model$model$coef[1,1][[1]],2)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "slope"] <- round(pgls_summary$coefficients[,'Estimate'][2][[1]],2)
            Corr_table[paste(SelectedVars,collapse=" ~ "), "R2"] <- round(pgls_summary$r.squared, 3)
          }
          Corr_table_PGLS <- cbind(rep("PGLS",times=nrow(Corr_table)),Corr_table)
          
      # Juntar OLS y PGLS
        if(exists("Corr_table_PGLS") && exists("Corr_table_OLS")){ 
          Corr_table <- rbind(Corr_table_OLS, Corr_table_PGLS)
          row.names(Corr_table) <- gsub("[.]", " ", row.names(Corr_table)); row.names(Corr_table) <- gsub("  ", " ", row.names(Corr_table)); row.names(Corr_table) <- gsub("X", "", row.names(Corr_table)); 
          Corr_table[,'P'][Corr_table[,'P']<0.1] <- paste(Corr_table[,'P'],"*",sep="")[Corr_table[,'P']<0.1]
          Corr_table[,'P'][Corr_table[,'P']<0.05] <-  paste(Corr_table[,'P'],"*",sep="")[Corr_table[,'P']<0.05] 
          Corr_table[,'P'][Corr_table[,'P']<0.001] <- paste(Corr_table[,'P'],"*",sep="")[Corr_table[,'P']<0.001] 
          Corr_table[,'P'][Corr_table[,'P']<0.001] <- '<0.001***'
    }
        # limpiar tabla general
        Corr_table <- Corr_table[order(Corr_table[,1]),]
        rownamesCorr_table <- row.names(Corr_table)
        rownamesCorr_table <- gsub("_"," ",rownamesCorr_table); rownamesCorr_table <- gsub("X","",rownamesCorr_table)
        rownamesCorr_table <- gsub("Mean complexity","EF",rownamesCorr_table);
        rownamesCorr_table <- gsub("Atotal mm2","OTA",rownamesCorr_table)
        row.names(Corr_table) <- rownamesCorr_table
        colnames(Corr_table)[1] <- 'method'
        write.csv(Corr_table,"Corr_table.csv")
        print(Corr_table)
    }
    
    save(corr_data, file='corr_data_hypsodonty.R')
    
# Test 2A: resolution test
  if(WD == "Test-resolution"|| WD == "Test-resolution2"){
   # Sacar los valores sin hacer resampling por tamano (con las resoluciones originales) para comparar
    if(WD == "Test-resolution"){ 
       DataMatrix2plotREF <- read.csv(paste(getwd(),"2/AllTeeth.csv",sep=""))
       DataMatrix2plotREF <- DataMatrix2plotREF[2:nrow(DataMatrix2plotREF),]
       row.names(DataMatrix2plotREF) <- DataMatrix2plotREF[,1]; DataMatrix2plotREF <- DataMatrix2plotREF[,2:ncol(DataMatrix2plotREF)]
       DataMatrix2plotREF <- DataMatrix2plotREF[,ColOrder]
       for(k in 1:nrow(DataMatrix2plotREF)){
         if(any(row.names(DataMatrix2plotREF)[k]==LogTransformedVar)){
           DataMatrix2plotREF[k,] <- log10(DataMatrix2plotREF[k,])
           row.names(DataMatrix2plotREF)[k] <- paste("log10 ",row.names(DataMatrix2plotREF)[k],sep="")
       }}
    }
   # Plotear tabla combinada
     DataMatrix2plot <- DataMatrix
     ColNames <- sapply(strsplit(colnames(DataMatrix2plot),"_"), "[[", length(strsplit(colnames(DataMatrix2plot),"_")[[1]]))
     ColOrder <- as.numeric(gsub("px","",ColNames))
     ColOrder <- order(ColOrder); ColNames <- ColNames[ColOrder]
     DataMatrix2plot <- DataMatrix2plot[,ColOrder]
     par(mfrow=c(3,3), mai = c(0.8, 0.8, 0.5, 0.3))
   # Activar para test de resolucion
     for(n in 1:nrow(DataMatrix2plot)){
       ImageRes <- as.numeric(gsub("px","",matrix(unlist(strsplit(colnames(DataMatrix2plot),"_")),ncol=4,byrow=T)[,4]))
       maxY <- max(unlist(c(DataMatrix2plot[n,],DataMatrix2plotREF[n,])))
       minY <- min(unlist(c(DataMatrix2plot[n,],DataMatrix2plotREF[n,])))
       plot(ImageRes, DataMatrix2plot[n,], ylab=rownames(DataMatrix2plot)[n], type="l", ylim=c(minY,maxY), pch=16, xlab='Resolution', main=SpeciesName)
       if(WD == "Test-resolution"){ lines(ImageRes, DataMatrix2plotREF[n,], pch=16, type="l",col="red") }
   }
     par(mfrow=c(1,1))
   # test estadisticos
     # summary(t(DataMatrix))
   
   # Robust t-test of Yuen [Yuen and Dixon (1973), Yuen (1974)] in place of the standard Welch t test(t.test stats). The yuen t test makes no assumption of normality.
   # The function computes the robust TOST for a sample of paired differences or for two samples. 
   #library(equivalence)
   #rtost(x, y = NULL, alpha = 0.05, epsilon = 0.31, tr = 0.15) # tr the proportion (percent/100) of the data set to be trimmed
  
  }
 
# Test 2B: resolution test
  if(WD == "Test-resolution3"){
   # Sacar los valores sin hacer resampling por tamano (con las resoluciones originales) para comparar
   # Plotear tabla combinada
   DataMatrix2plot <- DataMatrix
   ColNames <- matrix(unlist(strsplit(colnames(DataMatrix2plot),"_")),ncol=4,byrow=T)
   ColNames <- cbind(paste(ColNames[,1], ColNames[,2]), ColNames[,4])
   ColNames[,2] <- gsub("px","", ColNames[,2])
   var2plotVector <- c('X2D.OPC','total.complexity')
   par(mfrow=c(1,2), mai = c(0.8, 0.8, 0.5, 0.3))
   for(n in 1:length(var2plotVector)){
     var2plot <- var2plotVector[n]
     var2plot2 <- gsub("log10","",var2plot)
     Matrix2plot <- matrix(DataMatrix2plot[var2plot,],ncol=length(unique(ColNames[,1])),byrow=F)
     row.names(Matrix2plot) <- unique(ColNames[,2]); colnames(Matrix2plot) <- unique(ColNames[,1])
     Matrix2plot <- Matrix2plot[order(as.numeric(row.names(Matrix2plot))),]
     # plot
     xlimPlot <- as.numeric(c(50,max(as.numeric(row.names(Matrix2plot)))))
     ylimPlot <- c(1.5,ceiling(max(Matrix2plot)))
     cl <- rep("red",times=ncol(Matrix2plot))
     plot(0,0, xlim = xlimPlot, ylim = ylimPlot, ylab=var2plot2, xlab="resolution", type = "l")
     for (i in 1:ncol(Matrix2plot)){par(new=TRUE); plot(as.numeric(row.names(Matrix2plot)), Matrix2plot[,i], axes = FALSE, ylab=NA, xlab=NA,
          xlim = xlimPlot, ylim = ylimPlot, col = cl[i], type = 'o', pch=16) }
   }
   par(mfrow=c(1,1))
   }
 
# Test 3: wear
  if(WD == "Test wear1"){
   PlotAll <- TRUE
   Hsteps <- as.numeric(as.numeric(sapply(strsplit(colnames(DataMatrix),"_"), "[[", 1)))
   Data2plot <- DataMatrix[,order(Hsteps)]
   Hsteps <- sort(Hsteps)
   # Plotear tabla combinada
   VarsToPlot <- c(1:nrow(Data2plot))
   par(mfrow=c(3,3), mai = c(0.8, 0.8, 0.5, 0.3))
   if(PlotAll){for(n in VarsToPlot){
     RowToPlot <- n
     plot(Hsteps, Data2plot[RowToPlot,],type="l", ylab=row.names(Data2plot)[RowToPlot], ylim=c(min(Data2plot[RowToPlot,]),max(Data2plot[RowToPlot,])), xlab="tooth's crown height")
     # axis(side=1, at=Hsteps, labels = Hsteps)
    }}
   par(mfrow=c(1,1))
   # Plotear variables seleccionadas
   VarsToPlot <- c("Mean.thickness.from.image.", "log10 total.complexity", "log10 X2D.OPC")
   par(mfrow=c(3,1), mai = c(0.8, 0.8, 0.5, 0.3))
   for(n in 1:length(VarsToPlot)){
     RowToPlot <- VarsToPlot[n]
     plot(Hsteps, Data2plot[RowToPlot,],type="l", ylab=RowToPlot, ylim=c(max(Data2plot[RowToPlot,]),max(Data2plot[RowToPlot,])), xlab="tooth's crown height")
     dev.copy(pdf, file = paste("wear_test_1.pdf",sep=""), width=12, height=4); dev.off()
     # axis(side=1, at=Hsteps, labels = Hsteps)
   }
   par(mfrow=c(1,1))
   # OLS
   # OLS Regression in R programming is a type of statistical technique, that is used for modeling. It is also used for the analysis of linear relationships
   # between a response variable. If the relationship between the two variables is linear, a straight line can be drawn to model their relationship.
   
  }
   
# Test 3B: wear2
  if(WD == "Test wear1" || WD == "Test wear2"){
    PlotAll <- TRUE
    colnamesDataMatrix <- colnames(DataMatrix)
    Hsteps <- as.numeric(as.numeric(sapply(strsplit(colnamesDataMatrix,"_"), "[[", 1)))
    Data2plot <- DataMatrix[,order(Hsteps)]
    Hsteps <- rev(sort(Hsteps))
    # Plotear tabla combinada
    VarsToPlot <- c(1:nrow(Data2plot))
    par(mfrow=c(3,3), mai = c(0.8, 0.8, 0.5, 0.3))
    if(PlotAll){for(n in VarsToPlot){
      RowToPlot <- n
      plot(Hsteps, Data2plot[RowToPlot,],type="l", ylab=row.names(Data2plot)[RowToPlot], ylim=c(min(Data2plot[RowToPlot,]),max(Data2plot[RowToPlot,])), xlab="tooth's crown height")
      # axis(side=1, at=Hsteps, labels = Hsteps)
    }}
    par(mfrow=c(1,1))
    # Plotear variables seleccionadas
    print(row.names(Data2plot))
    VarsToPlot <- c("Mean.thickness", "total.complexity", "X2D.OPC")
    par(mfrow=c(1,3), mai = c(0.8, 0.8, 0.5, 0.3))
    for(n in 1:length(VarsToPlot)){
      RowToPlot <- VarsToPlot[n]
      maxInterval <- Data2plot[RowToPlot,]+5
      minInterval <- Data2plot[RowToPlot,]-5
      plot(Hsteps, Data2plot[RowToPlot,],type="o", ylab=RowToPlot, ylim=c(0,max(Data2plot[RowToPlot,])), 
           xlim=c(max(Hsteps),min(Hsteps)), xlab="tooth's crown height")
      # polygon(c(Hsteps,rev(Hsteps)),c(maxInterval,rev(minInterval)),col = alpha("grey75",0.5), border = FALSE)
      # axis(side=1, at=Hsteps, labels = Hsteps)
    }
    par(mfrow=c(1,1))
    dev.copy(pdf, file = paste("wear_test_2.pdf",sep=""), width=12, height=4); dev.off()
    
    # OLS
    # OLS Regression in R programming is a type of statistical technique, that is used for modeling. It is also used for the analysis of linear relationships
    # between a response variable. If the relationship between the two variables is linear, a straight line can be drawn to model their relationship.
  }
  
# Test 5A: rotation test
  if(WD == "Test-rotation"){
    # Plotear tabla combinada
    DataMatrix2plot <- DataMatrix
    ColNames <- sapply(strsplit(colnames(DataMatrix2plot),"_"), "[[", length(strsplit(colnames(DataMatrix2plot),"_")[[1]]))
    ColOrder <- as.numeric(gsub("px","",ColNames))
    ColOrder <- order(ColOrder); ColNames <- ColNames[ColOrder]
    DataMatrix2plot <- DataMatrix2plot[,ColOrder]
    par(mfrow=c(2,2), mai = c(0.8, 0.8, 0.5, 0.3))
    # Activar para test de resolucion
    for(n in 1:nrow(DataMatrix2plot)){
      plot(DataMatrix2plot[n,],type="l", ylab=rownames(DataMatrix2plot)[n], ylim=c(min(DataMatrix2plot[n,]),max(DataMatrix2plot[n,])), xlab='Resolution', xaxt="n")
      axis(side=1, at=c(1:ncol(DataMatrix2plot)), labels = FALSE)
      text(x=c(1:ncol(DataMatrix2plot)), y=par()$usr[3]-0.02*(par()$usr[4]-par()$usr[3]), labels = ColNames, srt=45, adj=1, xpd = TRUE)
    }
    par(mfrow=c(1,1))
    # test estadisticos
    # Robust t-test of Yuen [Yuen and Dixon (1973), Yuen (1974)] in place of the standard Welch t test(t.test stats). The yuen t test makes no assumption of normality.
    
    # The function computes the robust TOST for a sample of paired differences or for two samples. 
    # library(equivalence)
    # rtost(x, y = NULL, alpha = 0.05, epsilon = 0.31, tr = 0.15) # tr the proportion (percent/100) of the data set to be trimmed
    
  }

# Test 5B: rotation test2
  if(WD == "Test-rotation3"){
   # Sacar los valores sin hacer resampling por tamano (con las resoluciones originales) para comparar
   # Plotear tabla combinada
   DataMatrix2plot <- DataMatrix
   ColNames <- matrix(unlist(strsplit(colnames(DataMatrix2plot),"_")),ncol=4,byrow=T)
   ColNames <- cbind(paste(ColNames[,1], ColNames[,2]), ColNames[,4])
   ColNames[,2] <- gsub("px","", ColNames[,2])
   var2plotVector <- c('X2D.OPC','total.complexity')
   par(new=TRUE); par(mfrow=c(1,2), mai = c(0.8, 0.8, 0.5, 0.3))
   for(n in 1:length(var2plotVector)){
     var2plot <- var2plotVector[n]
     var2plot2 <- gsub("log10","_",var2plot)
     Matrix2plot <- matrix(DataMatrix2plot[var2plot,],ncol=length(unique(ColNames[,1])),byrow=F)
     row.names(Matrix2plot) <- unique(ColNames[,2]); colnames(Matrix2plot) <- unique(ColNames[,1])
     Matrix2plot <- Matrix2plot[order(as.numeric(row.names(Matrix2plot))),]
     # plot
     xlimPlot <- as.numeric(c(0,max(as.numeric(row.names(Matrix2plot)))))
     ylimPlot <- c(1.5,ceiling(max(Matrix2plot)))
     cl <- rep("red",times=ncol(Matrix2plot))
     plot(0,0, xlim = xlimPlot, ylim = ylimPlot, ylab=var2plot2, xlab="resolution", type = "l")
     for (i in 1:ncol(Matrix2plot)){par(new=TRUE); plot(as.numeric(row.names(Matrix2plot)), Matrix2plot[,i], axes = FALSE, ylab=NA, xlab=NA,
                                                        xlim = xlimPlot, ylim = ylimPlot, col = cl[i], type = 'o', pch=16) }
     dev.copy(pdf, file = paste(paste(var2plotVector,collapse="_"),'_rotation.pdf',sep="_"), width=10, height=7); dev.off()
   }
   par(mfrow=c(1,1))
 }
  
# Test 6: Suidae
  if(WD == "Tooth_Pigs"){
   # Sacar los valores sin hacer resampling por tamano (con las resoluciones originales) para comparar
   # Plotear tabla combinada
   print(row.names(DataMatrix))
   DataMatrix2plot <- DataMatrix
   ColNames <- matrix(unlist(strsplit(colnames(DataMatrix2plot),"_")),ncol=4,byrow=T)
   ColNames <- cbind(paste(ColNames[,1], ColNames[,2]),ColNames[,3], ColNames[,4])
   PigNames <- c("Babyrousa babyrussa", "Hylochoerus meinertzhageni","Phacochoerus africanus","Sus scrofa","Potamochoerus porcus")
   TeethOrder <- c("p2","p3","p4","m1","m2","m3","combined"); TeethOrder <- cbind(TeethOrder,1:length(TeethOrder))
   # var2plotVector <- c('log10 X2D.OPC','log10 total.complexity')
   # var2plotVector <- c('log10 X2D.OPC','log10 Mean.complexity')
   var2plotVector <- c('X2D.OPC','total.complexity')
   TraitPalette_microwear <- c("#778c93","#0d5888","#831e11","#776813","#ceab0d")
   
  # load microwear data (Souron et al 2014).  Asfc / epLsar (×10^-3) / HAsfc-81
   MicroWear_complexity_mean <-    c(NA,1.71,2.12,3.57,3.03)
   MicroWear_complexity_sd <-      c(NA,0.89,1.80,2.84,3.18)
   MicroWear_anisotropy_mean <-    c(NA,4.81,4.09,2.81,3.34)
   MicroWear_anisotropy_sd <-      c(NA,2.97,2.39,1.33,1.84)
   MicroWear_heterogeneity_mean <- c(NA,1.34,0.72,0.97,0.98)
   MicroWear_heterogeneity_sd <-   c(NA,0.72,1.03,0.86,1.04)
   MicroWear <- t(rbind(MicroWear_complexity_mean,MicroWear_complexity_sd, MicroWear_anisotropy_mean, MicroWear_anisotropy_sd, MicroWear_heterogeneity_mean, MicroWear_heterogeneity_sd))
   MicroWear <- as.data.frame(cbind(species=PigNames, MicroWear))
   MicroWear <- MicroWear[complete.cases(MicroWear), ]
   
   # var2plotVector <- c("Mean.complexity.thick","Mean.thickness.from.image.")
   par(mfrow=c(1,2), mai = c(0.8, 0.8, 0.5, 0.3))
   for(n in 1:length(var2plotVector)){
     var2plot <- var2plotVector[n]
     DataMatrix2plot2 <- DataMatrix2plot[,ColNames[,1]%in%PigNames]
     Matrix2plot <- as.data.frame(cbind(ColNames[,1][ColNames[,1]%in%PigNames],ColNames[,2][ColNames[,1]%in%PigNames],DataMatrix2plot2[var2plot,]))
     Matrix2plot <- Matrix2plot[order(Matrix2plot[,1],match(Matrix2plot[,2],TeethOrder[,1])),]
     Matrix2plot2 <- matrix(NA,nrow=nrow(TeethOrder),ncol=length(unique(Matrix2plot[,1]))); row.names(Matrix2plot2) <- TeethOrder[,1]; colnames(Matrix2plot2) <- unique(Matrix2plot[,1])
     for(j in 1:ncol(Matrix2plot2)){ Matrix2plot2[row.names(Matrix2plot2)%in%Matrix2plot[Matrix2plot[,1]==unique(Matrix2plot[,1])[j],2],j] <- Matrix2plot[Matrix2plot[,1]==unique(Matrix2plot[,1])[j],3] }
     Matrix2plot2 <- Matrix2plot2[rowSums(is.na(Matrix2plot2)) != ncol(Matrix2plot2),]
     class(Matrix2plot2)<-"numeric"
     # Matrix2plot2 <- rbind(Matrix2plot2,total=colSums(Matrix2plot2,na.rm=TRUE))
     # plot
     xlimPlot <- as.numeric(c(1,length(unique(TeethOrder[,1]))))
     ylimPlot <- c(min(as.numeric(Matrix2plot[,3])),max(as.numeric(Matrix2plot[,3])))
     plot(0,0, xlim = xlimPlot, ylim = ylimPlot, ylab=var2plot, xlab="tooth position", xaxt = "n", type = "l"); axis(1, at=1:nrow(Matrix2plot2), labels=row.names(Matrix2plot2))
     cl <- TraitPalette_microwear
     for (i in 1:ncol(Matrix2plot2)){
       par(new=TRUE); plot(c(1:nrow(Matrix2plot2))[!is.na(Matrix2plot2[,i])],Matrix2plot2[,i][!is.na(Matrix2plot2[,i])], axes=F, ylab=NA, xlab=NA, xlim = xlimPlot, ylim = ylimPlot, col = cl[i], type = 'o', pch=16)
       legend("topleft", legend=colnames(Matrix2plot2), col=cl, pch=16, bty = "n") # optional legend
      }
   }
   par(mfrow=c(1,1))
   # plotear microwear
   cl <- TraitPalette_microwear[2:length(TraitPalette_microwear)]
   par(mfrow=c(1,3), mai = c(0.8, 0.8, 0.5, 0.3))
   for(n in c(2,4,6)){
     MicroWear_selected <- MicroWear[,c(1,n,n+1)]
     MicroWear_selected <- cbind(MicroWear_selected, sdmax = as.numeric(MicroWear_selected[,2])+as.numeric(MicroWear_selected[,3]), 
                                 sdmin = as.numeric(MicroWear_selected[,2])-as.numeric(MicroWear_selected[,3]))
     MicroWear_ylim <- c(min(as.numeric(MicroWear_selected[,2])-as.numeric(MicroWear_selected[,3])),
                         max(as.numeric(MicroWear_selected[,2])+as.numeric(MicroWear_selected[,3])))
     X <- seq(from=0,to=(0.1*nrow(MicroWear_selected)),by=0.1)[2:length(seq(from=0,to=(0.1*nrow(MicroWear_selected)),by=0.1))]
     plot(X,MicroWear_selected[,2], pch=16,xlab="",ylab="",xaxt="n", ylim=MicroWear_ylim, cex=2, col=cl)
     X2 <- rbind(cbind( rep(X[1],times=2),rep(X[2],times=2),rep(X[3],times=2),rep(X[4],times=2)), c(NA,NA,NA,NA))
     Y2 <- rbind(t(MicroWear_selected[,c(4:5)]),rep(NA, times=ncol(t(MicroWear_selected[,c(4:5)]))))
     lines(X2,Y2,col=cl, lwd=2)
     legend("topleft", legend=MicroWear_selected[,1], col=cl, pch=16, bty = "n") # optional legend
     mtext(colnames(MicroWear_selected)[2])
   }
   axis(side=1,at=1:3,labels=MicroWear$color)
   par(mfrow=c(1,1))
 }
 
# Test 7: Prolagus
  if(WD == "Tooth-Prolagus"){
    print(row.names(DataMatrix))
    var2plotVector <- c('log10_X2D.OPC','log10_Mean.complexity')
    # var2plotVector <- c('X2D.OPC','Mean.complexity')
    # var2plotVector <- c('log10_X2D.OPC/Atotal..mm2.','log10_Mean.complexity/Atotal..mm2.')
    # var2plotVector <- c('Mean.thickness', 'X2D.OPC','Mean.complexity')
    # Sacar los valores sin hacer resampling por tamano (con las resoluciones originales) para comparar
    # Plotear tabla combinada
    DataMatrix2plot <- as.data.frame(t(DataMatrix[var2plotVector,]))
    ColNames <- matrix(unlist(strsplit(row.names(DataMatrix2plot),"_")),ncol=3,byrow=T)
    ColNames <- apply(ColNames[,1:2],1,paste,collapse="_")
    DataMatrix2plot <- cbind(species=ColNames, DataMatrix2plot)
    Hypsodonty_names <- c("Prolagus_tobieni", "Prolagus_vasconiensis", "Prolagus_oeningensis", "Prolagus_michauxi", "Prolagus_affdepereti", "Prolagus_depereti", 
                          "Prolagus_apricenicus", "Prolagus_imperialis", "Prolagus_figaro", "Prolagus_cfsardus","Prolagus_sardus", "Prolagus_crusafonti","Prolagus_major",
                          "Prolagus_calpensis")
    Hypsodonty_value <- c(4.31,4.57,4.69,4.53,NA,NA,4.73,4.83,4.53,6.33,6.33, NA,NA,NA)
    Insular_value <- c("C","C","C","C","I","I","I","I","I","I","I","C","C","C")
    Hypsodonty <- cbind(Hypsodonty_names, Hypsodonty_value, Insular_value)
    
    # unir hipsodoncia e insularidad a las especies
    attachedMatrix <- matrix(NA,ncol=2,nrow=nrow(DataMatrix2plot))
    for(j in 1:nrow(attachedMatrix)){ attachedMatrix[j,] <- Hypsodonty[match(DataMatrix2plot[j,1],Hypsodonty[,1]),2:3]}
    DataMatrix2plot <- cbind(DataMatrix2plot, attachedMatrix) 
    TraitPalette_microwear <- c("#778c93","#0d5888","#831e11","#776813","#ceab0d", "", "", "", "")
    # plotear microwear
    cl <- TraitPalette_microwear[2:length(TraitPalette_microwear)]
    par(mfrow=c(1,1), mai = c(0.8, 0.8, 0.5, 0.3))
    plot(DataMatrix2plot[,2],DataMatrix2plot[,3],xlab=colnames(DataMatrix2plot)[2],ylab=colnames(DataMatrix2plot)[3],pch=16)
    text(DataMatrix2plot[,2],DataMatrix2plot[,3], labels=DataMatrix2plot[,1],cex=0.9)
    dev.copy(pdf, file = paste("Prolagus.pdf",sep=""), width=10, height=7); dev.off()
    # Body size
    SelectedVars1 <- "Atotal..mm2."
    SelectedVars2 <- c("log10_Mean.complexity", "log10_X2D.OPC")
    corMethod <- "kendall" 
    DataMatrix2plot <- as.data.frame(t(DataMatrix[c(SelectedVars1,SelectedVars2),]))
    for(h in 1:length(SelectedVars2)){
      windows(height = 7, width = 7)
      ggscatterPlot1 <- ggscatter(DataMatrix2plot, x=SelectedVars1, y=SelectedVars2[h], add = "reg.line", conf.int = TRUE, cor.coef = TRUE, main="")
      print(ggscatterPlot1)
      dev.copy(pdf, file = paste(gsub("[.]","_",SelectedVars1),"_",gsub("[.]","_",SelectedVars2[h]),".pdf",sep=""), width=7, height=7); dev.off()
    }
  }
  
# Test 8: proboscideans
  if(WD == "Tooth-proboscideans"){
    print(row.names(DataMatrix))
    var1 <- 'Mean.thickness'
    var2 <- c('X2D.OPC','Mean.complexity')
    OPC_data <- OPC_list[[1]]
    # Sacar los valores sin hacer resampling por tamano (con las resoluciones originales) para comparar
    DataMatrix2 <- as.data.frame(DataMatrix)
    DataMatrixPhyloNames <- apply( matrix(unlist(strsplit(colnames(DataMatrix2), "_")), ncol=3,byrow=T)[,1:2],1,paste,collapse="_")
    colnames(DataMatrix2) <- DataMatrixPhyloNames
    var2plotVector <- c(var1, var2)
    DataMatrix2 <- t(DataMatrix2[var2plotVector,])
    # plotear
    par(mfrow=c(1,1), mai = c(0.8, 0.8, 0.5, 0.3))
    for(n in 1:2){#length(var2)){
      plot(DataMatrix2[,var1],DataMatrix2[,var2[n]], ylab=var2[n], xlab=var1, col=rainbow(nrow(DataMatrix2)), pch=16, cex=2)
      # text(DataMatrix2[,var1],DataMatrix2[,var2[n]], labels=row.names(DataMatrix2),cex=0.9)
      legend("bottomright", legend=row.names(DataMatrix2), col=rainbow(nrow(DataMatrix2)), pch=16, bty = "n") # optional legend
      dev.copy(pdf, file = paste(gsub('[.]',"_",var1),"_vs_",gsub('[.]',"_",var2[n]),'.pdf',sep=""), width=10, height=7); dev.off()
    }
    # plotear las areas de los OPC
    barplot_col <- c("#E92210","#E0AF12","#BEDC14","#2AD316","#18C89F","#1A8FC2","#0A31AE","#2F0277")
    barplot_col <- c("#2F0277", "#0A31AE", "#1A8FC2", "#18C89F", "#2AD316","#BEDC14","#E0AF12","#E92210")
    par(mfrow=c(2,2), mai = c(0.8, 0.8, 0.5, 0.3))
    for(n in 1:length(OPC_list)){
      barplot(OPC_list[[n]][,'Area..pixels.2.'], names.arg=OPC_list[[n]][,'Mean'], col=barplot_col, main=names(OPC_list)[n], space=0, border=NA)
    }
  }
  
# Test 9: 2D OPC vs 3D OPC
  if(WD == "Test OPC"){
    print(row.names(DataMatrix))
    var2plotVector <- c('X2D.OPC')
    # var2plotVector <- c('X2D.OPC','Mean.complexity')
    # var2plotVector <- c('log10_X2D.OPC/Atotal..mm2.','log10_Mean.complexity/Atotal..mm2.')
    # var2plotVector <- c('Mean.thickness', 'X2D.OPC','Mean.complexity')
    # Sacar los valores sin hacer resampling por tamano (con las resoluciones originales) para comparar
    OPC3D_names <- c("Eohippus_sp", "Miohippus_sp", "Parahippus_sp", "Merychippus_sp", "Hipparion_sp", "Equus_sp")
    OPC3D_values <- c(35, 45, 25, 41, 31, 28,
                      88.625, 99.25, 96.5, 110.25, 112.625, 109.5,
                      150, 160, 167, 188, 180, 173,
                      193.25, 204.875, 223.25, 254.375, 271.375, 259.75,
                      239.25, 262.25, 292.375, 317.625, 339.625, 357)
    
    OPC3D <- t(rbind(species=OPC3D_names, matrix(OPC3D_values,nrow=5,byrow=T)))
    row.names(OPC3D) <- OPC3D[,1]
    OPC3D <- OPC3D[,2:ncol(OPC3D)]
    colnames(OPC3D) <- c('3DOPC 25rows','3DOPC 50rows','3DOPC 75rows','3DOPC 100rows','3DOPC 125rows')
      
    # unir OPC 3D y 2D
    DataMatrix2 <- as.data.frame(DataMatrix)
    DataMatrixPhyloNames <- apply( matrix(unlist(strsplit(colnames(DataMatrix2), "_")), ncol=3,byrow=T)[,1:2],1,paste,collapse="_")
    colnames(DataMatrix2) <- DataMatrixPhyloNames
    DataMatrix2 <- t(DataMatrix2[var2plotVector,])
    DataMatrix2 <- DataMatrix2[match(row.names(OPC3D), row.names(DataMatrix2)),]
    DataMatrix2 <- cbind(OPC3D, DataMatrix2)
    colnames(DataMatrix2)[ncol(DataMatrix2)] <- '2D OPC (1000px)'
    TraitPalette_OPC <- c("skyblue","skyblue1","skyblue2","skyblue3","skyblue4","black")
    # plotear
    cl <- TraitPalette_OPC
    par(mfrow=c(1,1), mai = c(0.8, 0.8, 0.5, 0.3))
    xlimPlot <- as.numeric(c(1,nrow(DataMatrix2)))
    ylimPlot <- c(min(as.numeric(DataMatrix2)),max(as.numeric(DataMatrix2)))
    plot(0,0, xlim = xlimPlot, ylim = ylimPlot, ylab='2D/3D OPC', xlab="species", xaxt = "n", type = "l"); axis(1, at=1:nrow(DataMatrix2), labels=row.names(DataMatrix2))
    for (i in 1:ncol(DataMatrix2)){
      par(new=TRUE); plot(c(1:nrow(DataMatrix2))[!is.na(DataMatrix2[,i])],DataMatrix2[,i][!is.na(DataMatrix2[,i])], axes=F, ylab=NA, xlab=NA, xlim = xlimPlot, ylim = ylimPlot, col = cl[i], type = 'o', pch=16)
      legend("topleft", legend=colnames(DataMatrix2), col=cl, pch=16, bty = "n") # optional legend
    }
    text(DataMatrix2[,],DataMatrix2[,], labels=DataMatrix2[,1],cex=0.9)
    # dev.copy(pdf, file = paste("Prolagus.pdf",sep=""), width=10, height=7); dev.off()
   
  }
  
  
#____PHYLOGENETIC MODELLING________________________________________________________________________________________________________________________________________________________________________________

  PhenogramPlot <- TRUE
  PGLS <- FALSE
  OUWIE <- TRUE
  PhenoSpreadLabels <- FALSE
  
  print(row.names(DataMatrix))
  # VarsToPlot <- c("Mean.thickness","log10_Mean.complexity","log10_X2D.OPC")
  # VarsToPlot <- c("Mean.thickness.LocTh.Perim.","log10 Mean.complexity","log10 X2D.OPC")
  VarsToPlot <- c("log10_Mean.complexity/Atotal..mm2.","log10_X2D.OPC/Atotal..mm2.")
  # VarsToPlot <- c("Mean.thickness/Atotal..mm2.","log10_Mean.complexity/Atotal..mm2.","log10_X2D.OPC/Atotal..mm2.")
  # VarsToPlot <- c("Mean.thickness/Atotal..mm2.","log10_Mean.complexity/Atotal..mm2.","log10_X2D.OPC/Atotal..mm2.")
  # VarsToPlot <- c("Mean.thickness","log10_Mean.complexity","log10_X2D.OPC")
  # VarsToPlot <- c("Mean.thickness.Aenam.At.","Mean.complexity/Atotal..mm2.","X2D.OPC/Atotal..mm2.")
  # VarsToPlot <- c("Mean.thickness.LocTh.Perim.","Mean.complexity","X2D.OPC")
  
  GrousToPlot <- c("Elasmotheriini", "Rhinocerotini")
  cols_pheno <- c("orangered","chartreuse3","gold2")
  FAD_file <- "C:/Users/Oscar/Dropbox/R Scripts Files/Orientation Analysis 2D/Rhinocerotidae/Refs/FAD.csv"
  FAD_names <- c("FAD","LAD") 
  Tree_file <- "C:/Users/Oscar/Dropbox/R Scripts Files/Orientation Analysis 2D/Rhinocerotidae/Refs/Elasmotheriini+Rhinocerotina.nex"
  # VarIndep <- "Mean.thickness"
  # VarIndep <- "log10_Mean.thickness/Atotal..mm2."
  VarIndep <- "log10_Mean.thickness/Atotal..mm2."
  s_aic_threshold <- 6
  { 
  # explorar datos
    library(GGally)
    ggpaisData <- DataMatrix[VarsToPlot,]
    ggpairs(as.data.frame(t(DataMatrix[VarsToPlot,])))
  
  # Transponer matriz de datos
     DataMatrixPhylo <- t(DataMatrix)
     DataMatrixPhyloNames <- apply( matrix(unlist(strsplit(row.names(DataMatrixPhylo), "_")), ncol=3,byrow=T)[,1:2],1,paste,collapse="_")
     row.names(DataMatrixPhylo) <- DataMatrixPhyloNames
     
    # leer los FAD
     FADs <- read.csv(FAD_file, sep=";")
     FADs[,1] <- gsub("__","", FADs[,1])
     ages <- FADs
     row.names(ages) <- FADs[,1]
     ages[,1] <- ages[,2]
     colnames(ages) <- FAD_names
       
    # Leer el arbol 
     Tree <- read.nexus(Tree_file)
     # sacar nombre de las especies y comprobar que todas estan en el CSV y viceversa
     Tree <- drop.tip(Tree, Tree$tip.label[!Tree$tip.label%in%row.names(DataMatrixPhylo)], trim.internal = TRUE)
     TreeTipLabel <- Tree$tip.label
     # volver a seleccionar ages y fad/lad
     ages <- ages[row.names(ages) %in% TreeTipLabel,]
     FADs <- FADs[FADs[,1] %in% TreeTipLabel,]
     
    # sacar la lista de Elasmotheriina y Rhinocerotina
     tribe <- rep("Elasmotheriini",times=nrow(DataMatrixPhylo))
     tribe[1:23] <- "Rhinocerotini"
     names(tribe) <- TreeTipLabel
     names(cols_pheno)<- c(unique(tribe),"NA")
    
   # Select taxa
     SelectedTaxa <- names(tribe[tribe%in%GrousToPlot])
     DataMatrixTree <-  DataMatrixPhylo[row.names(DataMatrixPhylo) %in% SelectedTaxa,]
     CSVnames <- row.names(DataMatrixTree)
     Tree2 <- drop.tip(Tree, setdiff(Tree$tip.label, SelectedTaxa));
     
  # Datar el arbol con todas las especies
     ages <- ages[rownames(ages) %in% Tree2$tip.label,]
     Tree2 <- timePaleoPhy(ladderize(Tree2),ages,type="equal",vartime=2,ntrees=1,add.term=T,randres=TRUE, inc.term.adj = FALSE, dateTreatment = "randObs", plot=FALSE)
     par(mfrow=c(1,1))
     plotTree(Tree2,pts=F,node.numbers=T)
     
  # colores por tribu
    Tree_tribe <- Tree2
    if((length(unique(tribe[tribe %in% GrousToPlot])))>1){
      Tree_tribe <- paintSubTree(Tree_tribe,node=48,state="NA")
      Tree_tribe <- paintSubTree(Tree_tribe,node=49,state="Elasmotheriini")
      Tree_tribe <- paintSubTree(Tree_tribe,node=71,state="Rhinocerotini")
      # Tree_tribe <- make.simmap(Tree2, tribe[tribe %in% GrousToPlot], model="SYM", nsim=1) 
    }else{Tree_tribe <- Tree2}
    
  # Fenograma: plotear las variables en la filogenia
    PhenogramPlot <- TRUE
    if(PhenogramPlot){
   # colores para los caracteres
    library(wesanderson)
    # TraitPalette <- "YlOrBr"
    TraitPalette <- wes_palette("Zissou1", 100, type = "continuous")
   
   # plotear fenograma
     for(k in 1:length(VarsToPlot)){
       windows(height = 7, width = 12)
       par(mfrow=c(1,2))
        # fenograma
          PhenoYlabel <- VarsToPlot[k]
          PhenoData <- DataMatrixTree[,PhenoYlabel]
          phenogram(Tree_tribe, PhenoData, ylab=PhenoYlabel, spread.labels=PhenoSpreadLabels, fsize=0.8, ftype="i", 
                    colors = cols_pheno, main=VarsToPlot[k], xlab="time (Ma)", axes=list("time","trait"))
        # plotear caracteres
          objTree <- contMap(Tree2, PhenoData, plot=FALSE)
          objTree_palette <- setMap(objTree, TraitPalette)
        # cambiar colores
          n_cols <- n_distinct(PhenoData)
          plot(objTree_palette, fsize=c(0.8,1), outline=FALSE, lwd=c(3,7), mtext=VarsToPlot[k])
          mylims <- par("usr")
          # plotear hipsodoncia en caso de ser complejidad o 2D OPC
          if(AddExternalVars && any(VarsToPlot[k]==c("log10 Mean.complexity","log10 X2D.OPC"))){
            Ordered_hypsodonty <- Hypsodonty[objTree_palette$tree$tip.label %in% names(Hypsodonty)]
            Ordered_hypsodonty <- rev(Ordered_hypsodonty[match(objTree_palette$tree$tip.label,names(Ordered_hypsodonty))])
            TraitPalette2 <- wes_palette(name="Zissou1", max(Ordered_hypsodonty), type = "continuous")
            Ordered_hypsodonty_cols <- rep(NA,times=length(Ordered_hypsodonty))
            for(j in 1:length(Ordered_hypsodonty_cols)){Ordered_hypsodonty_cols[j] <- TraitPalette2[Ordered_hypsodonty[j]]}
            par(oma=c(0,65,0,6))
            points(rep(1,times=length(objTree_palette$tree$tip.label)), c(1:length(objTree_palette$tree$tip.label)), 
                   pch=15, cex=1.5, col=alpha(Ordered_hypsodonty_cols,1))
          }
          pdf_name <- paste("Phenogram_",gsub("[.]","_",VarsToPlot[k]),".pdf",sep="")
          pdf_name <- paste("Phenogram_",gsub("/","_div_",VarsToPlot[k]),".pdf",sep="")
          dev.copy(pdf, file = pdf_name, width=12, height=7); dev.off()
      }
     par(mfrow=c(1,1))
  }
   
# PHYLOGENETIC MODELLING: SURFACE
   library(surface) 
   library(OUwie)
   library(extendedSurface) # devtools::install_github("mongiardino/extendedSurface") # DECIR NO A ACTUALIZAR PACKAGES ACCESORIOS!!!!
   # library(devtools); install_github("https://github.com/Schraiber/pulsR")
   library(pulsR)
   library(pals)#colores
   phylomorphospace_labels <- "horizontal" # "horizontal o off
   VarDep <- VarsToPlot[VarsToPlot!=VarIndep]
   cols2_list <- list()
   t_list <- list()
   OptimaCoords_list <- list()
   
   # plotear filomorfoespacio
   for(h in VarDep){
   # seleccionar variables
     SelectedVars <- c(VarIndep,h)
   
   # Dvidir por clados?
     splitClades <- FALSE
     if(splitClades){n <- 1 }else{n=c(1,2)}
     GroupSurf <- unique(tribe)[n]
     SelectedTaxa <- names(tribe[tribe%in%GroupSurf])
     
   # arbolico
     # nombrar todos los nodos, requisito de surface
     surface_tree <- Tree2
     surface_tree$node.label <- c(1:(length(surface_tree$tip.label)-1))
     
   # datos
     surface_data <- as.data.frame(DataMatrixTree)
     surface_data <- surface_data[,SelectedVars,drop=FALSE]
     colnames(surface_data) <- gsub("/",".",colnames(surface_data)); colnames(surface_data) <- gsub("[.]","_",colnames(surface_data)); colnames(surface_data) <- gsub(" ","_",colnames(surface_data))
       
   # run surface
      s <- runSurface(surface_tree, surface_data, error_skip=F, plotaic=F, verbose=T, aic_threshold = s_aic_threshold)
   
   # resultados del mejor modelo de surface
      sum <-	surfaceSummary(s$bwd)
      kk <- length(s$bwd)
      kk2 <- length(s$fwd)
      hansenfit <- s$bwd[[kk]]
      fit <- hansenfit$fit[[1]]
      otree <- as(fit, "data.frame")
      otree <- data.frame(otree, shifts = rep(NA, length(otree$nodes)))
      otree$shifts[match(names(hansenfit$savedshifts), otree$nodes)] <- 1:length(hansenfit$savedshifts)
      otree2 <- otree[match(c(surface_tree$tip.label, surface_tree$node.label), otree$labels), ]
      otree2 <- otree2[surface_tree$edge[, 2], ]
      OptimaCoords <- sum$theta
      OptimaCoords_list[[h]] <- OptimaCoords
      branch_regimes <- list()
      for(i in 1:length(s$bwd[[kk]]$fit)){branch_regimes[[i]] <- s$bwd[[kk]]$fit[[i]]@regimes$regs}; names(branch_regimes) <- SelectedVars
      RegimesVector <- levels(unlist(branch_regimes))
      for(i in 1:length(s$bwd[[kk]]$fit)){print(s$bwd[[kk]]$fit[[i]])}
      print(sum$n_regimes)
      sum$n_regimes["c"]/(sum$n_regimes["k"]-1) # based on shifts
      sum$n_regimes["kprime"]/(sum$n_regimes["k"]) # based on regimes
         # k (the number of regime shifts, counting the basal regime as 1)
         # kprime, (the number of regimes, some of which may be reached by multiple shifts)
         # deltak (k-kprime, a measure of convergence)
         # c (the number of shifts to convergent regimes, another measure of convergence)
         # kprime_conv (the number of convergent regimes shifted to multiple times)
         # kprime_nonconv (the number of nonconvergent regimes only shifted to once)
     # plot regimes
       cols <- brewer.pal(n = length(levels(branch_regimes[[1]])), name = "Dark2") 
       branch.cols <- cols[as.numeric(factor(otree2[, 5]))]
      # cols <- wes_palette("Darjeeling1", n = length(levels(branch_regimes[[1]])), type = "continuous")
       if(ncol(surface_data)==1){maxtraits <- 1}else{maxtraits <- 2}
       if(ncol(surface_data)==1){
         par(mfrow=c(1,2), mai=c(1,1,1,1), mar=c(6,2,4,2))
         surfaceTreePlot(surface_tree, hansenfit=s$bwd[[kk]], labelshifts = T, edge.width=3, cex=0.75, convcol=T)
         surfaceTraitPlot(surface_data, s$bwd[[kk]], whattraits = c(1,maxtraits), optellipses=T, plotoptima=T, cols = NULL, cex.opt = 3)
         mtext(colnames(surface_data), side = 3, line = - 2, outer = TRUE)
         par(mfrow=c(1,1))
       }
       if(ncol(surface_data)>1){
         t <- surface_tree
         windows(height = 7, width = 12)
         par(mfrow=c(1,2), mar=c(0,0,0,0))
         surfaceTreePlot(t, hansenfit=s$bwd[[kk]], labelshifts = T, cols=cols, edge.width=3, cex=.7, convcol=T)
         # hacer filomoroespacio
         t <- paintSubTree(t, 48, "a", anc.state="a", stem=F)
         equiv.int.nodes <- as.numeric(otree$labels[1:surface_tree$Nnode])
         equiv.int.nodes <- equiv.int.nodes+length(surface_tree$tip.label)
         equiv.nodes <- c(equiv.int.nodes, length(surface_tree$tip.label):1)
         shift.nodes.equiv <- equiv.nodes[as.numeric(names(sum$shifts))]
         ordered.shift.nodes <- c(sort(shift.nodes.equiv[shift.nodes.equiv>length(t$tip.label)]),sort(shift.nodes.equiv[shift.nodes.equiv<=length(t$tip.label)]))
         ordered.shift.regimes <- sum$shifts[match(ordered.shift.nodes, shift.nodes.equiv)]
         h <- gsub("[.]","_",h); h <- gsub(" ","_",h); h <- gsub("/","_",h, fixed=TRUE)
         regime <- as.vector(unlist(s$bwd[[kk]]$fit[[h]]@regimes))
         for(i in 2:length(ordered.shift.nodes)){
           nnod <- ordered.shift.nodes[i]
           regime.to.paint <- ordered.shift.regimes[i]
           t <- paintSubTree(t, nnod, regime.to.paint, anc.state="a", stem=T)
         }
         surf_nodes <- unlist(s$bwd[[kk]]$fit[[1]]@regimes)
         names(surf_nodes) <- s$bwd[[kk]]$fit[[1]]@nodelabels
         cols2 <- setNames(cols,levels(surf_nodes))
         cols2_list[[h]] <- cols2
         t_list[[h]] <- t
         lim_plot_x <- c(min(min(OptimaCoords[,1]),min(surface_data[,1])),max(max(OptimaCoords[,1]), max(surface_data[,1])))
         # lim_plot_y <- c(min(min(OptimaCoords[,2]),min(surface_data[,2])),max(max(OptimaCoords[,2]), max(surface_data[,2])))
         lim_plot_y <- c(min(surface_data[,2]),max(surface_data[,2]))
         par(mar=c(6,5,3,1))
         # surfaceTraitPlot(surface_data, s$bwd[[kk]], whattraits = c(1,maxtraits), optellipses=T, plotoptima=T, cols = alpha(cols,0.5), pch=16, cex.opt = 3) # 
         plot(NULL, xlim=lim_plot_x, ylim=lim_plot_y, xlab=colnames(surface_data)[1], ylab=colnames(surface_data)[2], bty="n")
         points(OptimaCoords, pch=16, cex=11, col=alpha(cols2,0.5))
         phylomorphospace(t, as.matrix(surface_data), node.size=c(0,1.5), pch=16, lwd=2, cex=0.6, colors=cols2,
                          xlim=lim_plot_x, ylim=lim_plot_y, label=phylomorphospace_labels, node.by.map=TRUE, add=T)
         mtext(paste(gsub("_"," ",colnames(surface_data)),collapse=" vs "), side = 3, line = - 2, outer = TRUE)
         dev.copy(pdf, file = paste(paste(colnames(surface_data),collapse=" vs "),".pdf",sep=""), width=12, height=7); dev.off()
         
         # plotear lo mismo pero pintandolo por grupos
         windows(height = 7, width = 12)
         par(mfrow=c(1,2), mar=c(0,0,0,0))
         surfaceTreePlot(surface_tree, hansenfit=s$bwd[[kk]], labelshifts = T, cols=cols, edge.width=3, cex=.7, convcol=T)
         par(mar=c(6,5,3,1))
         plot(NULL, xlim=lim_plot_x, ylim=lim_plot_y, xlab=colnames(surface_data)[1], ylab=colnames(surface_data)[2], bty="n")
         points(OptimaCoords, pch=16, cex=11, col=alpha(cols2,0.5))
         phylomorphospace(Tree_tribe, as.matrix(surface_data), node.size=c(0,1.5), pch=16, lwd=2, cex=0.6, colors=cols_pheno,
                          xlim=lim_plot_x, ylim=lim_plot_y, label=phylomorphospace_labels, node.by.map=TRUE, add=T)
         mtext(paste(gsub("_"," ",colnames(surface_data)),collapse=" vs "), side = 3, line = - 2, outer = TRUE)
         dev.copy(pdf, file = paste(paste(colnames(surface_data),collapse=" vs "),"_groups.pdf",sep=""), width=12, height=7); dev.off()
         par(mfrow=c(1,1))
       }
    }  
  }
  
  # Export data
  # save(DataMatrixTree, VarsToPlot, Tree2, file="caper_data2.RData")
  
  # plotear fenogramas con el color de los regimenes evolutivos
    DataMatrixTree2 <- DataMatrixTree
    colnames(DataMatrixTree2) <- gsub("[.]","_", colnames(DataMatrixTree2))
    colnames(DataMatrixTree2) <- gsub("/","_", colnames(DataMatrixTree2))
    windows(height = 7, width = 12)
    par(mfrow=c(1,3))
    for(h in 1:length(t_list)){
      PhenoData <- DataMatrixTree2[,names(t_list)[h]]
      phenogram(t_list[[h]], PhenoData, ylab=PhenoYlabel, spread.labels=PhenoSpreadLabels, fsize=0.8, ftype="i", 
          colors = cols2_list[[h]], main=names(t_list)[h], xlab="time (Ma)", axes=list("time","trait"))
  }
    # sacar el arbol de la variable indepdentiente
    phenogram(t_list[[1]], DataMatrixTree2[,gsub("[.]","_",'log10_Mean_thickness_Atotal__mm2_')], ylab=PhenoYlabel, spread.labels=PhenoSpreadLabels, fsize=0.8, ftype="i", 
          colors = cols2_list[[1]], main='log10_Mean_thickness_Atotal__mm2_', xlab="time (Ma)", axes=list("time","trait"))
    dev.copy(pdf, file = paste("phenograms_regimes.pdf",sep=""), width=12, height=7); dev.off()
    par(mfrow=c(1,1))
  
    
    
    