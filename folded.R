# Supplementary data:
#
# Oscar Sanisidro, Ignacio Arganda Carreras, and Juan L. Cantalapiedra. 2022. Folded: a toolkit to describe mammalian herbivore dentition from 2D images. Methods in Ecology and Evolution.
#
require(ape)
require(paleotree)
require(stats)
require(caper)
require(circular) 
require(colorspace) 
require(nlme)
require(surface) 
require(OUwie)
require(pulsR)
require(pals)
require(phytools)
require(dplyr)
require(extendedSurface) 
# devtools::install_github("mongiardino/extendedSurface")
'%!in%' <- function(x,y)!('%in%'(x,y))
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
options(stringsAsFactors = FALSE)
options(scipen=999)
  
# Get Working directory (SELECT ONE)

  # Set the working folder. Depending on the analysis, the options would be:
  # WD <- "Rhinocerotidae dataset" # the example dataset from the main text
  # WD <- "Resolution dataset" # resolution test with all the Rhinocerotidae included in the main manuscript. Takes a while to process...
  # WD <- "Rotation dataset" # rotation test with all the Rhinocerotidae included in the main manuscript.
  # WD <- "Wear dataset 1" # wear test of Diceros bicornis in the Supplementary dataset
  # WD <- "Wear dataset 2" # wear test of Ceratotherium simum in the Supplementary dataset
  # WD <- "Suidae dataset" # Suidae in the Supplementary dataset
  # WD <- "Lagomorpha dataset" # Prolagus in the Supplementary dataset
  # WD <- "Proboscidea dataset" # Proboscidea in the Supplementary dataset
  # WD <- "2D vs 3D OPC dataset" # Comparison between 2D and 3D OPC in the Supplementary dataset
  
  working_folder <- "C:/Users/Rhinocerotidae dataset/Rhinocerotidae dataset" # ADD HERE THE PATH TO THE FOLDER WHERE ALL FOLDED IMAGES AND OUTPUT ARE. E.G.: 'C:/Users/Rhinocerotidae dataset'
  
  PlotInMM <- TRUE # plot ROI's in milimeters? (otherwise in pixels)
  PlotSpecimen <- FALSE # plot the tooth?
  ImageFormat <- 'jpg' # image format (no period). CAUTION: do not place other images in the same format than the data sample in the working folder
  
### Start the analysis ####
  wd_backup <- working_folder
  setwd(working_folder)
  ImageFileVector <- list.files(pattern=paste(".",ImageFormat,sep=""))

# loop for all teeth included in the WD folder defined above
  DataList <- list()
  OPC_list <- list()
  DataMatrix <- matrix(NA,ncol=length(ImageFileVector),nrow=29)
  colnames(DataMatrix) <- tools::file_path_sans_ext(ImageFileVector)
  for(m in 1:length(ImageFileVector)){ 
    
    ImageFile <- ImageFileVector[m]
    print(paste("loading folder data: ",ImageFile,sep=""))
    
    # extract fiji data
      setwd( paste(getwd(),"/",tools::file_path_sans_ext(ImageFileVector)[m],"_analysis/",sep="") )
      options(warn=-1)
      csv_FileNames <- list.files(pattern=".csv")
      csv_FileNames <- tools::file_path_sans_ext(csv_FileNames)
      csv_FileNames <- csv_FileNames[grepl( tools::file_path_sans_ext(ImageFile), csv_FileNames, fixed = TRUE)]
      csv_FileNamesMatrix <- t(sapply(strsplit(csv_FileNames, "_", fixed=T), "[", i = seq_len(max(sapply(strsplit(csv_FileNames, "_", fixed=T), length)))))
      csv_FileNamesMatrix[csv_FileNamesMatrix==""] <- NA
      csv_FileNamesMatrix <- csv_FileNamesMatrix[,colSums(is.na(csv_FileNamesMatrix))<nrow(csv_FileNamesMatrix)]
      if(ncol(csv_FileNamesMatrix)<4){ csv_FileNamesMatrix <- cbind(matrix(NA,ncol=4-ncol(csv_FileNamesMatrix),nrow=nrow(csv_FileNamesMatrix)),csv_FileNamesMatrix)}
      csv_SummaryNames <- csv_FileNames[which(csv_FileNamesMatrix[,ncol(csv_FileNamesMatrix)]%in%"summary")]
      csv_EFDataNames <- csv_FileNames[which(csv_FileNamesMatrix[,ncol(csv_FileNamesMatrix)]%in%"complexity")]
      csv_OPCDataNames <- csv_FileNames[which(csv_FileNamesMatrix[,ncol(csv_FileNamesMatrix)]%in%"kmeans")]
      csv_FD_DataNames <- csv_FileNames[grepl('FractalDimension', csv_FileNames, fixed = TRUE)]
      csv_EnamelDataNames <- csv_FileNames[which(csv_FileNamesMatrix[,ncol(csv_FileNamesMatrix)]%in%"enamel")]
      csv_DentineDataNames <- csv_FileNames[which(csv_FileNamesMatrix[,ncol(csv_FileNamesMatrix)]%in%"dentine")]
      csv_ToothDataNames <- csv_FileNames[which(csv_FileNamesMatrix[,ncol(csv_FileNamesMatrix)]%in%"remainingTooth")]
      options(warn=0)
      SpeciesName <- strsplit(csv_EFDataNames,'_')[[1]]
      SpeciesName[SpeciesName==""] <- NA
      SpeciesName <- paste(SpeciesName[!is.na(SpeciesName)][1:2],collapse="_")
    
    # Read orientation parameters (meanOrient, medianOrient, sdOrient, varOrient) 
      FijiData <- read.csv(paste(csv_EFDataNames,".csv",sep=""))
      DistanceData <- as.numeric(FijiData[,'Accumulated.distance..mm.'])
      SpatialData <- FijiData[,c('X','Y')]
      OrientationData <- as.numeric(FijiData[,'Orientation'])
      ThicknessData <- as.numeric(FijiData[,'Thickness..pix.'])
      ThicknessData2 <- as.numeric(FijiData[,'Thickness..mm.'])
      CoherData <- as.numeric(FijiData[,'Folding..1.coherency.'])
      ROIData <- as.numeric(FijiData[,'Branch'])
      
    # Read the remaining parameters
      OPC_Data <- read.csv(paste(csv_OPCDataNames,".csv",sep=""))
      FD_Data <- read.csv(paste(csv_FD_DataNames,".csv",sep=""))
      total_Data <- read.csv(paste(csv_SummaryNames,".csv",sep=""))
      tooth_Data <-  read.csv(paste(csv_ToothDataNames,".csv",sep=""))
      enamel_Data <- read.csv(paste(csv_EnamelDataNames,".csv",sep=""))
      dentine_Data <- read.csv(paste(csv_DentineDataNames,".csv",sep=""))
    
  # Read Orientation
      # Transform orientation data from radians to degrees (-90 to +90).
        OrientationData2 <- scales::rescale(OrientationData, to = c(-90,90), from = (range(c(-1.6,1.6),na.rm=T, finite=F)))
      # Transform data into circular
        OrientationData.circular <- circular(as.numeric(OrientationData2), units = "degrees", zero=pi, rotation="counter", modulo = "pi")
      # Additional circular parameters
        meanOrient <- mean.circular(OrientationData.circular)[[1]] # 'mean angle' or 'preferred direction'
        medianOrient <- median.circular(OrientationData.circular)[[1]]
        sdOrient <- sd.circular(OrientationData.circular) #See equation (2.33) at pag. 36 in Fisher (1993) for its definition. 
        varOrient <- var.circular(OrientationData.circular)[[1]]
        
  # Read enamel folding
      # Assemble EF_Data
        EF_Data <- cbind('distance'=DistanceData, 'enamel folding'=CoherData, 'ROI'=ROIData)
      # Additional thickness parameters
        totalEF <- sum(as.numeric(EF_Data[,'enamel folding']))
        maxEF <- max(as.numeric(EF_Data[,'enamel folding']))
        meanEF <- mean(as.numeric(EF_Data[,'enamel folding']))
        sdEF <- sd(as.numeric(EF_Data[,'enamel folding']))
        varEF <- var(as.numeric(EF_Data[,'enamel folding']))
        accumEF <- (cumsum(as.numeric(EF_Data[,'enamel folding']))-as.numeric(EF_Data[,'enamel folding'])[1])[length(EF_Data[,'enamel folding'])]
    
  # Read thickness
      # Additional thickness parameters
        totalET <- sum(as.numeric(ThicknessData2))
        meanET <- mean(as.numeric(ThicknessData2))
        sdET <- sd(as.numeric(ThicknessData2))
        varET <- var(as.numeric(ThicknessData2))
    
  # Read the number of patches and fractality
      OPC2D <- total_Data['X2D.OPC']
      FD <- FD_Data['D'] 
      OPC_list[[m]] <- OPC_Data
      names(OPC_list)[[m]] <- ImageFile
      
  # Read area and perimeter
      At <- tooth_Data[3]+enamel_Data[3]+dentine_Data[3]
        Abasetooth <- At
      Aen <- enamel_Data[3]
      Adentine <- dentine_Data[3]
        Aocc <- Aen + Adentine
        AocclusRel <- Aocc/At
      PerimTotal <- DistanceData[length(DistanceData)] # in mm
    
  # Additional measurements
      totalEF_At <- totalEF/At
      meanEF_At <- meanEF/At
      totalEF_Aocc <- totalEF/Aocc
      meanEF_Aocc <- meanEF/Aocc
      relEF_Perim <- totalEF/PerimTotal
      EF_Aocc <- totalEF/Aocc
      relET <- totalET/At
      mean_anis_x_thick <- mean(CoherData*ThicknessData)
      meanET_Perim <- meanET/PerimTotal # media de LocalThickness / Perimetro. LA ORIGINAL!!!!
      meanET_At <- meanET/At # media de LocalThickness /  area total
      meanET_Aocc <- meanET/Aocc # media de LocalThickness /  area occlusal
      meanET2A <- Aen/PerimTotal # area de esmalte dividido por perimetro esmalte
      meanET2B <- Aen/Aocc # area de esmalte dividido por el area total
      RELenamel <- Aen/Aocc
      OPC2D_At <- OPC2D/At
      OPC2D_Aocc <- OPC2D/Aocc
      
 # Include all the generated parameters in a list
      DataListSpecimen <- list(OrientationData, CoherData, ThicknessData)
      names(DataListSpecimen) <- c(paste(SpeciesName,"_Orientation",sep=""),paste(SpeciesName,"_Coherency",sep=""),paste(SpeciesName,"_Thickness",sep=""))
      DataList[[m]] <- DataListSpecimen
      AreaPerimData <- c(At, Aocc, Aen, RELenamel, Adentine, Abasetooth, Aocc, AocclusRel, PerimTotal) 
      names(AreaPerimData) <- c('At (mm2)', 'Aocc (mm2)', 'Aen (mm2)', 'AenRel (Aen/Aocc)', 'Adentine (mm2)', 'Abasetooth (mm2)', 'Aocc (mm2)', 'AocclusRel (Aocc/At)', 'PerimTotal (mm)')
      AnalysisData <- c(At, Aocc, meanOrient, medianOrient, sdOrient, varOrient, totalEF, meanEF, sdEF, varEF, totalEF_At, meanEF_At, totalEF_Aocc, meanEF_Aocc, EF_Aocc, accumEF, 
                        mean_anis_x_thick, totalET, meanET, meanET_At, meanET_Perim, meanET_Aocc, sdET, varET, relET, FD, OPC2D, OPC2D_At, OPC2D_Aocc)
      names(AnalysisData) <- c("At", "Aocc", "Mean orient(\BA)", "Median orient(\BA)", "SD orient", "variance orient", "totalEF", "meanEF", "sdEF", "varianceEF", "totalEF/At","meanEF/At","totalEF/Aocc","meanEF/Aocc", "EF/Aocc", "EFgain", 
                                "EF*ET", "totalET", "meanET", "meanET/At", "meanET/Perim", "meanET/Aocc", "sdET", "varianceET", "relativeET", "FD", "2D OPC", "2D OPC/At", "2D OPC/Aocc")
      AnalysisData <- unlist(AnalysisData)
    # Export tooth's parameters as csv
      write.csv(c(AreaPerimData, AnalysisData), file=paste(tools::file_path_sans_ext(ImageFile),".csv",sep=""))
    # Add to DataMatrix
      DataMatrix[,m] <- unlist(AnalysisData)
      
    # Go back to the parent directory
      setwd(wd_backup) 
      
    # Clean environment for the next loop
    }
  row.names(DataMatrix) <- names(AnalysisData)
  WD <- unlist(strsplit(working_folder, "/"))[[length(unlist(strsplit(working_folder, "/")))]]
  
# Test 1: RHINOCEROTIDAE PHYLOGENETIC MODELLING
  Tree_file <- "C:/Users/Elasmotheriini+Rhinocerotina.nex" # ADD HERE THE PATH TO THE TREE FILE E.G.: C:/Users/Elasmotheriini+Rhinocerotina.nex
  if(WD == "Rhinocerotidae dataset"){
  
  # log transform some of the variables
    LogTransform <- TRUE
    rownames(DataMatrix)
    LogTransformedVar <- c("meanEF", "2D OPC", "meanET", "At", "meanEF/At", "2D OPC/At", "meanET/At")
    if(LogTransform){
      for(k in 1:length(LogTransformedVar)){
        DataMatrix <- rbind(DataMatrix,log10(DataMatrix[LogTransformedVar[k],]))
        row.names(DataMatrix)[nrow(DataMatrix)] <- paste("log10_",LogTransformedVar[k],sep="")
       } 
    }
  
  # Read Rhino variables
    Rhino_data <- rbind(
      c('Bugtirhinus_praecursor', 'Elasmotheriini', 1, 18),
      c('Caementodon_fangxiense', 'Elasmotheriini', 2, 16),
      c('Ceratotherium_efficax', 'Rhinocerotini', 2, 4.5),
      c('Ceratotherium_mauritanicum', 'Rhinocerotini', 2, 7.5),
      c('Ceratotherium_simum', 'Rhinocerotini', 2, 5.3),
      c('Chilotheridium_pattersoni', 'Elasmotheriini', 2, 17.7),
      c('Coelodonta_antiquitatis', 'Rhinocerotini', 3, 0.46),
      c('Coelodonta_thibetana', 'Rhinocerotini', 3, 3.7),
      c('Coelodonta_tologoijensis', 'Rhinocerotini', 3, 0.7),
      c('Dicerorhinus_sumatrensis', 'Rhinocerotini', 1, 2.6),
      c('Diceros_bicornis', 'Rhinocerotini', 1, 2.5),
      c('Diceros_douariensis', 'Rhinocerotini', 1, 7),
      c('Diceros_gansuensis', 'Rhinocerotini', 1, 9),
      c('Diceros_praecox', 'Rhinocerotini', 1, 4.3),
      c('Pliorhinus_megarhinus', 'Rhinocerotini', 2, 7.5),
      c('Dihoplus_pikermiensis', 'Rhinocerotini', 2, 7.5),
      c('Dihoplus_schleiermacheri', 'Rhinocerotini', 2, 11),
      c('Elasmotherium_caucasicum', 'Elasmotheriini', 5, 1.5),
      c('Elasmotherium_chaprovicum', 'Elasmotheriini', 5, 2.6),
      c('Elasmotherium_peii', 'Elasmotheriini', 5, 2.2),
      c('Elasmotherium_primigenium', 'Elasmotheriini', 5, 9),
      c('Elasmotherium_sibiricum', 'Elasmotheriini', 5, 0.8),
      c('Eoazara_xerrii', 'Elasmotheriini', 2, 11),
      c('Hispanotherium_beonense', 'Elasmotheriini', 2, 18),
      c('Hispanotherium_corcolense', 'Elasmotheriini', 2, 18),
      c('Hispanotherium_matritense', 'Elasmotheriini', 2, 16),
      c('Hispanotherium_tungurense', 'Elasmotheriini', 2, 15),
      c('Iranotherium_morgani', 'Elasmotheriini', 3, 9.6),
      c('Lartetotherium_sansaniense', 'Rhinocerotini', 1, 17.7),
      c('Menoceras_arikarense', 'Elasmotheriini', 1, 20),
      c('Miodiceros_neumayri', 'Rhinocerotini', 3, 1.8),
      c('Ningxiatherium_euryrhinus', 'Elasmotheriini', 5, 11),
      c('Ningxiatherium_longirhinus', 'Elasmotheriini', 3, 11),
      c('Parelasmotherium_linxiaense', 'Elasmotheriini', 4, 11),
      c('Parelasmotherium_simplum', 'Elasmotheriini', 4, 11),
      c('Procoelodonta_borissiaki', 'Elasmotheriini', 2, 11),
      c('Procoelodonta_mongoliense', 'Elasmotheriini', 2, 15),
      c('Rhinoceros_platyrhinus', 'Rhinocerotini', 2, 2.6),
      c('Rhinoceros_sondaicus', 'Rhinocerotini', 1, 3),
      c('Rhinoceros_unicornis', 'Rhinocerotini', 1, 2.6),
      c('Stephanorhinus_etruscus', 'Rhinocerotini', 1, 3.3),
      c('Stephanorhinus_hemitoechus', 'Rhinocerotini', 1, 0.5),
      c('Stephanorhinus_hundsheimensis', 'Rhinocerotini', 1, 1.3),
      c('Stephanorhinus_jeanvireti', 'Rhinocerotini', 1, 3.6),
      c('Stephanorhinus_kirchbergensis', 'Rhinocerotini', 1, 0.5),
      c('Turkanatherium_acutirostratum', 'Elasmotheriini', 2, 17.7),
      c('Victoriaceros_kenyensis', 'Elasmotheriini', 2, 15)
    )
    colnames(Rhino_data) <- c('image file name', 'group', 'Hypsodonty', 'FAD')
    DataMatrix <- rbind(DataMatrix, Hypsodonty=as.numeric(Rhino_data[,'Hypsodonty']))
    write.csv(t(DataMatrix), file="AllTeeth.csv")
    
  # Read the tree
    # transform the Dataframe
      DataMatrix2 <- as.data.frame(DataMatrix)
      DataMatrixPhyloNames <- apply(matrix(unlist(strsplit(colnames(DataMatrix2), "_")), ncol=3,byrow=T)[-which(matrix(unlist(strsplit(colnames(DataMatrix2), "_")), ncol=3,byrow=T) == ""), ][,1:2],1,paste,collapse="_")
      colnames(DataMatrix2) <- DataMatrixPhyloNames
    # read ages and tree
      Tree <- read.nexus(Tree_file)
      TreeTipLabel <- Tree$tip.label
      Tree <- drop.tip(Tree, Tree$tip.label[!Tree$tip.label %in% colnames(DataMatrix2)], trim.internal = TRUE)
      ages <- data.frame(Rhino_data[,c('FAD', 'FAD')]); colnames(ages) <- c('X', 'FAD')
      ages$FAD <- as.numeric(as.character(ages$FAD))
      ages$X <- as.numeric(as.character(ages$X))
      row.names(ages) <- Rhino_data[,'image file name']; 
      Tree <- timePaleoPhy(ladderize(Tree),ages,type="equal",vartime=2,ntrees=1,add.term=T,randres=TRUE, inc.term.adj = FALSE, dateTreatment = "randObs", plot=FALSE)

  # OLS and PGLS
    {
  # Variables to test
      SelectedVars_Matrix <- rbind(
      c("log10_meanEF","At"),
      c("log10_2D OPC","At"),
      c("log10_meanET","At"),
      c("log10_meanEF/At","At"),
      c("log10_2D OPC/At","At"),
      c("log10_meanET/At","At"),
      c("log10_meanEF","log10_meanET"),
      c("log10_2D OPC","log10_meanET"),
      c("log10_2D OPC","log10_meanEF"),
      c("log10_meanEF/At","log10_meanET/At"),
      c("log10_2D OPC/At","log10_meanET/At"),
      c("log10_2D OPC/At","log10_meanEF/At"),
      c("log10_meanEF","Hypsodonty"),
      c("log10_2D OPC","Hypsodonty"),
      c("log10_meanET","Hypsodonty"),
      c("log10_meanEF/At","Hypsodonty"),
      c("log10_2D OPC/At","Hypsodonty"),
      c("log10_meanET/At","Hypsodonty")
      )
      
  # Regressions
        corr_data <- as.data.frame(t(DataMatrix2))
        corr_data$Species <- rownames(corr_data)
              
        # Generalized least squares (GLS) regression
          Corr_table <- matrix(NA, nrow=nrow(SelectedVars_Matrix),ncol=5)
          row.names(Corr_table) <- apply(SelectedVars_Matrix,1,paste,collapse=" ~ ")
          colnames(Corr_table) <- c("intercept","slope", "p-value", "RSE", "t-value")
          GLS_models <- list()
          for(m in 1:nrow(SelectedVars_Matrix)){
              SelectedVars <-  SelectedVars_Matrix[m,]
              corr_data <- as.data.frame(t(DataMatrix2))
              gls_data <- corr_data[match(Tree$tip.label,rownames(corr_data)),]
              gls_data <- gls_data[,SelectedVars[c(1,2)]]; colnames(gls_data) <- c('x','y')
              gls_data <- cbind(gls_data, Species=rownames(gls_data))
              gls_model <- nlme::gls(x ~ y, data=gls_data) 
              gls_summary <- summary(gls_model)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "t-value"] <- round(gls_summary$tTable[1,'t-value'],3)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "intercept"] <- round(gls_summary$tTable[1,'Value'],3)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "slope"] <- round(gls_summary$tTable[2,'Value'],3)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "p-value"] <- round(gls_summary$tTable[2,'p-value'],4)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "RSE"] <-  round(gls_summary$sigma,3)
              GLS_models[[m]] <- gls_model
            }
          Corr_table_GLS <- cbind(rep("GLS",times=nrow(Corr_table)),Corr_table)
          names(GLS_models) <- paste("GLS", apply(SelectedVars_Matrix,1,paste,collapse=" ~ "),sep=" ")
            
        # Ordinary least squares (OLS) regression
          Corr_table <- matrix(NA, nrow=nrow(SelectedVars_Matrix),ncol=5)
          row.names(Corr_table) <- apply(SelectedVars_Matrix,1,paste,collapse=" ~ ")
          colnames(Corr_table) <- c("intercept","slope", "p-value", "RSE", "t-value")
          OLS_models <- list()
          for(m in 1:nrow(SelectedVars_Matrix)){
              SelectedVars <-  SelectedVars_Matrix[m,]
              corr_data <- as.data.frame(t(DataMatrix2))
              ols_data <- corr_data[match(Tree$tip.label,rownames(corr_data)),]
              ols_data <- ols_data[,SelectedVars[c(1,2)]]; colnames(ols_data) <- c('x','y')
              ols_data <- cbind(ols_data, Species=rownames(ols_data))
              ols_model <- lm(formula = x ~ y, data=ols_data) 
              ols_summary <- summary(ols_model)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "t-value"] <- round(ols_summary$coefficients[2,'t value'],3)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "intercept"] <- round(ols_summary$coefficients[1,'Estimate'],3)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "slope"] <- round(ols_summary$coefficients[2,'Estimate'],3)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "p-value"] <- round(ols_summary$coefficients[2,'Pr(>|t|)'],4)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "RSE"] <-  round(ols_summary$sigma,3)
              OLS_models[[m]] <- ols_model
            }
          Corr_table_OLS <- cbind(rep("OLS",times=nrow(Corr_table)),Corr_table)
          names(OLS_models) <- paste("OLS", apply(SelectedVars_Matrix,1,paste,collapse=" ~ "),sep=" ")
            
        # Phylogenetic Generalized least squares (PGLS) regression
          Corr_table <- matrix(NA, nrow=nrow(SelectedVars_Matrix),ncol=5)
          row.names(Corr_table) <- apply(SelectedVars_Matrix,1,paste,collapse=" ~ ")
          colnames(Corr_table) <- c("intercept","slope", "p-value", "RSE", "t-value")
          PGLS_models <- list()
          for(m in 1:nrow(SelectedVars_Matrix)){
              SelectedVars <-  SelectedVars_Matrix[m,]
              corr_data <- as.data.frame(t(DataMatrix2))
              pgls_data <- corr_data[match(Tree$tip.label,rownames(corr_data)),]
              pgls_data <- pgls_data[,SelectedVars[c(1,2)]]; colnames(pgls_data) <- c('x','y')
              pgls_data <- cbind(pgls_data, Species=rownames(pgls_data))
              Species_pgls_data <- rownames(pgls_data)
              corBM <- corBrownian(phy=Tree, form= ~Species_pgls_data)
              pgls_model <- nlme::gls(x ~ y, data=pgls_data, correlation=corBM) 
              pgls_summary <- summary(pgls_model)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "t-value"] <- round(pgls_summary$tTable[1,'t-value'],3)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "intercept"] <- round(pgls_summary$tTable[1,'Value'],3)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "slope"] <- round(pgls_summary$tTable[2,'Value'],3)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "p-value"] <- round(pgls_summary$tTable[2,'p-value'],4)
              Corr_table[paste(SelectedVars,collapse=" ~ "), "RSE"] <-  round(pgls_summary$sigma,3)
              PGLS_models[[m]] <- pgls_model
            }
          Corr_table_PGLS <- cbind(model=rep("PGLS",times=nrow(Corr_table)),Corr_table)
          names(PGLS_models) <- paste("PGLS", apply(SelectedVars_Matrix,1,paste,collapse=" ~ "),sep=" ")
            
        # Merge OLS and PGLS
          if(exists("Corr_table_PGLS") && exists("Corr_table_OLS")){ 
              Corr_table <- rbind(Corr_table_OLS, Corr_table_PGLS)
              row.names(Corr_table) <- gsub("[.]", " ", row.names(Corr_table)); row.names(Corr_table) <- gsub("  ", " ", row.names(Corr_table)); row.names(Corr_table) <- gsub("X", "", row.names(Corr_table)); 
              Corr_table[,'p-value'][Corr_table[,'p-value']<0.1] <- paste(Corr_table[,'p-value'],"*",sep="")[Corr_table[,'p-value']<0.1]
              Corr_table[,'p-value'][Corr_table[,'p-value']<0.05] <-  paste(Corr_table[,'p-value'],"*",sep="")[Corr_table[,'p-value']<0.05] 
              Corr_table[,'p-value'][Corr_table[,'p-value']<0.001] <- paste(Corr_table[,'p-value'],"*",sep="")[Corr_table[,'p-value']<0.001] 
              Corr_table[,'p-value'][Corr_table[,'p-value']<0.001] <- '<0.001***'
      }
          Corr_table <- Corr_table[order(Corr_table[,1]),]
          rownamesCorr_table <- row.names(Corr_table)
          rownamesCorr_table <- gsub("_"," ",rownamesCorr_table); rownamesCorr_table <- gsub("X","",rownamesCorr_table)
          rownamesCorr_table <- gsub("Mean enamel folding","EF",rownamesCorr_table);
          rownamesCorr_table <- gsub("At mm2","At",rownamesCorr_table)
          row.names(Corr_table) <- rownamesCorr_table
          colnames(Corr_table)[1] <- 'method'
          write.csv(Corr_table,"Corr_table.csv")
          print(Corr_table)
    }
    
  # Phylogenetic modelling
    PhenogramPlot <- TRUE
    PGLS <- FALSE
    OUWIE <- TRUE
    PhenoSpreadLabels <- FALSE
    VarsToPlot <- c("log10_meanEF/At","log10_2D OPC/At")
    GrousToPlot <- c("Elasmotheriini", "Rhinocerotini")
    cols_pheno <- c("orangered","chartreuse3","gold2")
    FAD_names <- c("FAD","LAD") 
    VarIndep <- "log10_meanET/At"
    s_aic_threshold <- 6
    { 
    # Prepare the data matrix
       DataMatrixPhylo <- t(DataMatrix)
       row.names(DataMatrixPhylo) <- DataMatrixPhyloNames

      # extract the list of Elasmotheriina and Rhinocerotina
       tribe <- rep("Elasmotheriini",times=nrow(DataMatrixPhylo))
       tribe[1:23] <- "Rhinocerotini"
       names(tribe) <- TreeTipLabel
       names(cols_pheno)<- c(unique(tribe),"NA")
      
     # Select taxa
       SelectedTaxa <- names(tribe[tribe%in%GrousToPlot])
       DataMatrixTree <-  DataMatrixPhylo[row.names(DataMatrixPhylo) %in% SelectedTaxa,]
       CSVnames <- row.names(DataMatrixTree)
       Tree <- drop.tip(Tree, setdiff(Tree$tip.label, SelectedTaxa));

    # Colors by group
      Tree_tribe <- Tree
      if((length(unique(tribe[tribe %in% GrousToPlot])))>1){
        Tree_tribe <- paintSubTree(Tree_tribe,node=48,state="NA")
        Tree_tribe <- paintSubTree(Tree_tribe,node=49,state="Elasmotheriini")
        Tree_tribe <- paintSubTree(Tree_tribe,node=71,state="Rhinocerotini")
      }else{Tree_tribe <- Tree}
      
    # Plot phenograms
      PhenogramPlot <- TRUE
      if(PhenogramPlot){
     # trait palette
      require(wesanderson)
      TraitPalette <- wes_palette("Zissou1", 100, type = "continuous")
     # plot
      VarsToPlot <- c(VarsToPlot, VarIndep)
       for(k in 1:length(VarsToPlot)){
         windows(height = 7, width = 12)
         par(mfrow=c(1,2))
          # fenograma
            PhenoYlabel <- VarsToPlot[k]
            PhenoData <- DataMatrixTree[,PhenoYlabel]
            phenogram(Tree_tribe, PhenoData, ylab=PhenoYlabel, spread.labels=PhenoSpreadLabels, fsize=0.8, ftype="i", 
                      colors = cols_pheno, main=VarsToPlot[k], xlab="time (Ma)", axes=list("time","trait"))
          # plotear caracteres
            objTree <- contMap(Tree, PhenoData, plot=FALSE)
            objTree_palette <- setMap(objTree, TraitPalette)
          # cambiar colores
            n_cols <- n_distinct(PhenoData)
            plot(objTree_palette, fsize=c(0.8,1), outline=FALSE, lwd=c(3,7), mtext=VarsToPlot[k])
            mylims <- par("usr")
            # plotear hipsodoncia en caso de ser complejidad o 2D OPC
            if(any(VarsToPlot[k]==c("log10 Mean.enamel folding","log10 X2D.OPC"))){
              Ordered_Hypsodonty <- Hypsodonty[objTree_palette$tree$tip.label %in% names(Hypsodonty)]
              Ordered_Hypsodonty <- rev(Ordered_Hypsodonty[match(objTree_palette$tree$tip.label,names(Ordered_Hypsodonty))])
              TraitPalette2 <- wes_palette(name="Zissou1", max(Ordered_Hypsodonty), type = "continuous")
              Ordered_Hypsodonty_cols <- rep(NA,times=length(Ordered_Hypsodonty))
              for(j in 1:length(Ordered_Hypsodonty_cols)){Ordered_Hypsodonty_cols[j] <- TraitPalette2[Ordered_Hypsodonty[j]]}
              par(oma=c(0,65,0,6))
              points(rep(1,times=length(objTree_palette$tree$tip.label)), c(1:length(objTree_palette$tree$tip.label)), 
                     pch=15, cex=1.5, col=Ordered_Hypsodonty_cols)
            }
            pdf_name <- paste("Phenogram_",gsub("[.]","_",VarsToPlot[k]),".pdf",sep="")
            pdf_name <- paste("Phenogram_",gsub("/","_div_",VarsToPlot[k]),".pdf",sep="")
            dev.copy(pdf, file = pdf_name, width=12, height=7); dev.off()
        }
       par(mfrow=c(1,1))
    }
     
  # PHYLOGENETIC MODELLING: SURFACE
     phylomorphospace_labels <- "horizontal" # "horizontal or "off"
     VarDep <- VarsToPlot[VarsToPlot!=VarIndep]
     cols2_list <- list()
     t_list <- list()
     OptimaCoords_list <- list()
     
     # plot phylomorphospace
     for(h in VarDep){
     # select variables
       SelectedVars <- c(VarIndep,h)
     
     # divide by clades?
       splitClades <- FALSE
       if(splitClades){n <- 1 }else{n=c(1,2)}
       GroupSurf <- unique(tribe)[n]
       SelectedTaxa <- names(tribe[tribe%in%GroupSurf])
       
     # name nodes of the trees
       surface_tree <- Tree
       surface_tree$node.label <- c(1:(length(surface_tree$tip.label)-1))
       
     # data
       surface_data <- as.data.frame(DataMatrixTree)
       surface_data <- surface_data[,SelectedVars,drop=FALSE]
       colnames(surface_data) <- gsub("/",".",colnames(surface_data)); colnames(surface_data) <- gsub("[.]","_",colnames(surface_data)); colnames(surface_data) <- gsub(" ","_",colnames(surface_data))
         
     # run surface
        s <- runSurface(surface_tree, surface_data, error_skip=F, plotaic=F, verbose=T, aic_threshold = s_aic_threshold)
     
     # results of the best surface model
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
           # phylomorphospace
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
           lim_plot_y <- c(min(surface_data[,2]),max(surface_data[,2]))
           par(mar=c(6,5,3,1))
           plot(NULL, xlim=lim_plot_x, ylim=lim_plot_y, xlab=colnames(surface_data)[1], ylab=colnames(surface_data)[2], bty="n")
           points(OptimaCoords, pch=16, cex=11, col=cols2)
           phylomorphospace(t, as.matrix(surface_data), node.size=c(0,1.5), pch=16, lwd=2, cex=0.6, colors=cols2,
                            xlim=lim_plot_x, ylim=lim_plot_y, label=phylomorphospace_labels, node.by.map=TRUE, add=T)
           mtext(paste(gsub("_"," ",colnames(surface_data)),collapse=" vs "), side = 3, line = - 2, outer = TRUE)
           dev.copy(pdf, file = paste(paste(colnames(surface_data),collapse=" vs "),".pdf",sep=""), width=12, height=7); dev.off()
           
           # same plot, but colored by groups
           windows(height = 7, width = 12)
           par(mfrow=c(1,2), mar=c(0,0,0,0))
           surfaceTreePlot(surface_tree, hansenfit=s$bwd[[kk]], labelshifts = T, cols=cols, edge.width=3, cex=.7, convcol=T)
           par(mar=c(6,5,3,1))
           plot(NULL, xlim=lim_plot_x, ylim=lim_plot_y, xlab=colnames(surface_data)[1], ylab=colnames(surface_data)[2], bty="n")
           points(OptimaCoords, pch=16, cex=11, col=cols2)
           phylomorphospace(Tree_tribe, as.matrix(surface_data), node.size=c(0,1.5), pch=16, lwd=2, cex=0.6, colors=cols_pheno,
                            xlim=lim_plot_x, ylim=lim_plot_y, label=phylomorphospace_labels, node.by.map=TRUE, add=T)
           mtext(paste(gsub("_"," ",colnames(surface_data)),collapse=" vs "), side = 3, line = - 2, outer = TRUE)
           dev.copy(pdf, file = paste(paste(colnames(surface_data),collapse=" vs "),"_groups.pdf",sep=""), width=12, height=7); dev.off()
           par(mfrow=c(1,1))
         }
      }  
  }
    
    # Export data
      save(DataMatrixTree, VarsToPlot, Tree, file="caper_data2.RData")
    
    # Plot phenogram with the evolutionary regimes
      DataMatrixTree2 <- DataMatrixTree
      colnames(DataMatrixTree2) <- gsub("[.]","_", colnames(DataMatrixTree2))
      colnames(DataMatrixTree2) <- gsub("/","_", colnames(DataMatrixTree2))
      windows(height = 7, width = 12)
      par(mfrow=c(1,3))
      for(h in 1:length(t_list)){
        SelectedVar <- names(t_list)[h]; SelectedVar <- gsub("2D_OPC", "2D OPC", SelectedVar)
        PhenoData <- DataMatrixTree2[,SelectedVar]
        phenogram(t_list[[h]], PhenoData, ylab=PhenoYlabel, spread.labels=PhenoSpreadLabels, fsize=0.8, ftype="i", 
            colors = cols2_list[[h]], main=SelectedVar, xlab="time (Ma)", axes=list("time","trait"))
    }
      # extract the tree of the independent variable
      phenogram(t_list[[1]], DataMatrixTree2[,gsub("[.]","_",'log10_meanET_At')], ylab=PhenoYlabel, spread.labels=PhenoSpreadLabels, fsize=0.8, ftype="i", 
        colors = cols2_list[[1]], main='log10_Mean_thickness_At__mm2_', xlab="time (Ma)", axes=list("time","trait"))
      dev.copy(pdf, file = paste("phenograms_regimes.pdf",sep=""), width=12, height=7); dev.off()
      par(mfrow=c(1,1))
  
  }
    
# Test 2: resolution test
  if(WD == "Resolution dataset"){
   DataMatrix2plot <- DataMatrix
   ColNames <- matrix(unlist(strsplit(colnames(DataMatrix2plot),"_")),ncol=6,byrow=T)
   ColNames <- cbind(paste(ColNames[,3], ColNames[,4]), ColNames[,6])
   ColNames[,2] <- gsub("px","", ColNames[,2])
   var2plotVector <- c('2D OPC','total EF')
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
  if(WD == "Wear dataset 1" || WD == "Wear dataset 2"){
    PlotAll <- FALSE
    colnamesDataMatrix <- colnames(DataMatrix)
    Hsteps <- as.numeric(matrix(unlist(strsplit(colnamesDataMatrix,"_")),ncol=4,byrow=T)[,3])
    Data2plot <- DataMatrix[,order(Hsteps)]
    Hsteps <- rev(sort(Hsteps))
    VarsToPlot <- c(1:nrow(Data2plot))
    print(row.names(Data2plot))
    VarsToPlot <- c("mean ET", "total EF", "2D OPC")
    par(mfrow=c(1,3), mai = c(0.8, 0.8, 0.5, 0.3))
    for(n in 1:length(VarsToPlot)){
      RowToPlot <- VarsToPlot[n]
      maxInterval <- Data2plot[RowToPlot,]+5
      minInterval <- Data2plot[RowToPlot,]-5
      plot(Hsteps, Data2plot[RowToPlot,],type="o", ylab=RowToPlot, ylim=c(0,max(Data2plot[RowToPlot,])), 
           xlim=c(max(Hsteps),min(Hsteps)), xlab="tooth's crown height")
    }
    par(mfrow=c(1,1))
    dev.copy(pdf, file = paste("wear_test_2.pdf",sep=""), width=12, height=4); dev.off()
  }
  
# Test 5: rotation test
  if(WD == "Rotation dataset"){
   DataMatrix2plot <- DataMatrix
   ColNames <- matrix(unlist(strsplit(colnames(DataMatrix2plot),"_")),ncol=6,byrow=T)
   ColNames[ColNames==""] <- NA; ColNames <- ColNames[,colSums(is.na(ColNames))<nrow(ColNames)]
   ColNames <- cbind(paste(ColNames[,1], ColNames[,2]), ColNames[,4])
   ColNames[,2] <- gsub("px","", ColNames[,2])
   var2plotVector <- c('2D OPC','total EF')
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
  if(WD == "Suidae dataset"){
  # extract matrix
     DataMatrix2plot <- DataMatrix
     ColNames <- matrix(unlist(strsplit(colnames(DataMatrix2plot),"_")),ncol=6,byrow=T)
     ColNames[ColNames==""] <- NA; ColNames <- ColNames[,colSums(is.na(ColNames))<nrow(ColNames)]
     ColNames <- cbind(paste(ColNames[,1], ColNames[,2]),ColNames[,3], ColNames[,4])
     PigNames <- c("Hylochoerus meinertzhageni","Phacochoerus africanus","Sus scrofa","Potamochoerus porcus")
     TeethOrder <- c("p2","p3","p4","m1","m2","m3","combined"); TeethOrder <- cbind(TeethOrder,1:length(TeethOrder))
     var2plotVector <- c("total ET","2D OPC")
     Palette_Suidae <- c("#778c93","#0d5888","#831e11","#776813","#ceab0d")
   # plot
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
       # plot
       xlimPlot <- as.numeric(c(1,length(unique(TeethOrder[,1]))))
       ylimPlot <- c(min(as.numeric(Matrix2plot[,3])),max(as.numeric(Matrix2plot[,3])))
       plot(0,0, xlim = xlimPlot, ylim = ylimPlot, ylab=var2plot, xlab="tooth position", xaxt = "n", type = "l"); axis(1, at=1:nrow(Matrix2plot2), labels=row.names(Matrix2plot2))
       cl <- Palette_Suidae
       for (i in 1:ncol(Matrix2plot2)){
         par(new=TRUE); plot(c(1:nrow(Matrix2plot2))[!is.na(Matrix2plot2[,i])],Matrix2plot2[,i][!is.na(Matrix2plot2[,i])], axes=F, ylab=NA, xlab=NA, xlim = xlimPlot, ylim = ylimPlot, col = cl[i], type = 'o', pch=16)
         legend("topleft", legend=colnames(Matrix2plot2), col=cl, pch=16, bty = "n") # optional legend
        }
    }
     par(mfrow=c(1,1))
 }
 
# Test 7: Lagomorpha
  if(WD == "Lagomorpha dataset"){
    require(ggpubr)
    # log transform some of the variables
      LogTransform <- TRUE
      LogTransformedVar <- c("mean EF","2D OPC")
      if(LogTransform){
        for(k in 1:length(LogTransformedVar)){
          DataMatrix <- rbind(DataMatrix,log10(DataMatrix[LogTransformedVar[k],]))
          row.names(DataMatrix)[nrow(DataMatrix)] <- paste("log10_",LogTransformedVar[k],sep="")
        } 
      }
    
    # plot
      print(row.names(DataMatrix))
      var2plotVector <- c('log10_2D OPC','log10_mean EF')
      DataMatrix2plot <- as.data.frame(t(DataMatrix[var2plotVector,]))
      ColNames <- matrix(unlist(strsplit(row.names(DataMatrix2plot),"_")),ncol=5,byrow=T)
      ColNames[ColNames==""] <- NA; ColNames <- ColNames[,colSums(is.na(ColNames))<nrow(ColNames)]
      ColNames <- apply(ColNames[,1:3],1,paste,collapse="_")
      DataMatrix2plot <- cbind(species=ColNames, DataMatrix2plot)
      Hypsodonty_names <- c("Prolagus_tobieni", "Prolagus_vasconiensis", "Prolagus_oeningensis", "Prolagus_michauxi", "Prolagus_affdepereti", "Prolagus_depereti", 
                            "Prolagus_apricenicus", "Prolagus_imperialis", "Prolagus_figaro", "Prolagus_cfsardus","Prolagus_sardus", "Prolagus_crusafonti","Prolagus_major",
                            "Prolagus_calpensis")
      Hypsodonty_value <- c(4.31,4.57,4.69,4.53,NA,NA,4.73,4.83,4.53,6.33,6.33, NA,NA,NA)
      Insular_value <- c("C","C","C","C","I","I","I","I","I","I","I","C","C","C")
      Hypsodonty <- cbind(Hypsodonty_names, Hypsodonty_value, Insular_value)
      attachedMatrix <- matrix(NA,ncol=2,nrow=nrow(DataMatrix2plot))
      for(j in 1:nrow(attachedMatrix)){ attachedMatrix[j,] <- Hypsodonty[match(DataMatrix2plot[j,1],Hypsodonty[,1]),2:3]}
      DataMatrix2plot <- cbind(DataMatrix2plot, attachedMatrix) 
      TraitPalette_microwear <- c("#778c93","#0d5888","#831e11","#776813","#ceab0d", "", "", "", "")
      cl <- TraitPalette_microwear[2:length(TraitPalette_microwear)]
      par(mfrow=c(1,1), mai = c(0.8, 0.8, 0.5, 0.3))
      plot(DataMatrix2plot[,2],DataMatrix2plot[,3],xlab=colnames(DataMatrix2plot)[2],ylab=colnames(DataMatrix2plot)[3],pch=16)
      text(DataMatrix2plot[,2],DataMatrix2plot[,3], labels=DataMatrix2plot[,1],cex=0.9)
      dev.copy(pdf, file = paste("Prolagus.pdf",sep=""), width=10, height=7); dev.off()
      SelectedVars1 <- "At"
      SelectedVars2 <- c("log10_mean EF", "log10_2D OPC")
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
  if(WD == "Proboscidea dataset"){
    print(row.names(DataMatrix))
    var1 <- 'mean ET'
    var2 <- c('2D OPC','mean EF')
    OPC_data <- OPC_list[[1]]
    DataMatrix2 <- as.data.frame(DataMatrix)
    DataMatrixPhyloNames <- apply( matrix(unlist(strsplit(colnames(DataMatrix2), "_")), ncol=5,byrow=T)[,3:4],1,paste,collapse="_")
    colnames(DataMatrix2) <- DataMatrixPhyloNames
    var2plotVector <- c(var1, var2)
    DataMatrix2 <- t(DataMatrix2[var2plotVector,])
    par(mfrow=c(1,1), mai = c(0.8, 0.8, 0.5, 0.3))
    for(n in 1:2){
      plot(DataMatrix2[,var1],DataMatrix2[,var2[n]], ylab=var2[n], xlab=var1, col=rainbow(nrow(DataMatrix2)), pch=16, cex=2)
      legend("bottomright", legend=row.names(DataMatrix2), col=rainbow(nrow(DataMatrix2)), pch=16, bty = "n") 
      dev.copy(pdf, file = paste(gsub('[.]',"_",var1),"_vs_",gsub('[.]',"_",var2[n]),'.pdf',sep=""), width=10, height=7); dev.off()
    }
    barplot_col <- c("#E92210","#E0AF12","#BEDC14","#2AD316","#18C89F","#1A8FC2","#0A31AE","#2F0277")
    barplot_col <- c("#2F0277", "#0A31AE", "#1A8FC2", "#18C89F", "#2AD316","#BEDC14","#E0AF12","#E92210")
    par(mfrow=c(2,2), mai = c(0.8, 0.8, 0.5, 0.3))
    for(n in 1:length(OPC_list)){
      barplot(OPC_list[[n]][,'Area..pixels.2.'], names.arg=OPC_list[[n]][,'Mean'], col=barplot_col, main=names(OPC_list)[n], space=0, border=NA)
    }
  }
  
# Test 9: 2D OPC vs 3D OPC
  if(WD == "2D vs 3D OPC dataset"){
    var2plotVector <- c('2D OPC')
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
      
    DataMatrix2 <- as.data.frame(DataMatrix)
    DataMatrixPhyloNames <- apply( matrix(unlist(strsplit(colnames(DataMatrix2), "_")), ncol=5,byrow=T)[,3:4],1,paste,collapse="_")
    colnames(DataMatrix2) <- DataMatrixPhyloNames
    DataMatrix2 <- t(DataMatrix2[var2plotVector,])
    DataMatrix2 <- DataMatrix2[match(row.names(OPC3D), row.names(DataMatrix2)),]
    DataMatrix2 <- cbind(OPC3D, DataMatrix2)
    colnames(DataMatrix2)[ncol(DataMatrix2)] <- '2D OPC (1000px)'
    TraitPalette_OPC <- c("skyblue","skyblue1","skyblue2","skyblue3","skyblue4","black")
   
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
   
  }
