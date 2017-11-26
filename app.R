library (lubridate)
library (data.table)
library (utils)
library (shiny)
library (shinyjs)
library (shinydashboard)
library (magrittr)
library (ggplot2)
library (plotly)
library (DT)
#library (ega)
library (mgcv)
library (mcr)
library (rhandsontable)
library (rmarkdown)
library (shinyAce)


#utility functions
'%nin%' <- Negate('%in%')#to find elements not in vector

replacer <- function (df, replace=NA, with=0){
  h <- is.vector(df); i <- is.matrix(df); j <- is.data.frame(df)
  y <- as.matrix(df)
  if (is.na(replace)) {
    y[is.na(y)] <- with
  } else { 
    y[y==replace] <- with
  }
  if(h) y <- as.vector(y)
  if(i) y <- as.matrix(y)
  if(j) y <- as.data.frame(y)
  return (y)
} #replace values with others (with==) for different data types

#EGA functions

#finding tan coefficient for a line going trough 2 points
coef <- function (x, y, xend, yend){
  if (xend==x) {
    stop("Vertical line - function inapplicable")
  }
  return ( (yend - y) / (xend - x))
}

#setting y axis end point with tan coeff and x endpoint known
endy <- function (startx, starty, maxx, coef){
  return ( (maxx - startx) * coef + starty)
}

#setting x axis end point with tan coeff and y endpoint known
endx <- function (startx, starty, maxy, coef){
  return ( (maxy - starty)  / coef + startx)
}

getClarkeZones <- function (referenceVals, testVals, unit="gram"){
  if (unit != "mol" & unit != "gram") {
    stop("'unit' must be either 'mol' or 'gram'.")
  }
  if (unit == "mol") {
    n <- 18 #scaling factor for mmol/l conversion from mg/dl
  } else {
    n <- 1
  }
  zones <- vector (mode="character", length=length (referenceVals))
  bias <- testVals - referenceVals
  
  # absolute relative error = abs(bias)/reference*100
  are <- abs (bias) / referenceVals * 100
  eq1 <- (7 / 5) * (referenceVals - 130 / n)
  eq2 <- referenceVals + 110 / n
  
  # zone D: ref < 70 and (test > 70 and test < 180) or
  #   ref > 240 and (test > 70 and test < 180)
  test_D <- testVals >= 70 / n & testVals < 180 / n#error corrected >=70 instead of >70
  zoneD <- (referenceVals < 70 / n & test_D) |
    (referenceVals > 240 / n & test_D)
  zones[zoneD] <- "D"
  
  # assign A after D, since part of A will overwrite D
  
  # zone C: (ref >= 130 and ref <= 180 and test < eq1) or
  #   (ref > 70 and ref > 180 and ref > eq2)
  zoneC <- (referenceVals >= 130 / n & referenceVals <= 180 / n & testVals < eq1) |
    (referenceVals > 70 / n & testVals > 180 / n & testVals > eq2)
  zones[zoneC] <- "C"
  
  #Assign A after C, since part of C will override A
  
  # zone A: are <= 20  or (ref < 58.3 and test < 70)
  zoneA <- (are <= 20) |
    (referenceVals < 70 / n & testVals < 70 / n)#error solved
  zones[zoneA] <- "A"
  
  # zone E: (ref <= 70 and test >= 180) or (ref >=180 and test <=70)
  zoneE <- (referenceVals <= 70 / n & testVals >= 180 / n) |
    (referenceVals >= 180 / n & testVals <= 70 / n)
  zones[zoneE] <- "E"
  
  # the rest are zone B
  zones <- replace (zones, zones == "", "B")
  return (zones)
}

getParkesZones <- function (ref, test, type=1, unit="gram"){
  if (type != 1 & type != 2){
    stop("'type' must be 1 or 2.")
  }
  if (unit != "mol" & unit != "gram") {
    stop("'unit' must be either 'mol' or 'gram'.")
  }
  if (unit == "mol") {
    n <- 18 #scaling factor for mmol/l conversion from mg/dl
  } else {
    n <- 1
  }
  
  #setting the graph limits with some space to accomondate all the datapoints
  maxX <- max (max (ref) + 20 / n, 550 / n)#better solution
  maxY <- max ( (test + 20 / n), maxX, 550 / n)
  
  #common block
  testdf <- as.matrix (cbind (as.numeric (ref), as.numeric (test)))
  zones <- vector (mode = "character", length = length (ref))
  zones[1:length (ref)] <- "A" #all datapoints are A by default
  
  #diffirenciation by diabites type
  if (type==1){
    #determining line coefficients for final segments
    ce <- coef (35, 155, 50, 550)
    cdu <- coef (80, 215, 125, 550)
    cdl <- coef (250, 40, 550, 150)
    ccu <- coef (70, 110, 260, 550)
    ccl <- coef (260, 130, 550, 250)
    cbu <- coef (280, 380, 430, 550)
    cbl <- coef (385, 300, 550, 450)
    
    #diabetes type 1 zones - creates polygons of a size dependent on data
    limitE1 <- matrix (data= c (0, 35 / n, endx (35 / n, 155 / n, maxY, ce), 0, 0, #x limits E upper
                                150 / n, 155 / n, maxY, maxY, 150 / n),#y limits E upper
                       ncol=2, byrow=FALSE)
    
    limitD1L <- matrix (data= c (250 / n, 250 / n, maxX, maxX, 250 / n,#x limits D lower
                                 0, 40 / n, endy (410 / n, 110 / n, maxX, cdl), 0, 0),#y limits D lower
                        ncol=2, byrow=FALSE)
    
    limitD1U <- matrix (data= c (0, 25 / n, 50 / n, 80 / n, endx (80 / n, 215 / n, maxY, cdu), 0, 0,#x limits D upper
                                 100 / n, 100 / n, 125 / n, 215 / n , maxY, maxY, 100 / n),#y limits D upper
                        ncol=2, byrow=FALSE)
    
    limitC1L <- matrix (data= c (120 / n, 120 / n, 260 / n, maxX, maxX, 120 / n, #x limits C lower
                                 0, 30 / n, 130 / n, endy (260 / n, 130 / n, maxX, ccl) , 0, 0),#y limits C lower
                        ncol=2, byrow=FALSE)
    
    limitC1U <- matrix (data= c (0, 30 / n, 50 / n, 70 / n, endx (70 / n, 110 / n, maxY, ccu), 0, 0, #x limits C upper
                                 60 / n, 60 / n, 80 / n, 110 / n , maxY, maxY, 60 / n),#y limits C upper
                        ncol=2, byrow=FALSE)
    
    limitB1L <- matrix (data= c (50 / n, 50 / n, 170 / n, 385 / n, maxX, maxX, 50 / n, #x limits B lower
                                 0, 30 / n, 145 / n, 300 / n , endy (385 / n, 300 / n, maxX, cbl), 0, 0),#y limits B lower
                        ncol=2, byrow=FALSE)
    
    
    limitB1U <- matrix (data= c (0, 30 / n, 140 / n, 280 / n, endx (280 / n, 380 / n, maxY, cbu), 0, 0, #x limits B upper
                                 50 / n, 50 / n, 170 / n, 380 / n , maxY, maxY, 50 / n),#y limits B upper
                        ncol=2, byrow=FALSE)
    
    #labelling zones using in.out function from mgcv package
    zones[which (in.out (limitB1L, testdf))] <- "B"
    zones[which (in.out (limitB1U, testdf))] <- "B"
    zones[which (in.out (limitC1L, testdf))] <- "C"
    zones[which (in.out (limitC1U, testdf))] <- "C"
    zones[which (in.out (limitD1L, testdf))] <- "D"
    zones[which (in.out (limitD1U, testdf))] <- "D"
    zones[which (in.out (limitE1, testdf))] <- "E"
    
  }else{ #type 2 diabetes
    ce <- coef (35, 200, 50, 550)
    cdu <- coef (35, 90, 125, 550)
    cdl <- coef (410, 110, 550, 160)
    ccu <- coef (30, 60, 280, 550)
    ccl <- coef (260, 130, 550, 250)
    cbu <- coef (230, 330, 440, 550)
    cbl <- coef (330, 230, 550, 450)
    #diabetes type 2 zones - creates polygons of a size dependent on data
    limitE2 <- matrix (data= c (0, 35 / n, endx (35 / n, 200 / n, maxY, ce), 0, 0, #x limits E upper
                                200 / n, 200 / n, maxY, maxY, 200 / n),#y limits E upper
                       ncol=2, byrow=FALSE)
    
    limitD2L <- matrix (data= c (250 / n, 250 / n, 410 / n, maxX, maxX, 250 / n,#x limits D lower
                                 0, 40 / n, 110 / n, endy (410 / n, 110 / n, maxX, cdl), 0, 0),#y limits D lower
                        ncol=2, byrow=FALSE)
    
    limitD2U <- matrix (data= c (0, 25 / n, 35 / n, endx (35 / n, 90 / n, maxY, cdu), 0, 0, #x limits D upper
                                 80 / n, 80 / n, 90 / n, maxY, maxY, 80 / n),#y limits D upper
                        ncol=2, byrow=FALSE)
    
    limitC2L <- matrix (data= c (90 / n, 260 / n, maxX, maxX, 90 / n, #x limits C lower
                                 0, 130 / n, endy (260 / n, 130 / n, maxX, ccl), 0, 0),#y limits C lower
                        ncol=2, byrow=FALSE)
    
    limitC2U <- matrix (data= c (0, 30 / n, endx (30 / n, 60 / n, maxY, ccu), 0, 0, #x limits C upper
                                 60 / n, 60 / n, maxY, maxY, 60 / n),#y limits C upper
                        ncol=2, byrow=FALSE)
    
    limitB2L <- matrix (data= c (50 / n, 50 / n, 90 / n, 330 / n, maxX, maxX, 50 / n, #x limits B lower
                                 0, 30 / n, 80 / n, 230 / n , endy (330 / n, 230 / n, maxX, cbl), 0, 0),#y limits B lower
                        ncol=2, byrow=FALSE)
    
    limitB2U <- matrix (data= c (0, 30 / n, 230 / n, endx (230 / n, 330 / n, maxY, cbu), 0, 0, #x limits B upper
                                 50 / n, 50 / n, 330 / n , maxY, maxY, 50 / n),#y limits B upper
                        ncol=2, byrow=FALSE)
    
    #labelling zones using in.out function from mgcv package
    zones[which (in.out (limitB2L, testdf))] <- "B"
    zones[which (in.out (limitB2U, testdf))] <- "B"
    zones[which (in.out (limitC2L, testdf))] <- "C"
    zones[which (in.out (limitC2U, testdf))] <- "C"
    zones[which (in.out (limitD2L, testdf))] <- "D"
    zones[which (in.out (limitD2U, testdf))] <- "D"
    zones[which (in.out (limitE2, testdf))] <- "E"
  }
  return (zones)
}

plotClarkeGrid <- function (referenceVals, testVals, title = "Clarke Error Grid", xlab="", ylab="",
                            linesize = 0.5, linetype = "solid", linecolor = "black",
                            linealpha = 0.6, pointsize = 2, pointalpha = 1, zones = NA, unit='gram')
{
  #unit selection and control
  if (unit != "mol" & unit != "gram") {
    stop("'unit' must be either 'mol' or 'gram'.")
  }
  #axis named depending on unit
  if (unit == "mol") {
    n <- 18 #where n is a scaling factor for unit conversion
    if (xlab==""){
      xlab="Reference Glucose Concentration (mmol/L)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mmol/L)"
    }
  } else {
    n <- 1
    if (xlab==""){
      xlab="Reference Glucose Concentration (mg/dL)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mg/dL)"
    }
  }
  
  # use default zone assignment if none is provided
  if (is.na (zones)) {
    zones <- getClarkeZones (referenceVals, testVals, unit)
  }
  tolerance <- 0.2
  
  # create a df for ggplot
  data <- data.frame (ref=referenceVals, test=testVals, zones=zones)
  
  #better solution for scaling axis automatically with some extra space
  maxX <- max (max (data$ref) + 20 / n, 550 / n)
  maxY <- max ( (data$test + 20 / n), (1 + tolerance) * maxX, 650 / n)
  
  #labels with coordinats and colors
  labels <- data.frame (x=c (240 / n, 120 / n, 350 / n, 120 / n, 163 / n,
                             35 / n, 350 / n, 35 / n, 350 / n),
                        y=c (230 / n, 200 / n, 230 / n, 300 / n, 20 / n,
                             130 / n, 130 / n, 300 / n, 35 / n),
                        label=c ("A", "B", "B", "C", "C", "D", "D", "E", "E"),
                        color=c ("blue", "blue", "blue", "blue", "blue",
                                 "red", "red", "red", "red"))
  
  #segment endpoints for borders
  border <- data.frame (x1=c (58.3 / n, 70 / n, 70 / n, 70 / n, 0, 240 / n,
                              0, 70 / n, 180 / n, 240 / n, 180 / n, 130 / n),
                        y1=c (70 / n, 56 / n, 180 / n, 83 / n, 180 / n,
                              180 / n, 70 / n, 0, 70 / n, 70 / n, 0, 0),
                        xend=c (maxX, maxX, 550 / n, 70 / n, 70 / n, maxX,
                                58.3 / n, 70 / n, maxX, 240 / n,
                                180 / n, 180 / n),
                        yend=c ((1 + tolerance) * maxX, (1 - tolerance) * maxX,
                                660 / n, maxY, 180 / n, 180 / n, 70 / n,
                                56 / n, 70 / n, 180 / n, 70 / n, 70 / n))
  
  ceg <- ggplot(data, aes(x = ref, y = test)) +
    scale_x_continuous (breaks = c (round (70 / n, digits=1),
                                    round (100 / n, digits=1),
                                    round (150 / n, digits=1),
                                    round (180 / n, digits=1),
                                    round (240 / n, digits=1),
                                    round (300 / n, digits=1),
                                    round (350 / n, digits=1),
                                    round (400 / n, digits=1),
                                    round (450 / n, digits=1),
                                    round (500 / n, digits=1),
                                    round (550 / n, digits=1),
                                    round (600 / n, digits=1),
                                    round (650 / n, digits=1),
                                    round (700 / n, digits=1),
                                    round (750 / n, digits=1),
                                    round (800 / n, digits=1),
                                    round (850 / n, digits=1),
                                    round (900 / n, digits=1),
                                    round (950 / n, digits=1),
                                    round (1000 / n, digits=1)),
                        expand = c (0, 0)) +
    scale_y_continuous (breaks = c (round (70 / n, digits=1),
                                    round (100 / n, digits=1),
                                    round (150 / n, digits=1),
                                    round (180 / n, digits=1),
                                    round (240 / n, digits=1),
                                    round (300 / n, digits=1),
                                    round (350 / n, digits=1),
                                    round (400 / n, digits=1),
                                    round (450 / n, digits=1),
                                    round (500 / n, digits=1),
                                    round (550 / n, digits=1),
                                    round (600 / n, digits=1),
                                    round (650 / n, digits=1),
                                    round (700 / n, digits=1),
                                    round (750 / n, digits=1),
                                    round (800 / n, digits=1),
                                    round (850 / n, digits=1),
                                    round (900 / n, digits=1),
                                    round (950 / n, digits=1),
                                    round (1000 / n, digits=1)),
                        expand = c (0, 0)) +
    geom_point(aes(color = zones), size = pointsize, alpha = pointalpha) +
    geom_segment (aes (x = x1, y = y1, xend = xend, yend = yend),
                  data=border, linetype=linetype) +
    annotate (geom="text", x = labels$x, y = labels$y, size = 6,
              label = labels$label, color=labels$color) +
    theme_bw () +
    #    theme (legend.position = "none") +
    ggtitle (title) +
    xlab (xlab) +
    ylab (ylab)
  ceg
}

plotParkesGrid <- function (referenceVals, testVals, type=1, title="", xlab="", ylab="", linesize = 0.5, 
                            linetype = "solid", linecolor = "black", linealpha = 0.6, pointsize = 2, pointalpha = 1, zones = NA, unit='gram') {
  
  if (type != 1 & type != 2) {
    stop("'type' must be 1 or 2.")
  }
  if (unit != "mol" & unit != "gram") { 
    stop("'unit' must be either 'mol' or 'gram'.")
  }
  if (title==""){ 
    title <- paste("Parkes (Consensus) Error Grid for Type", type, "Diabetes")
  }
  if (unit == "mol") { 
    n <- 18 
    if (xlab==""){
      xlab="Reference Glucose Concentration (mmol/L)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mmol/L)"
    }
  } else {
    n <- 1
    if (xlab==""){
      xlab="Reference Glucose Concentration (mg/dL)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mg/dL)"
    }
  }
  
  if (is.na (zones)) {
    zones <- getParkesZones (referenceVals, testVals, type, unit)
  }
  
  data <- data.frame (ref=referenceVals, test=testVals, zones=zones)
  
  
  maxX <- max (max (data$ref) + 20 / n, 550 / n)#better solution
  maxY <- max (max (data$test) + 20 / n, 550 / n)
  
  labels <- data.frame (x=c (320 / n, 220 / n, 385 / n, 140 / n, 405 / n, 415 / n, 75 / n, 21 / n), 
                        y=c (320 / n, 360 / n, 235 / n, 375 / n, 145 / n,  50 / n, 383 / n, 383 / n), 
                        label=c ("A", "B", "B", "C", "C", "D", "D", "E"),
                        color=c ("green", "blue", "blue", "red", "red", "red", "red", "red"))
  
  
  if (type==1){
    ce <- coef (35, 155, 50, 550)
    cdu <- coef (80, 215, 125, 550)
    cdl <- coef (250, 40, 550, 150)
    ccu <- coef (70, 110, 260, 550)
    ccl <- coef (260, 130, 550, 250)
    cbu <- coef (280, 380, 430, 550)
    cbl <- coef (385, 300, 550, 450)
    border <- data.frame (x1=c (0 / n, 
                                0 / n, 30 / n, 140 / n, 280 / n,
                                50 / n, 50 / n, 170 / n, 385 / n,
                                0 / n, 30 / n, 50 / n, 70 / n,
                                120 / n, 120 / n, 260 / n,
                                0 / n, 25 / n, 50 / n, 80 / n,
                                250 / n, 250 / n,
                                0 / n, 35 / n),
                          y1=c (0 / n, 
                                50 / n, 50 / n, 170 / n, 380 / n,
                                0 / n, 30 / n, 145 / n, 300 / n,
                                60 / n, 60 / n, 80 / n, 110 / n,
                                0 / n, 30 / n, 130 / n,
                                100 / n, 100 / n, 125 / n, 215 / n,
                                0 / n, 40 / n,
                                150 / n, 155 / n),
                          xend=c (min (maxX, maxY), 
                                  30 / n, 140 / n, 280 / n, endx (280 / n, 380 / n, maxY, cbu),
                                  50 / n, 170 / n, 385 / n, maxX,
                                  30 / n, 50 / n, 70 / n, endx (70 / n, 110 / n, maxY, ccu),
                                  120 / n, 260 / n, maxX,
                                  25 / n, 50 / n, 80 / n, endx (80 / n, 215 / n, maxY, cdu),
                                  250 / n, maxX,
                                  35 / n, endx (35 / n, 155 / n, maxY, ce) ),
                          yend=c (min (maxX, maxY), 
                                  50 / n, 170 / n, 380 / n, maxY,
                                  30 / n, 145 / n, 300 /n, endy (385 / n, 300 / n, maxX, cbl),
                                  60 / n, 80 / n, 110 / n, maxY,
                                  30 / n, 130 / n, endy (260 / n, 130 / n, maxX, ccl),
                                  100 / n, 125 / n, 215 / n, maxY,
                                  40 / n, endy (410 / n, 110 / n, maxX, cdl),
                                  155 / n, maxY))
    
  } else {#type 2
    ce <- coef (35, 200, 50, 550)
    cdu <- coef (35, 90, 125, 550)
    cdl <- coef (410, 110, 550, 160)
    ccu <- coef (30, 60, 280, 550)
    ccl <- coef (260, 130, 550, 250)
    cbu <- coef (230, 330, 440, 550)
    cbl <- coef (330, 230, 550, 450)
    
    border <- data.frame (x1=c (0 / n, 
                                0 / n, 30 / n, 230 / n,
                                50 / n, 50 / n, 90 / n, 330 / n,
                                0 / n, 30 / n,
                                90 / n, 260 / n,
                                0 / n, 25 / n, 35 / n,
                                250 / n, 250 / n, 410 / n, 
                                0 / n, 35 / n),
                          y1=c (0 / n, 
                                50 / n, 50 / n, 330 / n,
                                0 / n, 30 / n, 80 / n, 230 / n,
                                60 / n, 60 / n,
                                0 / n, 130 / n,
                                80 / n, 80 / n, 90 / n,
                                0 / n, 40 / n, 110 / n, 
                                200 / n, 200 / n),
                          xend=c (min (maxX, maxY), 
                                  30 / n, 230 / n, endx (230 / n, 330 / n, maxY, cbu),
                                  50 / n, 90 / n, 330 / n, maxX,
                                  30 / n, endx (30 / n, 60 / n, maxY, ccu),
                                  260 / n, maxX,
                                  25 / n, 35 / n, endx (35 / n, 90 / n, maxY, cdu),
                                  250 / n, 410 / n, maxX,
                                  35 / n, endx (35 / n, 200 / n, maxY, ce) ),
                          yend=c (min (maxX, maxY), 
                                  50 / n, 330 / n, maxY,
                                  30 / n, 80 / n, 230 / n, endy (330 / n, 230 / n, maxX, cbl),
                                  60 / n, maxY,
                                  130 / n, endy (260 / n, 130 / n, maxX, ccl),
                                  80 / n, 90 / n, maxY,
                                  40 / n, 110 / n, endy (410 / n, 110 / n, maxX, cdl),
                                  200 / n, maxY))
    
  }
  
  ref <- test <- NULL
  peg <- ggplot(data, aes(x = ref, y = test)) + 
    scale_x_continuous (breaks = c (round (70 / n, digits=1), 
                                    round (100 / n, digits=1), 
                                    round (150 / n, digits=1), 
                                    round (180 / n, digits=1), 
                                    round (240 / n, digits=1), 
                                    round (300 / n, digits=1), 
                                    round (350 / n, digits=1), 
                                    round (400 / n, digits=1), 
                                    round (450 / n, digits=1), 
                                    round (500 / n, digits=1), 
                                    round (550 / n, digits=1), 
                                    round (600 / n, digits=1), 
                                    round (650 / n, digits=1), 
                                    round (700 / n, digits=1), 
                                    round (750 / n, digits=1), 
                                    round (800 / n, digits=1), 
                                    round (850 / n, digits=1), 
                                    round (900 / n, digits=1), 
                                    round (950 / n, digits=1), 
                                    round (1000 / n, digits=1)), expand = c (0, 0)) + 
    scale_y_continuous (breaks = c (round (70 / n, digits=1), 
                                    round (100 / n, digits=1), 
                                    round (150 / n, digits=1), 
                                    round (180 / n, digits=1), 
                                    round (240 / n, digits=1), 
                                    round (300 / n, digits=1), 
                                    round (350 / n, digits=1), 
                                    round (400 / n, digits=1), 
                                    round (450 / n, digits=1), 
                                    round (500 / n, digits=1), 
                                    round (550 / n, digits=1), 
                                    round (600 / n, digits=1), 
                                    round (650 / n, digits=1), 
                                    round (700 / n, digits=1), 
                                    round (750 / n, digits=1), 
                                    round (800 / n, digits=1), 
                                    round (850 / n, digits=1), 
                                    round (900 / n, digits=1), 
                                    round (950 / n, digits=1), 
                                    round (1000 / n, digits=1)), expand = c (0, 0)) + 
    geom_point(aes(color = zones), size=pointsize, alpha=pointalpha) +  
    geom_segment (aes (x=x1, y=y1, xend=xend, yend=yend), data=border, linetype=linetype) +
    annotate (geom="text", x=labels$x, y=labels$y, size=6, label=labels$label, color=labels$color) +
    theme_bw () + 
    #    theme (legend.position = "none") + 
    ggtitle (title) + 
    xlab (xlab) + 
    ylab (ylab)
  peg
}

#}##comment out end

# Input mehtod choices
inputChoices <- c ("Test grid", "Example dataset", "File (.csv) upload", "Manual data entry")

# Data measured in the following units
inputUnits <- c ("mg/dL", "mmol/L")

#glucose data - example dataset
glucose_data <- readRDS (file="glucose_data.Rda")

# which fields are mandatory
fieldsMandatory <- c("ref", "test")

# add an asterisk to an input label
labelMandatory <- function(label) {
  tagList(
    label,
    span("*", class = "mandatory_star")
  )
}


##################
ui <- dashboardPage( 
  dashboardHeader(title = strong (em ("Error Grid Analysis"))), 
  dashboardSidebar( 
    sidebarMenu( 
      menuItem("Data Entry", tabName = "input", icon=icon("table")),
      menuItem("Data Analysis", tabName = "analysis", icon=icon("cogs")),
      menuItem("Result Interpretation", tabName = "interpretation", icon=icon("commenting-o")),
      menuItem("Export Data", tabName = "export", icon=icon("sign-out"))
    ),
    br (),
    p(class = "text-muted",
      br(),
      "Computation with EGA package"
    ),
    br (),
    fluidRow(
      column (12,
        textOutput("currentTime", container = span, inline = TRUE)
      )
    ),
    br (),
    br (),
    p(class = "text-muted",
      "Code by Sergei Mihhailov"
    )
  ), 
  dashboardBody( 
    tabItems( 
      tabItem("input", 
        box (width=12, title="Input", status = "warning", solidHeader=TRUE,
          column (12,
            radioButtons (inputId="inputMethod", label="Select input method", 
              width="100%", choices=inputChoices, inline=TRUE,
              selected=inputChoices[1]
            ),
            radioButtons (inputId="inputUnit", label="Select the unit your data is measured in", 
                          width="100%", choices=inputUnits, inline=TRUE,
                          selected=inputUnits[1]
            )
          )
        ),
        box (width=12, title="Details", status="info", solidHeader=TRUE,
          fluidRow(
            column (12, offset=0,
              htmlOutput("inputComment", container=span, inline=TRUE)
            )
          ),
          conditionalPanel(condition="input.inputMethod == 'File (.csv) upload'",       
            fileInput('file1', 'Choose CSV File', 
              accept=c('text/csv', 
              'text/comma-separated-values,text/plain', '.csv'
              )
            ),
            tags$hr (),
            fluidRow (
              column (4,
                checkboxInput ('header', 'Header', TRUE),
                radioButtons ('mark', 'Decimal Mark', c ('Point'='.', 'Comma'=','))
              ),
              column (4,
                radioButtons ('sep', 'Separator',
                  c (Comma=',', Semicolon=';', Tab='\t'),
                ','
                )
              ),
              column (4, 
                radioButtons('quote', 'Quote', 
                  c('None'='', 'Double Quote'='"', 'Single Quote'="'"), ''
                )
              )
            )
          ),
          conditionalPanel(condition="input.inputMethod == 'Manual data entry'",       
            #input fields
            box(title = "Enter Data", status = 'info',
                #rHandsontableOutput("hot")
                aceEditor("responses", value="ref\ttest\n1\t1\n2\t2.2\n3\t3.1\n4\t3.6\n5\t3.5", mode="r", theme="cobalt")
            )
          )
        ),
        box (width=12, title="Dataset", status="primary", solidHeader=TRUE,
          p ("Below you can see the selected dataset."),
          DT::dataTableOutput('rawdata')
        )
      ),
      tabItem("analysis",
        box (width=12, title="Input overview", status = "warning", solidHeader=TRUE,
          column (3,
            radioButtons (inputId="inputGrid", label="Select Grid", 
                                   width="100%", choices=c ("Clarke", "Parkes (consensus) Type 1 Diabetes", 
                                                            "Parkes (consensus) Type 2 Diabetes"), inline=FALSE,
              selected="Parkes (consensus) Type 1 Diabetes"
            )
          ),
          column (9,
            valueBoxOutput("inputbox"),
            valueBoxOutput("unitbox"),
            valueBoxOutput("lengthbox")
          )
        ), 
        box (width=12, title="Error Grid Analysis Plot", status="primary", solidHeader=TRUE,
             plotlyOutput("egaPlot", height=650)
        )
        
      ), 
      
      tabItem("interpretation", 
        box (width=12, title="Interpretation", status="primary", solidHeader=TRUE,
          column (12,
            valueBoxOutput("abzonebox"),
            valueBoxOutput("accuracybox"),
            valueBoxOutput("lengthokbox")
          ), 
          h1("ISO 15197:2013 accuracy criteria"),
          HTML ("<p>Current In vitro diagnostic test systems - Requirements for
              blood-glucose monitoring systems for self-testing in
              managing diabetes mellitus (ISO 15197:2013) standard lays down 
              new requirements on accuracy of glucose meters.</p>"),
          HTML ("<p>Minimum system accuracy performance criterium regarding Consensus Error grid is
             <b>99 % of individual glucose measured values shall fall within zones A and B</b> of 
              the Consensus Error Grid (CEG) for type 1 diabetes.</p>"),
          DT::dataTableOutput('zones'),
          htmlOutput("absum", container = span, inline = TRUE),
          HTML ("<p><b>95 %</b> of the measured glucose values shall fall within either &plusmn; 0,83 mmol/L (&plusmn; 15 mg/dL)
              of the average measured values of the reference measurement procedure at glucose 
              concentrations < 5,55 mmol/L (<100 mg/dL) or within &plusmn; 15 % at glucose concentrations 
             ā„ 5,55 mmol/L (ā„100 mg/dL).</p>"),
          htmlOutput("accuracyComment", container = span, inline = TRUE),
          HTML ("<p>There should be minimum 100 subjects studied with dublicte measurements taken 
             with at least 3 reagent lots (i.e. <b>600 measurements</b>).<p>"),
          htmlOutput("lengthComment", container = span, inline = TRUE)
        )
      ),
      tabItem("export",
        h1 ("Plot Export"),
        HTML ("Plot can be exported as a .png image file through context menu in 'Data Analysis' tab (on the left of this screen).
              There you can select a specific part of the plot (zoom in/out, pan, select etc.) or choose a subset of data points
              from selected zones (click on zones legend)."),
        h1 ("Zones export in a dataset"),
        HTML ("You can download a data source selected in 'Data Entry' tab (on the left of this screen) with additional column of zones
              based on your grid selection (in menu 'Grid Select' on 'Data Analysis' tab). To do so just press 'Dataset Download' 
              button below. The data is in .csv format, field separator is ',', fraction part separator is '.'"),
        br(),
        downloadButton ("downloadCsv", "Dataset Download")
        
      )
    )
  ), skin='red'
)


#################
server <- function(input, output, session) {
  
  dataset <- reactive ({
    if (unit () == "mol"){
      n <- 18#scaling factor for mmol/L compared to mg/dL
    } else {
      n <- 1
    }
    if (input$inputMethod==inputChoices[1]){
      ref <- seq (1, 1000, length.out=50)
      test <- seq (1, 1000, length.out=50)
      dataset <- expand.grid (ref=ref / n, test=test / n)
    } else if (input$inputMethod==inputChoices[2]){
      dataset <- glucose_data / n
    } else if (input$inputMethod==inputChoices[3]){
      inFile <- input$file1
      
      if (is.null(inFile))
        return("Dataset not loaded yet")
      
      dataset <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                                 quote=input$quote)
    } else {
      dataset <- read.csv(text=gsub("\\,", "\\.", input$responses), 
                          sep="", 
                          na.strings=c("","NA","."))
    }
  })
  
  unit <- reactive ({
    if (input$inputUnit=="mmol/L"){
      return ("mol")
    } else {
      return ("gram")
    }
  })
  
  refdata <- reactive ({
    data <- dataset ()
    if (length (data$ref) > 0){
      return (as.vector (data$ref))
    } else {
      return ("NULL")
    }
  })
  
  testdata <- reactive ({
    data <- dataset ()
    if (length (data$test) > 0){
      return (as.vector (data$test))
    } else {
      return ("NULL")
    }
  })
  
  #clock to show
  output$currentTime <- renderText({
    invalidateLater(1000, session)
    paste ("Current time: ", Sys.time ())
  })

  output$inputComment <- renderText({
    if (input$inputMethod == inputChoices[1]){
      ("<p>Test grid consists of simulated data ranged from 0 to 1000 mg/dL (0 to 55.56 mmol/L respectively). 
       It contains all possible combinations of reference method glucose results (column 'ref') and test method results (column 'test').</p> 
       <p>This choice is reccomended to get the feeling of this application and evaluate its accuracy. 
       Examine if zones are labelled correctly near the borders.</p>
       <p>Press 'Data Analysis' tab on the left to see the plot.</p>")
    } else if (input$inputMethod == inputChoices[2]){
      ("<p>Example dataset is glucose_data from package {ega}. It contains 5072 paired paired reference and test glucose values. 
        Reference method glucose value, in mg/dL or mmol/L, depending on your choice in input parameters is in column 'ref'.
        Test method glucose value, in mg/dL or mmol/L, depending on your choice in input parameters is in column 'test'.</p>
       <p>The data originates from a modified clinical dataset.</p>
       <p>Press 'Data Analysis' tab on the left to see the plot.</p>")
    } else if (input$inputMethod == inputChoices[3]){
      ("<p>The simpliest way to enter user data is by uploading it in plain csv format. It is easy to make one in MS Excel, LibreOffice 
       or similar software (or in a notepad). The data must be entered in two columns 'ref' for reference method results and 'test' for
       respective test method result. This application provides wide flexibility for data field and fraction part separators. Choose a file
       on your local drive for upload and look at details below.</p>
       <p>Press 'Data Analysis' tab on the left to see the plot.</p>")
    } else {
      ("<p>In this windows you can manually enter the data for analysis.</p>
        <p></p>
       <p><b>Note</b>: in case you run this application in RStudio this function stores the data in 'responses' dataset 
       in global environment. You can upload your data into this application by storing your dataset with the name 'responses' 
       (i.e. responses <- your_glucose_data).</p>
       <p>Press 'Data Analysis' tab on the left to see the plot.</p>")
    }
  })

  output$rawdata <- DT::renderDataTable({ 
    datatable (dataset(), 
               options = list(
                 dom = 'tip',
                 pageLength = 5
                 )
    )
  }) 
  

  #input method selected
  output$inputbox <- renderValueBox({
    valueBox(
      value = input$inputMethod,
      subtitle = paste ("Data source selected"),
      icon = icon("database"),
      color = "green"
    )
  })
  
  #units selected
  output$unitbox <- renderValueBox({
    valueBox(
      value = input$inputUnit,
      subtitle = paste ("Units selected"),
      icon = icon("hand-lizard-o"),
      color = "green"
    )
  })
  
  #number of comparisons in dataset
  output$lengthbox <- renderValueBox({
    valueBox(
      value = length(dataset()$ref),
      subtitle = paste ("Number of comparisons in the dataset"),
      icon = icon("arrows-v"),
      color = "green"
    )
  })
  
  output$egaPlot <- renderPlotly({
    if (input$inputGrid=="Clarke"){
      ggplotly (plotClarkeGrid (referenceVals=refdata (), testVals=testdata (), unit=unit()))
    } else if (input$inputGrid=="Parkes (consensus) Type 1 Diabetes"){
      ggplotly (plotParkesGrid (referenceVals=refdata (), testVals=testdata (), unit=unit()))
    } else {
      ggplotly (plotParkesGrid (referenceVals=refdata (), testVals=testdata (), unit=unit(), type=2))
    }
  })
  
  zones <- reactive ({
    zones <- getParkesZones (refdata (), testdata (), unit(), type=1)
    zones <- factor (zones)
    return (data.frame (round (table (zones) / length (zones) * 100, digits=2)))
  })
  
  output$zones <- DT::renderDataTable({
    datatable (zones (), 
      options = list (
        pageLength=6,
        dom = 't',#no extra functionality
        rownames = FALSE
      ),
      caption="Results for the current dataset Parkes (consensus) Type 1 Diabetes Error Grid"
    )
  })
  
  accuracy <- reactive ({
    data <- dataset ()
    if (unit ()=="gram"){
      datalow <- data [data$ref < 100, ]
      outlow <- (which ( (datalow$test > datalow$ref + 15) | 
                         (datalow$test < datalow$ref - 15)))
      datahigh <- data [data$ref >= 100, ]
      outhigh <- (which ( (datahigh$test > datahigh$ref * 1.15) | 
                           (datahigh$test < datahigh$ref * 0.85)))
      return (100 - (length (outlow) + length (outhigh))/ length (data$ref) * 100)
    } else {#mmol/l
        datalow <- data [data$ref < 5.55, ]
        outlow <- (which ( (datalow$test > datalow$ref + 0.83) | 
                             (datalow$test < datalow$ref - 0.83)))
        datahigh <- data [data$ref >= 5.55, ]
        outhigh <- (which ( (datahigh$test > datahigh$ref * 1.15) | 
                              (datahigh$test < datahigh$ref * 0.85)))
        return (100 - (length (outlow) + length (outhigh))/ length (data$ref) * 100)
    }
  })
  
  #share of datapoints in zones A and B adequate?
  output$abzonebox <- renderValueBox({
    valueBox(
      value = if (abS () > 99) {
        paste0 ("OK: ", round (abS (), digits=1), "%")
      } else {
        paste0 ("NOK: ", round (abS (), digits=1), "%")
      },
      subtitle = paste ("Data points in zones A and B"),
      icon = if (abS () > 99) icon("check-square") else icon ("exclamation-circle"),
      color = if (abS () > 99) "green" else "red"
    )
  })
  
  #95% within required accuracy?
  output$accuracybox <- renderValueBox({
    valueBox(
      value = if (accuracy () > 95) {
        paste0 ("OK: ", round (accuracy (), digits=1), "%")
      } else {
        paste0 ("NOK: ", round (accuracy (), digits=1), "%")
      },
      subtitle = paste ("Within required accuracy"),
      icon = if (accuracy () > 95) icon("check-square") else icon ("exclamation-circle"),
      color = if (accuracy () > 95) "green" else "red"
    )
  })
  
  #number of comparisons in dataset adequate?
  output$lengthokbox <- renderValueBox({
    valueBox(
      value = if (length(dataset()$ref)>=600) {
        paste ("OK:", length (dataset ()$ref))
        } else {
          paste ("NOK:", length (dataset ()$ref))
        },
      subtitle = paste ("Number of comparisons in the dataset"),
      icon = if (length(dataset()$ref)>=600) icon("check-square") else icon ("exclamation-circle"),
      color = if (length(dataset()$ref)>=600) "green" else "red"
    )
  })
  
  
  #counting points in zones A and B
  abS <- reactive ({#sum of zone A and B
    zone <- getParkesZones (refdata (), testdata (), unit(), type=1)
    zones <- factor (zone)
    table <- ((table (zones) / length (zones) * 100))
    a <- tryCatch(#handling an error when A or B zones are empty
      {
        table[["A"]] 
      },
      error=function(cond) {
        return (0)
      }
    ) 
    
    b <- tryCatch(
      {
        table[["B"]] 
      },
      error=function(cond) {
        return (0)
      }
    )    
    return (a + b)
    
  })
  
  #interpretation comment if datapoints in A and B zones is enough for validation
  output$absum <- renderText(
    if (abS () > 99){
      paste ("The summ of A and B zone percentages is<b>", round (abS (), digits=2), "</b>thus result is <b>acceptable</b>.<br>")
    } else {
      paste ("The summ of A and B zone percentages is<b>", round (abS (), digits=2), "</b>thus result is <b>unacceptable</b>.<br>")
    }
  )
  
  #interpretation comment if datapoints in A and B zones is enough for validation
  output$accuracyComment <- renderText(
    if (accuracy () > 95){
      paste ("This parameter for the current dataset is<b>", round (accuracy (), digits=2), "%</b> thus result is <b>acceptable</b>.<br>")
    } else {
      paste ("This parameter for the current dataset is<b>", round (accuracy (), digits=2), "%</b> thus result is <b>unacceptable</b>.<br>")
    }
  )
  
  #interpretation comment if datapoints in A and B zones is enough for validation
  output$lengthComment <- renderText(
    if (length (dataset ()$ref) >= 600){
      paste ("The number of reference method and test method comparisons is<b>", length (dataset ()$ref), 
             "</b>which is <b>acceptable</b>.<br>")
    } else {
      paste ("The number of reference method and test method comparisons is <b>", length (dataset ()$ref), 
            "</b>which is <b>not acceptable</b>.<br>")
    }
  )
  
  #input dataset with additional column vector of zones depending on selected grid
  egaout <- reactive({
      if (input$inputGrid == "Clarke"){
        out <- cbind (dataset(), getClarkeZones (refdata (), testdata (), unit ()))
        colnames (out) <- c ("Reference", "Test", "Clarke zone")
        return (out)
      } else if (input$inputGrid == "Parkes (consensus) Type 1 Diabetes"){
        out <- cbind (dataset(), getParkesZones (refdata (), testdata (), unit (), type=1))
        colnames (out) <- c ("Reference", "Test", "Parkes T1D zone")
        return (out)
      } else {
        out <- cbind (dataset(), getParkesZones (refdata (), testdata (), unit (), type=2))
        colnames (out) <- c ("Reference", "Test", "Parkes T2D zone")
        return (out)
      }
  })

  
  #data export as csv
  output$downloadCsv <- downloadHandler( 
    filename <- "ErrorGridAnalysis.csv",
    content = function (file) { 
      write.csv (egaout (), file) 
    }, 
    contentType = "text/csv" 
  ) 
  
}

shinyApp(ui = ui, server = server)
