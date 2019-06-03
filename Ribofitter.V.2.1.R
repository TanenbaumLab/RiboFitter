############################################################################################################
## SETTING UP THE MODEL Frame1; Model Parameters
############################################################################################################

Buildup_time_Frame1 <- 428.88 # Time it takes to translate the epitope array, in seconds
Plateau_time_Frame1 <- 74.19 # Time it takes to translate the remainder of the transcript, in seconds
single_ribosome_intensity <- 110.3 # Peak intensity of a single ribosome
acquisition_interval <- 30  #  Number of seconds between image acquisition in time-lapse
frames_per_minute <- 60/acquisition_interval

############################################################################################################
## SETTING UP THE MODEL Frame2; Model Parameters
############################################################################################################

Buildup_time_Frame2 <- 428.88# time it takes to translate the epitope array, in seconds
Plateau_time_Frame2 <- 74.65# time it takes to translate the remainder of the transcript, in seconds
single_ribosome_intensity_Frame2 <- 132.68# peak intensity of a single ribosome

############################################################################################################
## SETTING THE COMPUTATIONAL PARAMETERS 
############################################################################################################

parallel_runs <- 10         # Number of times the algorithm is run on the same signal
niterations <- 1000         # Number of iterations per run
change_probability <- 0.3   # The chance of a ribosome to relocate during each iteration 
change_size <- 25           # The distance a ribosome moves 
add_probability <- 0.1      # The chance to add a ribosome to the current fit during each iteration
remove_probability <- 0.1   # The chance to remove a ribosome to the current fit during each iteration 

############################################################################################################
bestfit=data.frame(mRNA=character(0),run=numeric(0),number_of_ribosomes=numeric(0),initiation_times=character(0),error=numeric(0),group=numeric(0))
############################################################################################################
##Included Functions:
############################################################################################################

list.as.matrix <- function(x, byrow=FALSE, filler=NA){
  maxlen <- max(sapply(x,length))
  xm <- sapply(x, 
               function(xs){fillen <- maxlen-length(xs)
               if (fillen>0) {c(xs,rep(filler,fillen))} else xs
               })
  if (byrow)return(t(xm)) else return(xm)
}

theoretical_signal <- function(initiation_time){ # initiation_time in seconds
  signal <- c(rep(0,initiation_time),    # before translation
              seq(0,1,1/Buildup_time_Frame1),    # translating the Frame1 array
              rep(1,Plateau_time_Frame1),    # translating the remainder of the transcript
              rep(0,length(timepoints_Frame1))  # after translation
  )
  return(single_ribosome_intensity*signal[1:length(timepoints_Frame1)])
}

update <- function(times, maxtime, change_probability, change_size, add_probability, remove_probability){
  updated_times <- times
  added_values = rnorm(length(times),mean=0,sd=change_size) * (runif(length(times)) < change_probability) # add a Gaussian value with probability change_probability 
  if(all(updated_times+added_values>=0) && all(updated_times+added_values<=maxtime)){ # check if proposed times are not too small or too large
    updated_times <- updated_times + added_values  
  }
  if(runif(1) < add_probability){  # add a new value to the list of initiation times
    x = sample(1:maxtime,size=1)
    updated_times <- sort(c(updated_times,x)) 
  } else if(runif(1) < remove_probability){  # remove a value from the list of initiation times
    number_of_mRNAs <- length(updated_times)
    if(number_of_mRNAs > 1){
      updated_times <- sort(sample(updated_times,size=number_of_mRNAs-1))
    }
  }
  return(updated_times)
}


# Functions to calculate the error and the RMS error between two consecutive 
calculate_error <- function(signal1,signal2){ 
  return(sum((signal1-signal2)^2))
}
calculate_RMS <- function(signal1,signal2){  # calculate the RMS error between two signals
  return(sqrt(calculate_error(signal1,signal2)/length(signal1)))
}


############################################################################################################
## LOADING THE Frame1 AND Frame2 DATA 
############################################################################################################
select_Frame1=file.choose(new = FALSE)
Frame1_Data <- read.csv(select_Frame1,header=T,stringsAsFactors=F,sep=";") # Read Frame1 Data 
Frame1_Data=Frame1_Data[rowSums(is.na(Frame1_Data)) != ncol(Frame1_Data),]

select_Frame2=tryCatch(file.choose(), error = function(e) "") 
if(select_Frame2==""){Frame2_Data=Frame1_Data[1,]
Frame2_Data[Frame2_Data!=0]<-0}
if(select_Frame2!=""){Frame2_Data <- read.csv(select_Frame2,header=T,stringsAsFactors=F, sep=";") # Read Frame2 Data
Frame2_Data=Frame2_Data[rowSums(is.na(Frame2_Data)) != ncol(Frame2_Data),]}

dif=(ncol(Frame1_Data)-ncol(Frame2_Data))

if(dif<0){Frame1_Data=cbind(Frame1_Data,matrix(rep(0,abs(dif)*nrow(Frame1_Data)),ncol=abs(dif)))
          colnames(Frame1_Data)=colnames(Frame2_Data)              }
if(dif>0){Frame2_Data=cbind(Frame2_Data,matrix(rep(0,dif*nrow(Frame2_Data)),ncol=dif))
          colnames(Frame2_Data)=colnames(Frame1_Data)}


############################################################################################################
## LOADING AND ADJUSTING Frame1_Data
############################################################################################################
fitwindow_Frame1=data.frame(seq(from=-round(Buildup_time_Frame1/acquisition_interval)+1,to=0,by=1),matrix(0,nrow=round(Buildup_time_Frame1/acquisition_interval),ncol=ncol(Frame1_Data)-1)) # Add zeros to the Frame1_Data-frame to allow fitting at the start of image acquisition. 
colnames(fitwindow_Frame1)<-colnames(Frame1_Data)
Frame1_Data=rbind(fitwindow_Frame1,Frame1_Data)

total_minutes_Frame1 <- nrow(Frame1_Data)/(60/acquisition_interval)    # Length of observation window (in minutes)
timepoints_Frame1 <- seq(1,total_minutes_Frame1*60)  # All timepoints (in seconds)

noise_level <- 1 # determines the acceptance probability

total_intensity_per_ribosome <- sum(theoretical_signal(0)) # Calculate the total intensity of a single translating ribosome based on the theoretical signal described above. 

mRNAnames <- colnames(Frame1_Data)[2:ncol(Frame1_Data)]
last_timepoint=integer()

for(i in 1:length(mRNAnames)){
  if(length(which(is.na(Frame1_Data[,mRNAnames[i]])==TRUE))==0){
    second.value=length(Frame1_Data[,mRNAnames[i]])+1
  }
  else if(length(split(which(is.na(Frame1_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame1_Data[,mRNAnames[i]]))) != 1))))==1 && max(split(which(is.na(Frame1_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame1_Data[,mRNAnames[i]]))) != 1)))$`1`)<length(Frame1_Data[,mRNAnames[i]])){
    second.value=length(Frame1_Data[,mRNAnames[i]])
  }
  else if(length(split(which(is.na(Frame1_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame1_Data[,mRNAnames[i]]))) != 1))))==1 && max(split(which(is.na(Frame1_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame1_Data[,mRNAnames[i]]))) != 1)))$`1`)==length(Frame1_Data[,mRNAnames[i]])){
    second.value=split(which(is.na(Frame1_Data[,mRNAnames[i]])), cumsum(c(1,diff(which(is.na(Frame1_Data[,mRNAnames[i]])))!=1)))
    second.value=second.value$`1`[1]
  }
  else{second.value=split(which(is.na(Frame1_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame1_Data[,mRNAnames[i]]))) != 1)))
  second.value=second.value$`2`[1]}
  last_timepoint=c(last_timepoint,second.value)
}

############################################################################################################
# ACTUAL SIMULATION
############################################################################################################

all_predictions <- list() # Final list of the initiation times 
error_list <- c() # Vector of the errors between the signal and the fit 

for(k in 1:length(mRNAnames)){  # select some mRNAs as an example
  observed_timepoints <- timepoints_Frame1[seq(acquisition_interval,(last_timepoint[k]-1)*acquisition_interval,acquisition_interval)] # the timepoints at which a signal has been recorded
  maxtime=(last_timepoint[k]-1)*acquisition_interval # The last time point (in seconds) at which the ribosome could initiate translation
  
  ## Start with the fluorescence signal
  real_intensity <- Frame1_Data[,mRNAnames[k]]  
  real_intensity[is.na(real_intensity)] <- 0        # Change NA-values to zeros    
  real_intensity[real_intensity<0]<-0               # Change negative intensity values for zeros 
  real_intensity[c(1:last_timepoint[k]-1)]          # Select the timewindow of interest from 1 to last_timepoint-1  
  
  real_intensity_copy <-real_intensity
  real_intensity_copy <-rep(real_intensity_copy,each=acquisition_interval)
  
  real_intensity <- rep(real_intensity,each=acquisition_interval)
  ## Estimate the number of ribosomes from the total intensity   
  estimated_number_of_ribosomes <- max(c(0,round(sum(real_intensity_copy)/total_intensity_per_ribosome)))
  
  timepoints=1:length(real_intensity)
  #If the number of estimated number of ribosomes is less than one, do not proceed with further fitting.  
  if(estimated_number_of_ribosomes==0){ prediction <- list(mRNA = k,
                                                           run = 0,
                                                           number_of_ribosomes = 0,
                                                           initiation_times = 0,
                                                           error = allRMSs_nulfit)
  all_predictions[[length(all_predictions)+1]] <- prediction
  }  
  
  else{
    ## Randomly initialize the set of initiation times.
    current_estimated_times_list <- replicate(parallel_runs,
                                              sort(sample(timepoints,size=estimated_number_of_ribosomes,
                                                          prob=c(real_intensity[(nrow(fitwindow_Frame1)*acquisition_interval+1):(2*nrow(fitwindow_Frame1)*acquisition_interval)]/sum(real_intensity[(nrow(fitwindow_Frame1)*acquisition_interval+1):length(timepoints)]),    real_intensity[(nrow(fitwindow_Frame1)*acquisition_interval+1):length(timepoints)]/sum(real_intensity[(nrow(fitwindow_Frame1)*acquisition_interval+1):length(timepoints)])))),
                                              simplify=F)
    current_estimated_times_list
    
    ## Run the algorithm 
    for(i in 0:niterations){
      current_estimated_times_list <- lapply(current_estimated_times_list,function(current_estimated_times){
        # Slightly alter initiation times
        new_estimated_times <- sort(update(current_estimated_times,maxtime,change_probability, change_size, add_probability, remove_probability))  
        # Calculate theoretical intensities and errors for both the old and the new initiation times
        current_theoretical_intensity <- rowSums(sapply(current_estimated_times,theoretical_signal))
        new_theoretical_intensity <- rowSums(sapply(new_estimated_times,theoretical_signal))
        current_error <- calculate_error(real_intensity[observed_timepoints[nrow(fitwindow_Frame1):length(observed_timepoints)]],current_theoretical_intensity[observed_timepoints[nrow(fitwindow_Frame1):length(observed_timepoints)]])
        new_error <- calculate_error(real_intensity[observed_timepoints[nrow(fitwindow_Frame1):length(observed_timepoints)]],new_theoretical_intensity[observed_timepoints[nrow(fitwindow_Frame1):length(observed_timepoints)]])
        
        # Calculate the acceptance probability 
        acceptance_probability <- exp((current_error - new_error)/(2*noise_level^2))
        
        # With probability acceptance_probability, keep the new initiation times 
        if(runif(1) < acceptance_probability){
          return(new_estimated_times)
        }
        else{
          return(current_estimated_times)
        }
      })
      
      if(i %% 100 == 0){
        print(paste(c("mRNA ",k,"; iteration ",i,"/",niterations),collapse=""))
      }
    }
    
    
    ## Calculate the RMS error for each of the found solutions
    all_RMSs <- sapply(1:parallel_runs,
                       function(x){calculate_RMS(rowSums(sapply(current_estimated_times_list[[x]],theoretical_signal))[observed_timepoints[1:length(observed_timepoints)]],
                                                 real_intensity[observed_timepoints])})
    
    allRMSs_nulfit <- calculate_RMS(real_intensity[observed_timepoints],rep(0,length(observed_timepoints)))
    
    
    
    
    ## Store the result
    if(min(all_RMSs)<allRMSs_nulfit){
      for(x in 1:parallel_runs){
        
        
        
        prediction <- list(mRNA = k,
                           run = x,
                           number_of_ribosomes = length(current_estimated_times_list[[x]]),
                           initiation_times = current_estimated_times_list[[x]],
                           error = all_RMSs[x])
        all_predictions[[length(all_predictions) + 1]] <- prediction
      }
      
      
    }
    else{
      prediction <- list(mRNA = k,
                         run = 0,
                         number_of_ribosomes = 0,
                         initiation_times = 0,
                         error = allRMSs_nulfit)
      all_predictions[[length(all_predictions)+1]] <- prediction
      
    }
    
  }
  error_list <- c(error_list,all_RMSs)
  
}



############################################################################################################
## STORING THE RESULT Frame1
############################################################################################################

simulation_results <- data.frame(mRNA=character(0),run=numeric(0),number_of_ribosomes=numeric(0),initiation_times=character(0),error=numeric(0),
                                 stringsAsFactors=F)
for(this_prediction in all_predictions){
  simulation_results[nrow(simulation_results)+1,] <- c(this_prediction$mRNA,
                                                       this_prediction$run,
                                                       this_prediction$number_of_ribosomes,
                                                       paste(round(this_prediction$initiation_times,2),collapse=","),
                                                       this_prediction$error
  )
}
############################################################################################################
## LOADING AND ADJUSTING Frame2_Data
############################################################################################################
fitwindow_Frame2=data.frame(seq(from=-round(Buildup_time_Frame2/acquisition_interval)+1,to=0,by=1),matrix(0,nrow=round(Buildup_time_Frame2/acquisition_interval),ncol=ncol(Frame2_Data)-1))

colnames(fitwindow_Frame2)<-colnames(Frame2_Data)
Frame2_Data=rbind(fitwindow_Frame2,Frame2_Data)
total_minutes_Frame2 <- nrow(Frame2_Data)/(60/acquisition_interval)    

timepoints_Frame2 <- seq(1,total_minutes_Frame2*60)
observed_timepoints <- timepoints[seq(acquisition_interval,length(timepoints_Frame2),acquisition_interval)] 
noise_level <- 1 # determines the acceptance probability

## The expected signal for a single ribosome

theoretical_signal_Frame2 <- function(initiation_time){ # initiation_time in seconds
  signal <- c(rep(0,initiation_time),    # before translation
              seq(0,1,1/Buildup_time_Frame2),    # translating the Frame1 array
              rep(1,Plateau_time_Frame2),    # translating the remainder of the transcript
              rep(0,length(timepoints_Frame2))  # after translation
  )
  return(single_ribosome_intensity_Frame2*signal[1:length(timepoints_Frame2)])
}
total_intensity_per_ribosome <- sum(theoretical_signal_Frame2(0))

last_timepoint_Frame2=integer()

for(i in 1:length(mRNAnames)){
  if(length(which(is.na(Frame2_Data[,mRNAnames[i]])==TRUE))==0){
    second.value=length(Frame2_Data[,mRNAnames[i]])+1
  }
  else if(length(split(which(is.na(Frame2_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame2_Data[,mRNAnames[i]]))) != 1))))==1 && max(split(which(is.na(Frame2_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame2_Data[,mRNAnames[i]]))) != 1)))$`1`)<length(Frame2_Data[,mRNAnames[i]])){
    second.value=length(Frame2_Data[,mRNAnames[i]])
  }
  else if(length(split(which(is.na(Frame2_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame2_Data[,mRNAnames[i]]))) != 1))))==1 && max(split(which(is.na(Frame2_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame2_Data[,mRNAnames[i]]))) != 1)))$`1`)==length(Frame2_Data[,mRNAnames[i]])){
    second.value=split(which(is.na(Frame2_Data[,mRNAnames[i]])), cumsum(c(1,diff(which(is.na(Frame2_Data[,mRNAnames[i]])))!=1)))
    second.value=second.value$`1`[1]
  }
  else{second.value=split(which(is.na(Frame2_Data[,mRNAnames[i]])), cumsum(c(1, diff(which(is.na(Frame2_Data[,mRNAnames[i]]))) != 1)))
  second.value=second.value$`2`[1]}
  last_timepoint_Frame2=c(last_timepoint_Frame2,second.value)
}

############################################################################################################
# ACTUAL SIMULATION
############################################################################################################

all_predictions_Frame2 <- list()
error_list_Frame2 <- c()

for(k in 1:length(mRNAnames)){  # select some mRNAs as an example
  observed_timepoints <- timepoints_Frame2[seq(acquisition_interval,(last_timepoint_Frame2[k]-1)*acquisition_interval,acquisition_interval)] # the timepoints at which we have an observation
  maxtime2=(last_timepoint_Frame2[k]-1)*acquisition_interval
  ## Start with the fluorescence signal
  real_intensity <- Frame2_Data[,mRNAnames[k]] 
  real_intensity[is.na(real_intensity)] <- 0
  real_intensity[real_intensity<0]<-0
  real_intensity[c(1:last_timepoint_Frame2[k]-1)]
  
  
  real_intensity_copy <-real_intensity
  real_intensity_copy <- rep(real_intensity_copy,each=acquisition_interval)
  
  
  real_intensity <- rep(real_intensity,each=acquisition_interval)
  timepoints=1:length(real_intensity)
  ## Estimate the number of ribosomes from the total intensity   
  estimated_number_of_ribosomes <- max(c(0,round(sum(real_intensity_copy)/total_intensity_per_ribosome)))
  
  if(estimated_number_of_ribosomes==0){ prediction <- list(mRNA = k,
                                                           run = 0,
                                                           number_of_ribosomes = 0,
                                                           initiation_times = 0,
                                                           error = allRMSs_nulfit)
  all_predictions_Frame2[[length(all_predictions_Frame2)+1]] <- prediction
  }  
  
  else{
    
    
    ## Randomly initialize the set of initiation times
    current_estimated_times_list <- replicate(parallel_runs,
                                              sort(sample(timepoints,size=estimated_number_of_ribosomes,
                                                          prob=c(real_intensity[(nrow(fitwindow_Frame2)*acquisition_interval+1):(2*nrow(fitwindow_Frame2)*acquisition_interval)]/sum(real_intensity[(17*acquisition_interval+1):length(timepoints)]),    real_intensity[(nrow(fitwindow_Frame2)*acquisition_interval+1):length(timepoints)]/sum(real_intensity[(17*acquisition_interval+1):length(timepoints)])))),
                                              simplify=F)
    
    ## Run the algorithm 
    for(i in 0:niterations){
      current_estimated_times_list <- lapply(current_estimated_times_list,function(current_estimated_times){
        # Slightly alter initiation times
        new_estimated_times <- sort(update(current_estimated_times,maxtime2,change_probability, change_size, add_probability, remove_probability))  
        # Calculate theoretical intensities and errors for both the old and the new initiation times
        current_theoretical_intensity <- rowSums(sapply(current_estimated_times,theoretical_signal_Frame2))
        new_theoretical_intensity <- rowSums(sapply(new_estimated_times,theoretical_signal_Frame2))
        current_error <- calculate_error(real_intensity[observed_timepoints[nrow(fitwindow_Frame2):length(observed_timepoints)]],current_theoretical_intensity[observed_timepoints[nrow(fitwindow_Frame2):length(observed_timepoints)]])
        new_error <- calculate_error(real_intensity[observed_timepoints[nrow(fitwindow_Frame2):length(observed_timepoints)]],new_theoretical_intensity[observed_timepoints[nrow(fitwindow_Frame2):length(observed_timepoints)]])
        
        # Calculate the acceptance probability 
        acceptance_probability <- exp((current_error - new_error)/(2*noise_level^2))
        
        # With probability acceptance_probability, keep the new initiation times 
        if(runif(1) < acceptance_probability){
          return(new_estimated_times)
        }
        else{
          return(current_estimated_times)
        }
      })
      
      if(i %% 100 == 0){
        print(paste(c("mRNA ",k,"; iteration ",i,"/",niterations),collapse=""))
      }
    }
    
    
    ## Calculate the RMS error for each of the found solutions
    all_RMSs <- sapply(1:parallel_runs,
                       function(x){calculate_RMS(rowSums(sapply(current_estimated_times_list[[x]],theoretical_signal_Frame2))[observed_timepoints[17:length(observed_timepoints)]],
                                                 real_intensity[observed_timepoints[17:length(observed_timepoints)]])})
    
    allRMSs_nulfit <- calculate_RMS(real_intensity[observed_timepoints],rep(0,length(observed_timepoints)))
    
    
    
    
    ## Store the result
    if(min(all_RMSs)<allRMSs_nulfit){
      for(x in 1:parallel_runs){
        
        
        
        prediction <- list(mRNA = k,
                           run = x,
                           number_of_ribosomes = length(current_estimated_times_list[[x]]),
                           initiation_times = current_estimated_times_list[[x]],
                           error = all_RMSs[x])
        all_predictions_Frame2[[length(all_predictions_Frame2) + 1]] <- prediction
      }
      
      
    }
    else{
      prediction <- list(mRNA = k,
                         run = 0,
                         number_of_ribosomes = 0,
                         initiation_times = 0,
                         error = allRMSs_nulfit)
      all_predictions_Frame2[[length(all_predictions_Frame2)+1]] <- prediction
      
    }
    
  }
  error_list_Frame2 <- c(error_list_Frame2,all_RMSs)
  
}

############################################################################################################
# STORING THE RESULT Frame2
############################################################################################################

simulation_results <- data.frame(mRNA=character(0),run=numeric(0),number_of_ribosomes=numeric(0),initiation_times=character(0),error=numeric(0),
                                 stringsAsFactors=F)
for(this_prediction in all_predictions_Frame2){
  simulation_results[nrow(simulation_results)+1,] <- c(this_prediction$mRNA,
                                                       this_prediction$run,
                                                       this_prediction$number_of_ribosomes,
                                                       paste(round(this_prediction$initiation_times,2),collapse=","),
                                                       this_prediction$error
  )
}
############################################################################################################
# PLOTTING THE Results
############################################################################################################
print(paste("Finished"))
select_save_location=choose.dir(default = "", caption = "Select folder")
number_of_mRNAs=as.numeric(unique(simulation_results$mRNA))
bestfit=data.frame()
pdf(paste(select_save_location,"Plots.pdf",sep="\\"), width=10,height=6)

for(k in 1:length(number_of_mRNAs)){
  z=number_of_mRNAs[k]
  #Load the predictions and the Frame1_Data for RNA Z
  observed_timepoints_Frame1 <- timepoints_Frame1[seq(acquisition_interval,(last_timepoint[z]-1)*acquisition_interval,acquisition_interval)] # the timepoints at which we have an observation
  real_intensity_Frame1 <- as.numeric(Frame1_Data[,mRNAnames[z]] )
  real_intensity_Frame1[is.na(real_intensity_Frame1)] <- 0
  real_intensity_Frame1=real_intensity_Frame1[c(1:last_timepoint[k]-1)]
  real_intensity_Frame1 <- rep(real_intensity_Frame1,each=acquisition_interval)
  predictions_Frame1 <- all_predictions[which(sapply(all_predictions,function(x){x$mRNA})==z)]
  
  observed_timepoints_Frame2 <- timepoints[seq(acquisition_interval,(last_timepoint_Frame2[z]-1)*acquisition_interval,acquisition_interval)] # the timepoints at which we have an observation
  real_intensity_Frame2 <- as.numeric(Frame2_Data[,mRNAnames[z]] )
  real_intensity_Frame2[is.na(real_intensity_Frame2)] <- 0
  real_intensity_Frame2[c(1:last_timepoint_Frame2[k]-1)]
  real_intensity_Frame2 <- rep(real_intensity_Frame2,each=acquisition_interval)
  predictions_Frame2 <- all_predictions_Frame2[which(sapply(all_predictions_Frame2,function(x){x$mRNA})==z)]
  
  #Order the predictions by error
  RMSorder <- order(sapply(predictions_Frame1,function(x){x$error}))
  
  
  #Plot the Frame1_Data and the fits
  ylimit=2500 
  if(max(real_intensity_Frame1[observed_timepoints_Frame1],na.rm=TRUE)>ylimit){
    ylimit=max(real_intensity_Frame1[observed_timepoints_Frame1])
  }
  if(max(real_intensity_Frame2[observed_timepoints_Frame2])>ylimit){
    ylimit=max(real_intensity_Frame2[observed_timepoints_Frame2])
  }
  
  if(nrow(fitwindow_Frame1)<nrow(fitwindow_Frame2)){
    Window_difference_Frame1=nrow(fitwindow_Frame2)-nrow(fitwindow_Frame1)
    Old_Time_Frame1=observed_timepoints_Frame1
    observed_timepoints_Frame1=acquisition_interval*Window_difference_Frame1+observed_timepoints_Frame1
    fitwindow=fitwindow_Frame2 
  }
  else{Old_Time_Frame1=observed_timepoints_Frame1
  Window_difference_Frame1=0}
  
  if(nrow(fitwindow_Frame2)<nrow(fitwindow_Frame1)){
    Window_difference_Frame2=nrow(fitwindow_Frame1)-nrow(fitwindow_Frame2)
    Old_Time_Frame2=observed_timepoints_Frame2
    observed_timepoints_Frame2=acquisition_interval*Window_difference_Frame2+observed_timepoints_Frame2
    fitwindow=fitwindow_Frame1 
  }
  else{Old_Time_Frame2=observed_timepoints_Frame2
  Window_difference_Frame2=0}
  
  if(nrow(fitwindow_Frame2)==nrow(fitwindow_Frame1)){
    Window_difference_Frame2=0
    Old_Time_Frame2=observed_timepoints_Frame2
    observed_timepoints_Frame2=acquisition_interval*Window_difference_Frame2+observed_timepoints_Frame2
    fitwindow=fitwindow_Frame1 
  }
  
  longer_signal=max(max(observed_timepoints_Frame2,na.rm=TRUE),max(observed_timepoints_Frame1,na.rm=TRUE))
  
  dt=.2*max(longer_signal)
  
  
  plot(observed_timepoints_Frame1,real_intensity_Frame1[Old_Time_Frame1],type="l",col="black",cex=0.5,
       main=paste("mRNA",z),xlab="Time (min)",ylab="Intensity (a.u.)",ylim=c(-2500,ylimit),xlim=c(0,max(longer_signal)+dt),xaxt="n",yaxt="n",
       cex.lab=1.5,cex.main=2)
  
  segments(x0=nrow(fitwindow)*acquisition_interval,y0=-500,x1=nrow(fitwindow)*acquisition_interval,y1=6000, lty=3)
  text(nrow(fitwindow)*acquisition_interval,ylimit-200,"T=0")
  
  abline(h=-500,lty=3)
  text(max(longer_signal)+dt*.25,-350,"#",col="gray")
  text(max(longer_signal)+dt*.75,-350,"RMSE",col="gray")
  
  
  
  if(length(predictions_Frame1)>1){
    a=rowSums(sapply(predictions_Frame1[[RMSorder[1]]]$initiation_times,theoretical_signal))[Old_Time_Frame1]
    #a[which(is.na(real_intensity_Frame1[observed_timepoints_Frame1]==TRUE))]=NA
    points(observed_timepoints_Frame1,a,
           type="l",lwd=2,col="green")  
    if(select_Frame2!=""){
      for(j in 1:10){
        i=RMSorder[j]                                   
        points(predictions_Frame1[[i]]$initiation_times+0*acquisition_interval+Window_difference_Frame1*acquisition_interval,rep(-400-200*j,predictions_Frame1[[i]]$number_of_ribosomes)+rep(c(-10,10),100)[1:predictions_Frame1[[i]]$number_of_ribosomes],
               bg="green",pch=24)
        abline(h=-500-200*j,lty=3)
        text(max(longer_signal)+dt*.125,-400-200*j,col="green",predictions_Frame1[[i]]$number_of_ribosomes)
        text(max(longer_signal)+dt*.625,-400-200*j,col="green",round(predictions_Frame1[[i]]$error,0))
        text(max(longer_signal)+dt*.25,-400-200*j,col="black","/")
        text(max(longer_signal)+dt*.75,-400-200*j,col="black","/")
      }
      
      
    }
  }
  if(length(predictions_Frame1)==1){points(observed_timepoints_Frame1,rep(0,length(observed_timepoints_Frame1)),type="l",lwd=2,col="green")
    for(j in 1:10){
      abline(h=-500-200*j,lty=3)
      text(max(longer_signal)+dt*.125,-400-200*j,col="green","0")
      text(max(longer_signal)+dt*.625,-400-200*j,col="green","0")
      text(max(longer_signal)+dt*.25,-400-200*j,col="black","/")
      text(max(longer_signal)+dt*.75,-400-200*j,col="black","/")
    }
  }
  
  #Order the predictions by error
  RMSorderm <- order(sapply(predictions_Frame2,function(x){x$error}))
  
  
  
  
  axis(side=2,at=seq(0,ylimit,1000),cex.axis=1.5)
  axis(side=1,at=seq(0,max(max(observed_timepoints_Frame1,na.rm=TRUE),max(observed_timepoints_Frame2,na.rm=TRUE)),acquisition_interval),round(seq(from=(-(nrow(fitwindow)/frames_per_minute)),to=((max(max(observed_timepoints_Frame1,na.rm=TRUE),max(observed_timepoints_Frame2,na.rm=TRUE))/acquisition_interval)-nrow(fitwindow))/frames_per_minute,by=1/frames_per_minute),1))
  
  if(select_Frame2!=""){
    
    lines(observed_timepoints_Frame2,real_intensity_Frame2[Old_Time_Frame2],col="red",type="l",cex=0.5)
    
    if(length(predictions_Frame2)>1){
      a=rowSums(sapply(predictions_Frame2[[RMSorderm[1]]]$initiation_times,theoretical_signal_Frame2))[Old_Time_Frame2]
      #a[which(is.na(real_intensity_Frame2[observed_timepoints_Frame2]==TRUE))]=NA
      points(observed_timepoints_Frame2,a,
             type="l",lwd=2,col="blue")
      abline(h=-500,lty=3)
      for(j in 1:10){
        i=RMSorderm[j]
        points(predictions_Frame2[[i]]$initiation_times+acquisition_interval*Window_difference_Frame2,rep(-400-200*j,predictions_Frame2[[i]]$number_of_ribosomes)+rep(c(-10,10),100)[1:predictions_Frame2[[i]]$number_of_ribosomes],
               bg="blue",pch=24)
        abline(h=-500-200*j,lty=3)
        text(max(longer_signal)+dt*.375,-400-200*j,col="blue",predictions_Frame2[[i]]$number_of_ribosomes)
        text(max(longer_signal)+dt*.875,-400-200*j,col="blue",round(predictions_Frame2[[i]]$error,0))
      }
      
      
    }
    if(length(predictions_Frame2)==1){points(observed_timepoints_Frame2,rep(0,length(observed_timepoints_Frame2)),type="l",lwd=2,col="blue")
      for(j in 1:10){
        abline(h=-500-200*j,lty=3)
        text(max(longer_signal)+dt*.375,-400-200*j,col="blue","0")
        text(max(longer_signal)+dt*.875,-400-200*j,col="blue","0")
      }
    }  
    legend(max(longer_signal)+dt*0.25, ylimit-100, legend=c("Frame 1","Frame 1 Fit","Frame 2","Frame 2 Fit"),
           col=c("black","green","red","blue"), lty=1, cex=0.58)
  }
  if(select_Frame2==""){
    for(j in 1:10){
      i=RMSorder[j]                                   
      points(predictions_Frame1[[i]]$initiation_times+0*acquisition_interval+Window_difference_Frame1*acquisition_interval,rep(-400-200*j,predictions_Frame1[[i]]$number_of_ribosomes)+rep(c(-10,10),100)[1:predictions_Frame1[[i]]$number_of_ribosomes],
             bg="green",pch=24)
      abline(h=-500-200*j,lty=3)
      text(max(longer_signal)+dt*.25,-400-200*j,col="green",predictions_Frame1[[i]]$number_of_ribosomes)
      text(max(longer_signal)+dt*.75,-400-200*j,col="green",round(predictions_Frame1[[i]]$error,0))
    }
    legend(max(longer_signal)+dt*0.25, ylimit-100, legend=c("Frame 1","Frame 1 Fit"),
           col=c("black","green"), lty=1, cex=0.58)}
  
  
  ############################################################################################################
  # Create Output File
  ############################################################################################################
  predictions_Frame1[[RMSorder[1]]]$Frame2_number_of_ribosomes<-predictions_Frame2[[RMSorderm[1]]]$number_of_ribosomes
  predictions_Frame1[[RMSorder[1]]]$Frame2_initiation_times<-predictions_Frame2[[RMSorderm[1]]]$initiation_times
  predictions_Frame1[[RMSorder[1]]]$Frame2_error<-predictions_Frame2[[RMSorderm[1]]]$error
  
  nofrm=NA  #number of ribosomes in Frame2 
  nofrs=NA  #number of ribosomes in Frame1 
  tm=NA     
  tm2=NA
  
  timeofone_Frame1=round((Buildup_time_Frame1+Plateau_time_Frame1)/acquisition_interval)
  timeofone_Frame2=round((Buildup_time_Frame2+Plateau_time_Frame2)/acquisition_interval)
  
  total_minutes=max(total_minutes_Frame1,total_minutes_Frame2)
  
  if(predictions_Frame1[[RMSorder[1]]]$Frame2_number_of_ribosomes>0){
    nofrm=round(predictions_Frame1[[RMSorder[1]]]$Frame2_initiation_times/acquisition_interval)
    nofrm=na.omit(nofrm)
    tm=matrix( rep( 0, len=(total_minutes*frames_per_minute*length(nofrm))), nrow = length(nofrm))
    tm2=matrix( rep( 0, len=(total_minutes*frames_per_minute*length(nofrs))), nrow = length(nofrs))
    for(j in 1:length(nofrm)){
      if(nofrm[j]+timeofone_Frame2<total_minutes*frames_per_minute){
        tm[j,nofrm[j]:(nofrm[j]+timeofone_Frame2)]=1
      }
      else(tm[j,nofrm[j]:(total_minutes*frames_per_minute)]=1)
    }
  }
  if(predictions_Frame1[[RMSorder[1]]]$number_of_ribosomes[1]>0){
    nofrs=round(predictions_Frame1[[RMSorder[1]]]$initiation_times/acquisition_interval)
    nofrs=na.omit(nofrs)
    tm2=matrix( rep( 0, len=(total_minutes*frames_per_minute*length(nofrs))), nrow = length(nofrs))
    for(m in 1:length(nofrs)){
      if(nofrs[m]+timeofone_Frame1<total_minutes*frames_per_minute){
        tm2[m,nofrs[m]:(nofrs[m]+timeofone_Frame1)]=1
      }
      else(tm2[m,nofrs[m]:(total_minutes*frames_per_minute)]=1)
    }
  }
  if(!is.na(tm2[1])){
    nofrs=colSums(tm2)}
  if(!is.na(tm[1])){
    nofrm=colSums(tm)}
  
  predictions_Frame1[[RMSorder[1]]]$time<-seq(from=(-(nrow(fitwindow)/frames_per_minute)),to=(nrow(Frame2_Data)-nrow(fitwindow))/frames_per_minute,by=1/frames_per_minute)
  predictions_Frame1[[RMSorder[1]]]$Frame1_number_of_ribosomes_per_frame<-nofrs
  predictions_Frame1[[RMSorder[1]]]$Frame2_number_of_ribosomes_per_frame<-nofrm
  
  datapoint=as.data.frame(list.as.matrix(predictions_Frame1[[RMSorder[1]]]))
  datapoint[,ncol(datapoint)+1]<-c(diff(datapoint$initiation_times)[1:length(datapoint$initiation_times)])
  datapoint[,ncol(datapoint)+1]<-c(diff(datapoint$Frame2_initiation_times)[1:length(datapoint$Frame2_initiation_times)])
  
  
  
  
  colnames(datapoint)[ncol(datapoint)]<-c("Frame2_Inititiation_Timing_Difference")
  colnames(datapoint)[ncol(datapoint)-1]<-c("Frame1_Inititiation_Timing_Difference")
  datapoint<-datapoint[,-c(2)]
  colnames(datapoint)[2:4]<-c("Frame1_number_of_ribosomes","Frame1_initiation_times","Frame1_error")
  datapoint$mRNA<-rep(datapoint$mRNA[1],nrow(datapoint))
  
  bestfit=rbind(bestfit,datapoint)
  
  
  
  
}
write.csv(bestfit,paste(select_save_location,"All_results.csv",sep="\\"))
dev.off()
