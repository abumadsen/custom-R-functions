#Mads F. Schou
#Functions for data formating

#1. Merge variable names
#2. Extract elements of a certain number of observations

#######################################
###--- 1. Merge variable names
#######################################

#merge n columns, c(camp,yr) -> camp_yr, except if any of the columns has an NA, then c(camp,yr) -> NA
MergeNoNA = function(dfToMerge){
  dfToMerge_char <- data.frame(lapply(dfToMerge, as.character)) #Change to character avoid spaces in the pasting below
  apply(dfToMerge_char,1,function(x){ifelse(any(is.na(x)), NA, paste(x,collapse="_"))})  #Paste by _
}

#######################################
###--- 2. Extract elements of a certain number of observations
#######################################

count.subset = function(vector_to_count = NULL, min_obs = NULL, max_obs = NULL){
  
  x_count = as.data.frame(table(vector_to_count)) #Count number of obs
  default_min = 0 #if no min is given
  default_max = max(x_count$Freq) #if no max is given
  
  print( ifelse( is.null(min_obs), paste('min_obs not specified, default =', default_min), paste('min_obs =',min_obs) ) )
  print( ifelse( is.null(max_obs), 'max_obs not specified: no max limit', paste('max_obs =',max_obs) ) )
  
  min_obs = ifelse(is.null(min_obs), default_min, min_obs) #if no min is given
  max_obs = ifelse(is.null(max_obs), default_max, max_obs) #if no max is given
  
  x_ok = x_count$vector_to_count[x_count$Freq >= min_obs & x_count$Freq <= max_obs] #filter
  return(as.character(x_ok))
  
}
