#john, 2016-10-11
#this code is helpful for understanding the logic of data that was the product
#of a SQL pull written by an analyst who does not document their work, does
#not respond to e-mails, does not pick up the phone, and gets mad when
#i visit in person. in other words, this code is useful for understanding
#spreadsheets from Sheriff's BI.

#"check" a single column
#by "check", i mean: given a "cause", i.e., a set of columns that look like
#they might identify certain "layers" in the data (like bookingID and
#ccoms case number), see if they "explain" a given "effect", e.g.,
#see if bond amount is invariant for all rows corresponding to a certain
#bookingID*caseNum combination.
checkUnitsCol <-
  function(
    data,              #a data.table
    cause,             #names of unit columns
    effect             #name of a column to be explained
  )   
  {
    require(data.table)
    
    #decide how many unit combinations we have to test
    
    numUnits <- length(cause)
    numCombos <- 2^numUnits
    
    #make a table of combinations
    
    comboTable <- data.frame(row.names=1:numCombos)
    for(i in 1:numUnits) {
      column <- rep(c(0,1), each = 2^(i-1), length = 2^numUnits)
      comboTable <- cbind(comboTable,column)
      colnames(comboTable)[i] <- cause[i]
    }
    
    #make a hierarchy of units
    
    #we prefer to use simple units.
    #given a choice between equally complex units,
    #we choose the units in the order given by the user.
    
    comboTable$complexity <- apply(comboTable,1,sum)
    comboTable <- comboTable[order(comboTable$complexity),]
    
    #define a metric of "disorderliness" or "unpredictability".
    #this one is pretty crude: it just counts whether the number
    #of unique values in (some subset of) the data is bigger than 1.
    #for example, if there is more than one unique value for mental health
    #status under a single bookingid, then we would say that mental health
    #status is unpredictable and disorderly for that bookingid.
    
      disorderliness <- function(z) length(unique(z))>1
    
    #test the disorderliness of each combination
    
    comboTable$messiness <- as.numeric(NA)
    for(row in 1:numCombos){
      thingsToTest <- unlist(comboTable[row,1:numUnits])
      if(max(thingsToTest)==0) {
        unitOfTheHour <- ''
      } else {
        unitOfTheHour <- cause[as.logical(thingsToTest)]
          #list(data.dt[,c(as.logical(thingsToTest),FALSE)])
      }
      messinessByPerson <- data[,
                                lapply(.SD, disorderliness),
                                by=unitOfTheHour,
                                .SDcols=which(colnames(data)==effect)
                                ]
      averageMessiness <- sum(messinessByPerson[[ncol(messinessByPerson)]]) /
        nrow(messinessByPerson)
      comboTable$messiness[row] <- averageMessiness
    }
    
    #return the results of the tests
    
    return(comboTable)
    
  }

#analyze checkUnitsCol results using occam's razor

occamsDullRazor <-
 function(
  resultsTable 	#the results from checkUnitCol
 )
 {

  #apply a crude version of occam's razor:
  #among explanations that predict the data equally well,
  #pick the simplest one
  
   #this ordering puts the "best" explanations at the top

    resultsTable <- resultsTable[order(resultsTable$messiness,resultsTable$complexity),]

   #i say an explanation is "occam indeterminate" if there is another with equal
   #predictive power and equal simplicity

    resultsTable$occamIndeterm <- 
     ave(x=1:nrow(resultsTable),
      resultsTable$messiness,resultsTable$complexity,
      FUN = function(x) length(x)>1)
    resultsTable$occamIndeterm <- as.logical(resultsTable$occamIndeterm)

   #just return the single best explanation
   
    return(resultsTable[1,])
 }

#apply checkUnitsCol and occamsDullRazor to an entire dataset
checkUnitsDF <-
  function(
    data,		#a whole dataframe
    cause   #names of columns that we'll use to explain the dataframe
  )
  {
    
    #make everything a factor and then a number (for speed)
    data <- mapply(FUN = function(x) as.numeric(as.factor(x)),
                   data)
    
    #first, get a data.table of the data
    #(independent package; faster for this sort of thing)
    require(data.table)
    data <- data.table(data)
    #print('hi') #for debug
    
    
    #just run checkUnitsCol2 in a loop, as many times as we need,
    #and store the results
    occamTable <- data.frame()
    for(col in colnames(data)){
      #print(paste0('started ',col)) #for debug
      startTime <- proc.time()
      bestUnit <- occamsDullRazor(checkUnitsCol(data,cause,col))
      bestUnit <- cbind(col,bestUnit)
      occamTable <- rbind(occamTable,bestUnit)
      stopTime <- proc.time()
      elapsed <- (stopTime - startTime)[3]
      #print(paste0('finished ',col,' in ',elapsed,' sec')) #for debug
    }
    return(occamTable)
  }

#"smoosh" a variable to fit a certain set of units
#e.g., if we sometimes have contradictory information about bond within
#specific bookingID*caseNum combinations, resolve this contradiction
#by keeping the highest known bond, and then compress the dataset to have a
#single row for each specific bookingID*caseNum combination

# updating 2022-06-03

smoosh <- function(
 fromData,   #origin dataset; something in it needs smooshing
 fromColumn, #column in fromData to smoosh (given by column name)
 byColumns,  #units to smoosh into (given by a vector of column names)
 toData,     #destination dataset; units are pre-smooshed
		         # (rows of toData are 1:1 w/ unique values of byColumns)
 aggFun,     #how to smoosh (e.g. add up all the numbers,
             # paste together all the strings...)
 newName     #what to call the column in the new dataset
 ) {
 require(data.table)
 fromData = data.table(fromData)
 smooshed <- fromData[,
  lapply(.SD, aggFun),
  by = byColumns,
  .SDcols = which(colnames(fromData) == fromColumn)
  ]
 colnames(smooshed)[which(colnames(smooshed) == fromColumn)] <- newName
 merged <- merge(
  x = toData,
  y = smooshed,
  all.x = TRUE,
  by = byColumns,
  sort = FALSE #original order will remain
 )
 return(merged)
}
