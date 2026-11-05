buildMatrix = function(tableau, isMax){
  ncol=ncol(tableau)
  nrow=nrow(tableau)
  
  #store the last col in variable to be appended at the end later
  if(isMax){ #max
    ansCol = tableau[,ncol]
  }else{ #min/dual of max
    ansCol = tableau[nrow,]
    
    #transpose
    tableau=t(tableau)
    
    #get new cols and rows
    nrow=nrow(tableau)
    ncol=ncol(tableau)
  }
  
  #create identity
  im = matrix(c(0),nrow=nrow,ncol=nrow)
  for(i in 1:nrow){
    for(j in 1:nrow){
      if(i==j){
        im[i,j]=1
      }
    }
  }
  
  #add slack cols
  tableau = cbind(tableau[,1:(ncol-1)],im,ansCol)
  
  #negate bottom row
  tableau[nrow,1:(ncol+nrow-2)] = tableau[nrow,1:(ncol+nrow-2)]*-1
  
  return(tableau)
}

#catch:
# no group selected or co2 cap is 0 or less or blank

Simplex = function(tableau, isMax){
  #build the matrix
  tableau = buildMatrix(tableau,isMax)
  
  #do maximization
  ncol = ncol(tableau)
  nrow = nrow(tableau)
  iter=0
  
  #iterate until no more negative values in the bottom row
  while(any(tableau[nrow, 1:(ncol-1)] < 0)){
    iter = iter+1
    #stores the larger number to be the pivot col
    largest = 0
    #store the test ratio
    testRatio = 0
    #store the current smallest test ratio
    STR = 0
    
    #store pivot col and pivot row
    PC = 0
    PR = 0
    PE = 0
    
    #get pivot col
    for(i in 1:(ncol-2)){
      #check the bottom and get the largest absolute number that is also negative
      if(abs(tableau[nrow,i])>largest && tableau[nrow,i]<0){
        largest = abs(tableau[nrow,i])
        PC = i
      }
    }
    
    #get an initial test ratio, STR, and PR
    STR = Inf
    
    #get pivot row
    for(j in 1:nrow){
      #if non negative, get the test ratio then 
      if(tableau[j,PC]>0){
        testRatio = (tableau[j,ncol]/tableau[j,PC])
        if(testRatio<=STR){
          STR = testRatio
          PR = j
        }
      }
    }
    
    if(STR>=Inf){
      print("The problem is infeasible!")
      return(list(finalTableau=NULL,basicSolution=NULL,numberOfIterations=iter))
    }
    
    #get PE
    PE = tableau[PR,PC]
    
    #normalize pivot row
    tableau[PR,] = tableau[PR,]/PE
    
    #use GJ to clear the col
    for(j in 1:nrow){
      if(j!=PR){
        tableau[j,] = tableau[j,] - (tableau[PR,]*tableau[j,PC])
      }
    }
  }
  
  #initialize basic sol matrix
  basicSolution = matrix(c(0),nrow=1,ncol=ncol)
  dimnames(basicSolution) <- dimnames(tableau)
  nonzeroRow = -1
  
  if(isMax){
    #get max sol
    for(i in 1:(ncol-1)){
      nonzeroRow = -1
      for(j in 1:nrow){
        #check if there is only 1 non-zero digit
        if(tableau[j,i]!=0){
          if(nonzeroRow == -1){
            nonzeroRow = j
          }else{
            #if there are multiple nonzero rows, reset value of nzr and break
            nonzeroRow = -1
            break
          }
        }
      }
      
      #there were no nonzero numbers in the row or there were multiple nonzero rows
      if(nonzeroRow==-1){
        basicSolution[1,i] = 0
      }else{
        basicSolution[1,i] = tableau[nonzeroRow,ncol]
      }
    }
  }else{
    #get min sol
    for(i in 1:ncol){
      basicSolution[1,i] = tableau[nrow,i]
    }
  }
  
  #put in a list
  solution = list(finalTableau=tableau,basicSolution=basicSolution,numberOfIterations=iter)
  
  #return basic sol mat
  solution
}