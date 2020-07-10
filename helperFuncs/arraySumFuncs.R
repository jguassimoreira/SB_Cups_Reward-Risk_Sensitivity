#######################################
## 3D array sum for initial grid search
colSumArray3 <- function(inArray) {
  
  #Sum of the columns across all dimensions in a 3D array, i.e., collapse across
  #first, and third dimensions
  #Only works for 3D arrays
  
  #Get all four dimensions of the matrix
  d1 = dim(inArray)[1]
  d2 = dim(inArray)[2]
  d3 = dim(inArray)[3]
  
  sumMat = matrix(NaN, d3, d2) #Num columns remains the same, Num rows is  = 3dim
  rcount = 1
  
  #loop over 3rd and 4th dimensions, collapse across rows in each & sum, add to sumMat
  for (i in 1:d3) {
    sumMat[rcount,] = .colSums(inArray[,,i], d1, d2)
    rcount = rcount + 1
  }
  
  colSumFinal = .colSums(sumMat, d3, d2)
  
  return(colSumFinal)
}


rowSumArray3 <- function(inArray) {
  
  #Sum of the columns across all dimensions in a 3D array, i.e., collapse across
  #second, and third dimensions
  #Only works for 3D arrays
  
  #Get all four dimensions of the matrix
  d1 = dim(inArray)[1]
  d2 = dim(inArray)[2]
  d3 = dim(inArray)[3]
  
  sumMat = matrix(NaN, d1, d3) #Num rows remains the same, Num cols is = 3dim
  ccount = 1
  
  #loop over 3rd and 4th dimensions, collapse across rows in each & sum, add to sumMat
  for (i in 1:d3) {
    sumMat[,ccount] = .rowSums(inArray[,,i], d1, d2)
    ccount = ccount + 1
  }
  
  rowSumFinal = .rowSums(sumMat, d1, d3)
  rowSumFinal = t(t(rowSumFinal))
  
  return(rowSumFinal)
}


d3SumArray3 <- function(inArray) {
  
  #Sum of all matrices across third dimension in a 3D array, i.e., collapse
  #across first, and second dimensions
  #only works for 3D arrays
  
  #Get all four dimensions of the matrix
  d1 = dim(inArray)[1]
  d2 = dim(inArray)[2]
  d3 = dim(inArray)[3]
  
  sumMat = matrix(NaN, 1, d3)
  
  #Loop through each 'page' (3rd dim) of the array
  for (i in 1:d3) {
    sumMat[1,i] = sum(inArray[,,i])
  }
  
  d3SumFinal = sumMat
  
  return(d3SumFinal)
  
}
