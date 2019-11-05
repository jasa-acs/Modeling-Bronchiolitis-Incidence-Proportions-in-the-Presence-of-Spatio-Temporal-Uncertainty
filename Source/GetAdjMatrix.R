##
## Code to calculate an (dim1 by dim2) adjacency matrix
##

AdjMat <- function(dim1,dim2){
  ## Who is my neighbor
  unit.num <- matrix(1:(dim1*dim2),nrow=dim1,byrow=TRUE)
  my.neigh <- vector("list",dim1*dim2)
  for(unit in 1:(dim1*dim2)){
    my.rc <- which(unit.num==unit,arr.ind=TRUE)
    ## Bottom, left, top, right
    if(my.rc[1]==1){
      if(my.rc[2]==1){
        my.neigh[[unit]] <- c(unit.num[my.rc[1]+1,my.rc[2]], #bottom
                              unit.num[my.rc[1],my.rc[2]+1]) #right
      } else if(my.rc[2]==dim2){
        my.neigh[[unit]] <- c(unit.num[my.rc[1]+1,my.rc[2]], #bottom
                              unit.num[my.rc[1],my.rc[2]-1]) #left
      } else {
        my.neigh[[unit]] <- c(unit.num[my.rc[1]+1,my.rc[2]], #bottom
                              unit.num[my.rc[1],my.rc[2]-1], #left
                              unit.num[my.rc[1],my.rc[2]+1]) #right
      }
    } else if(my.rc[1]==dim1){
      if(my.rc[2]==1){
        my.neigh[[unit]] <- c(unit.num[my.rc[1]-1,my.rc[2]], #top
                              unit.num[my.rc[1],my.rc[2]+1]) #right
      } else if(my.rc[2]==dim2){
        my.neigh[[unit]] <- c(unit.num[my.rc[1],my.rc[2]-1], #left
                              unit.num[my.rc[1]-1,my.rc[2]]) #top
      } else {
        my.neigh[[unit]] <- c(unit.num[my.rc[1],my.rc[2]-1], #left
                              unit.num[my.rc[1]-1,my.rc[2]], #top
                              unit.num[my.rc[1],my.rc[2]+1]) #right
      }
    } else {
      if(my.rc[2]==1){
        my.neigh[[unit]] <- c(unit.num[my.rc[1]+1,my.rc[2]], #bottom
                              unit.num[my.rc[1]-1,my.rc[2]], #top
                              unit.num[my.rc[1],my.rc[2]+1]) #right
      } else if(my.rc[2]==dim2){
        my.neigh[[unit]] <- c(unit.num[my.rc[1]+1,my.rc[2]], #bottom
                              unit.num[my.rc[1],my.rc[2]-1], #left
                              unit.num[my.rc[1]-1,my.rc[2]]) #top
      } else {
        my.neigh[[unit]] <- c(unit.num[my.rc[1]+1,my.rc[2]], #bottom
                              unit.num[my.rc[1],my.rc[2]-1], #left
                              unit.num[my.rc[1]-1,my.rc[2]], #top
                              unit.num[my.rc[1],my.rc[2]+1]) #right
      }
    }
    my.neigh[[unit]] <- sort(my.neigh[[unit]])
  }
  
  A <- matrix(0,nrow=dim1*dim2,ncol=dim1*dim2)
  for(unit in 1:(dim1*dim2)){
    A[unit,my.neigh[[unit]]] <- 1
  }
  
  return(list(A=A,unitlabel=unit.num))
  
}

AdjMat(5,5)

  
      