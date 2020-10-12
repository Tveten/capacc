# utility function to coerce data to an array structure
to_array<-function(X)
{
  if(!is.data.frame(X))
  {
     X<-as.array(X)
  }
  else
  {
    X<-as.array(array(X))
  }
  dims<-dim(X)
  if(length(dims) == 1)
  {
    X<-array(X,c(dims,1))
  }
  if(length(dims) > 2)
  {
    stop("data in array structures with dimension > 2 not supported")
  }
  return(X)
}
