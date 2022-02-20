# Euclidean distance 
dis = function(x)
{
  x = as.matrix(x)
  y = apply(x*x,1,sum) %*% matrix(1.0,1,nrow(x))
  sqrt(abs(y + t(y) - 2 * x %*% t(x)))
}

#dendrogram
iorder = function(m)
{
  N = nrow(m) + 1
  iorder = rep(0,N)
  iorder[1] = m[N-1,1]
  iorder[2] = m[N-1,2]
  loc = 2
  for(i in seq(N-2,1))
  {
    for(j in seq(1,loc))
    {
      if(iorder[j] == i)
      {
        iorder[j] = m[i,1]
        if(j==loc)
        {
          loc = loc + 1
          iorder[loc] = m[i,2]
        } else
        {
          loc = loc + 1
          for(k in seq(loc, j+2)) iorder[k] = iorder[k-1]
          iorder[j+1] = m[i,2]
        }
      }
    }
  }
  -iorder
}

#  hierarchical clustering implementation
hc = function(d, method=c("single","complete","average","centroid"))
{
  if(!is.matrix(d)) d = as.matrix(d)
  # clustering types:
  method_fn = switch(match.arg(method),
                     single   = min,
                     complete = max,
                     average  = mean,
                     centroid  = median)
  N = nrow(d)
  diag(d)=Inf
  n = -(1:N)                       
  m = matrix(0,nrow=N-1, ncol=2)   
  h = rep(0,N-1)                   
  for(j in seq(1,N-1))
  {
    # smallest distance 
    h[j] = min(d)
    
    i = which(d - h[j] == 0, arr.ind=TRUE)
    # take 1st.
    i = i[1,,drop=FALSE]
    p = n[i]
   
    p = p[order(p)]
    m[j,] = p
  
    # into the current jth group:
    grp = c(i, which(n %in% n[i[1,n[i]>0]]))
    n[grp] = j
    
    r = apply(d[i,],2,method_fn)
    
    # the distance matrix:
    d[min(i),] = d[,min(i)] = r
    d[min(i),min(i)]        = Inf
    d[max(i),] = d[,max(i)] = Inf
  }
  # Return Structure
  
  structure(list(merge = m, height = h, order = iorder(m),
                 labels = rownames(d), method = method, 
                 call = match.call(), dist.method = "euclidean"), 
            class = "hclust")
}
#loading data
ncidata = read.table("ncidata.txt",header = TRUE)  
ncidata1=t(ncidata)  

#calling and plotting
h11 = hc(dis(ncidata1), method="single")
plot(h11)

h21 = hc(dis(ncidata1), method="complete")
plot(h21)

h31 = hc(dis(ncidata1), method="average")
plot(h31)

h41 = hc(dis(ncidata1), method="centroid")
plot(h41)
