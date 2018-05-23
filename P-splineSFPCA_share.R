#P-spline smoothed FPCA developed on the paper: PLEASE CITE THE FOLLOWING PAPER IF YOU USE THIS RESULT

#Aguilera, A. nad Aguilera-Morilo, M.C. (2013). 
#Penalized PCA approaches for B-spline expansions of smooth functional data.
#Applied Mathematics and Computation, 219, pp. 7805–7819.

#This code is property of the authors of the paper: 
#M. Carmen Aguilera-Morillo (Universidad Carlos III de Madrid) 
#Ana Aguilera (Universiodad de Granada)

#(PLEASE CITE THIS PAPER IF YOU USE THIS RESULT)

  #data: row data matrix
  #times: vector with the observation knots
  #nbasis: number of basis functions
  #lambda: smoothing parameter
  #n: number of functional observations
  m=length(times)
  library(fda)
  #Basis representation
  bspl <- create.bspline.basis(c(min(times), max(times)), nbasis=nbasis)
  Xs<- Data2fd(data, argvals=times, basisobj=bspl)
  Xs_centered = center.fd(Xs)
  cs <- t(coef(Xs_centered))  
  
  #P-spline Penalty matrix
  P = inprod(bspl, bspl, 2, 2) #silverman penalty
  p=nbasis # number of basis functions
  D=diff(diff(diag(p+2)))[1:(p-2),1:p] 
  P2 = t(D)%*%D
  
  #Penalization: 
  Penalty = P + lambda*P2
  L = mroot(Penalty, method="svd")
  
  #In Aguilera and Aguilera-Morillo (2013) 
  #was proved that P-spline smoothed FPCA is equivalent to
  #clasical PCA of the following matrix: 
  matrixPCA <- cs%*%t(ginv(L)%*%P)
  
  #Clasical PCA
  SPCA<-princomp(matrizPCA, scores=TRUE, cor=FALSE)
  
  #scores
  scores<-SPCA$scores
  
  #loadings(SPCA): in the paper they are denoted by u_j, j=1,...,ncomponents (columns of matrix u)
  u<- SPCA$loadings
  
  #From the loadings matrix, you can obtain the PC curves or weight functions as esplained in the paper:
    #Basis coefficients of the weight functions: in the paper they are denoted by b_j= (L^-1)' * u_j 
    #One time you have the Basis coefficients of the weight functions, you can obtain the weight functions 
      #as usually in FDA (using the function eval.basis). You can evaluate the B-spline basis on a larger 
      #times vector for a correct graphical representation.
#-------------------------------------------------------------------------------