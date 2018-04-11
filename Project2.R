#Call data
data <- read.table("C://Onedrive/OneDrive - Knights - University of Central Florida/UCF/Courses/Statistical Computing/Project1/pb2.txt")

#Convert the data to matrix form
Data <- as.matrix(data, ncol=5)
N <- length(Y)
for (i in 1:N){
  if (Data[i,1] > 1){
    Data[i,1] <--1
  }
}


## Defining the Gaussian kernel
rbf_kernel <- function(x1,x2,gamma){
  K<-exp(-(1/gamma^2)*t(x1-x2)%*%(x1-x2))
  return(K)
}

## Conjugate Gradient Method to solving the problem Ax=B
conjugate_gradient_method<-function (A,B,N){
  M<-length(B)
  i<-1
  x<-matrix(0,N,M)
  r<-matrix(0,N,M)
  p<-matrix(0,N,M)
  beta<-numeric(N)
  lamda<-numeric(N)
  r[1,]<-B
  if (t(r[i,])%*%r[i,]>1e-5){
    i<-i+1
    if (i==2){
      p[i,]<-r[1,]
    }
    else {
      beta[i]<-t(r[i-1,])%*%r[i-1,]/t(r[i-2,])%*%r[i-2,]
      p[i,]<-r[i-1,]+beta[i]*p[i-1,]
    }
    lamda[i]<-t(r[i-1,])%*%r[i-1,]/t(p[i,])%*%A%*%p[i,]
    x[i,]<-x[i-1,]+lamda[i]*p[i,]
    r[i,]<-r[i-1,]-(lamda[i]*A)%*%p[i,]
  }
  return(x[i,])
}


## Training the LSSVM and predict the y for test dataset, then calculate the accuracy
lssvmmodel<-function(X, Y, Xt, Yt, step, gamma, C){
  N <- length(Y)
  Dm<-matrix(0,N,N)
  for(i in 1:N){
    for(j in 1:N){
      Dm[i,j]<-Y[i]*Y[j]*rbf_kernel(X[i,1:4],X[j,1:4],gamma)
    }
  }
  H<-Dm+diag(N)*(1/C)+diag(N)*1e-12 # adding a very small number to the diag, some trick
  d2<-as.vector(rep(1,N)) 
  # Call the Conjugate Gradient Method function to get the alpha and b
  nta<-conjugate_gradient_method(H,Y,step)
  vta<-conjugate_gradient_method(H,d2,step)
  s<-t(Y)%*%nta
  b<-(t(nta)%*%d2)/s
  alpha<-vta-nta*b
  alpha<-as.vector(alpha)
  
  # Predict the y for test dataset X_t
  M <- nrow(Xt)
  y_pred<-numeric(M)
  for(k in 1:M){
    ayK<-numeric(N)
    for (l in 1:N){
      ayK[l]<-alpha[l]*Y[l]*rbf_kernel(Xt[k,1:4],X[l,1:4],gamma)
    }
    y_pred[k] <- sign(sum(ayK)+b)
  }
  
  ### Evaluate the performance
  correct<-0
  for (m in 1:M){
    correct<-correct+(1-0.5*abs(Yt[m]-y_pred[m]))
  }
  accuracy<-correct/M
  return(accuracy)
}



################################
###  Problem 1-1   #############
################################



#Resample 60 rows from the raw dataset
Data<-Data[sample(1:62, 60, replace=FALSE), ]
#Choose 5-fold datasets split (5*12)- take the row number as the index of 5 fold datasets
folds <- cut(seq(1,nrow(Data)),breaks=5,labels=FALSE)


#Calculate the accuracy for different gamma and C value based 5-fold CV
optimal<-function(Data, C.range, gamma.range){
  N<-length(C.range)
  M<-length(gamma.range)
  tot.accuracy<-matrix(0,N,M)
  pb = txtProgressBar(min = 0, max = N, initial = 0) 
  for (i in 1:N){
    C<-C.range[i]
    setTxtProgressBar(pb,i)
    for (j in 1:M){
      gamma<-gamma.range[j]
      accuracy<-numeric(5)
      for (k in 1:5){
       #Segement the data by fold using the which() function 
       testIndexes <- which(folds==k,arr.ind=TRUE)
       testData <- Data[testIndexes, ]
       trainData <- Data[-testIndexes, ]
       #Use the test and train data partitions however you desire...
       train.X<-trainData[,2:5]
       train.Y<-trainData[,1]
       test.X<-testData[,2:5]
       test.Y<-testData[,1]
       accuracy[k]<-lssvmmodel(train.X, train.Y, test.X, test.Y, 500000, gamma, C)
     }
      tot.accuracy[i,j]<-sum(accuracy)/5
    }
  }
  return(tot.accuracy)
}

#Define the range of C value and gamma value
C.range<-seq(1,10, by=1)
gamma.range<-seq(0.1,4, by=0.1)


#Call the above function to generate the accuracy according to different C and gamma
test<-optimal(Data, C.range, gamma.range)

#index the maximum value in the accuracy matrix
valueindex<-which(test == max(test), arr.ind = TRUE)
valueindex

#Pick up the corresponding C value
C.list<-C.range[valueindex[,1]]
C.list

#Pick up the corresponding gamma value
gamma.list<-gamma.range[valueindex[,2]]
gamma.list


################################
###  Problem 1-2   #############
################################
## Training the LSSVM and predict the y for test dataset, then calculate the accuracy
lssvmmodel_LOOCV<-function(X, Y, Xt, Yt, step, gamma, C){
  N <- length(Y)
  Dm<-matrix(0,N,N)
  for(i in 1:N){
    for(j in 1:N){
      Dm[i,j]<-Y[i]*Y[j]*rbf_kernel(X[i,1:4],X[j,1:4],gamma)
    }
  }
  H<-Dm+diag(N)*(1/C)+diag(N)*1e-12 # adding a very small number to the diag, some trick
  d2<-as.vector(rep(1,N)) 
  # Call the Conjugate Gradient Method function to get the alpha and b
  nta<-conjugate_gradient_method(H,Y,step)
  vta<-conjugate_gradient_method(H,d2,step)
  s<-t(Y)%*%nta
  b<-(t(nta)%*%d2)/s
  alpha<-vta-nta*b
  alpha<-as.vector(alpha)
  
  # Predict the y for test dataset X_t
    ayK<-numeric(N)
    for (l in 1:N){
      ayK[l]<-alpha[l]*Y[l]*rbf_kernel(Xt,X[l,1:4],gamma)
    }
    y_pred <- sign(sum(ayK)+b)
  
  ### Evaluate the performance
  correct<-1-0.5*abs(Yt-y_pred)
  return(correct)
}



#Choose Leave-one-out cross validation
folds <- cut(seq(1,nrow(Data)),breaks=60,labels=FALSE)


#Calculate the accuracy for different gamma and C value based 5-fold CV
optimal_LOOV<-function(Data, C.range, gamma.range){
  N<-length(C.range)
  M<-length(gamma.range)
  tot.accuracy<-matrix(0,N,M)
  pb = txtProgressBar(min = 0, max = N, initial = 0) 
  for (i in 1:N){
    C<-C.range[i]
    setTxtProgressBar(pb,i)
    for (j in 1:M){
      gamma<-gamma.range[j]
      correct<-numeric(60)
      for (k in 1:60){
        #Segement the data by fold using the which() function 
        testIndexes <- which(folds==k,arr.ind=TRUE)
        testData <- Data[testIndexes, ]
        trainData <- Data[-testIndexes, ]
        #Use the test and train data partitions however you desire...
        train.X<-trainData[,2:5]
        train.Y<-trainData[,1]
        test.X<-testData[2:5]
        test.Y<-testData[1]
        correct[k]<-lssvmmodel_LOOCV(train.X, train.Y, test.X, test.Y, 500000, gamma, C)
      }
      tot.accuracy[i,j]<-sum(correct)/60
    }
  }
  return(tot.accuracy)
}

#Define the range of C value and gamma value
C.range<-seq(1,10, by=1)
gamma.range<-seq(0.1,4, by=0.1)

#Call the above function to generate the accuracy according to different C and gamma
test_LOOV<-optimal_LOOV(Data, C.range, gamma.range)
test_LOOV


#index the maximum value in the accuracy matrix
valueindex_LOOCV<-which(test_LOOV == max(test_LOOV), arr.ind=T)
valueindex_LOOCV

#Pick up the corresponding C value
C.list<-C.range[valueindex[,1]]
C.list

#Pick up the corresponding gamma value
gamma.list<-gamma.range[valueindex[,2]]
gamma.list



################################
###  Problem 2-1   #############
################################
library("rTensor")

################################
###  (a)   #############
################################

##Generate 1000 numbers from normal distribution
Tensor <- new("Tensor",3L,c(5L,5L,40L),data=rnorm(1000,0,1))

##Generate a matrix U with the column number equals to the 1th model of Tensor
U <- matrix(rep(1,100),ncol=5)

##Calculate teh Mode-1 product
A1<-ttm(Tensor,U,m=1)
A1
##Calculate teh Mode-2 product
A2<-ttm(Tensor,U,m=2)
A2

################################
###  (b)   #############
################################

##Generate a vector u with the column number equals to the 1th model of Tensor
u <- matrix(rep(1,5),ncol=5,nrow=1)

B1<-ttm(Tensor,u,m=1)
B1

B2<-ttm(Tensor,u,m=2)
B2

################################
###  (c)   #############
################################


tensor_reduce<-function(Tensor){
  N<-dim(Tensor)
  B<-Tensor
  N<-as.vector(N)
  M<-length(N)
  for (k in 1:M){
    u<-matrix(rep(1,N[k]),ncol=N[k],nrow=1)
    B<-ttm(B,u,m=k)
  }
  return(B)
}

A<-rand_tensor(modes = c(3, 4, 5), drop = FALSE)
tensor_reduce(A)

################################
###  2-2   #############
################################









