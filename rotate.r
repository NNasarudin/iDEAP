
#This file is part of DEAP.

#DEAP is free software: you can redistribute it and/or modify it under 
#the terms of the GNU Lesser General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#DEAP is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Lesser General Public License for more details.

#You should have received a copy of the GNU Lesser General Public License
#along with DEAP.  If not, see <http://www.gnu.org/licenses/>.

#This file was written by Roger Higdon. If you have any questions, 
#please contact me at: roger.higdon@seattlechildrens.org


#fit.resid should be run first only once to correct data for overall mean and covariates
#rotate.adjusted.y should be run as many times as necessary to generate rotation samples
#rotate.y is called by rotate.adjusted.y


fit.resid=function(y,x=NULL){
#calculates regression fitted values, basis matrix for residuals, and subspace (length(y)-ncol(x) vector for residuals
if (is.vector(y)) y=matrix(y,ncol=1)
if (!is.matrix(y)) stop("y must be a vector or matrix")
if (length(x)==0) return(list(yfit=NULL,yr=y,W=NULL)) else {
 if (is.vector(x)) x=matrix(x,ncol=1)
 if (!is.matrix(x)) stop("x must be a vector or matrix")
 if (nrow(x)!=nrow(y)) stop ("x and y must have the same number of rows")
 p=ncol(x)
 n=nrow(x)
 if (n<=p) stop("x must have more rows than columns")
 hat=x%*%solve(t(x)%*%x)%*%t(x)
 id=diag(rep(1,n))
 yfit=hat%*%y
 W=svd(id-hat)$u[,1:(n-p)]
 yr=t(W)%*%y
 return(list(yfit=yfit,W=W,yr=yr))}
}

rotate.y=function(y){
#generates a random rotation of a  vector or columns of a matrix
 if (is.vector(y)) y=matrix(y,ncol=1)
 if (!is.matrix(y)) stop("y must be a vector or matrix")
 n=nrow(y)
 x=matrix(rnorm(n*n),ncol=n)
 qrx=qr(x)
 qr.qty(qrx,y)}

rotate.adjusted.y=function(fit){
#creates random rotation of vector or columns of a matrix adusted by explanatory variable
#input is the list from fit.resid
if (length(fit$yfit)==0) return(rotate.y(fit$yr)) else return(fit$yfit+fit$W%*%rotate.y(fit$yr))
}

