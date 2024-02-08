# Functions from Humbert et al. 2009
# Define REML log-likelihoods
# REML objective function "negloglike.reml" is negative of log-likelihood
# for second differences of the log-scale observations.  The REML objective
# function uses equations A18-A25 from Humbert et al. (2009).  The three
# function arguments are:  theta, vector of parameters (transformed to the
# real line), yt, vector of time series observations (scaled), and
# tt, vector of observation times.  Function performs the differencing.

negloglike.reml = function(theta, yt, tt){
  sigsq = exp(theta[1]);         #  Constrains ssq > 0.
  tausq = exp(theta[2]);         #  Constrains tsq > 0.
  q = length(yt) - 1;
  qp1 = q + 1;
  vx = matrix(0, qp1, qp1);
  for (ti in 1:q)
  {
    vx[(ti + 1):qp1,(ti + 1):qp1] = matrix(1, 1, (qp1 - ti))*tt[ti + 1];
  }
  Sigma.mat = sigsq*vx;
  Itausq = matrix(rep(0, (qp1*qp1)), nrow = q + 1, ncol = q + 1);
  diag(Itausq) = rep(tausq, q + 1);
  V = Sigma.mat + Itausq;
  ss = tt[2:qp1] - tt[1:q];
  D1mat = cbind(-diag(1/ss), matrix(0, q, 1)) + cbind(matrix(0, q, 1), diag(1/ss));
  D2mat = cbind(-diag(1, q - 1), matrix(0,q - 1,1)) +
    cbind(matrix(0, q - 1, 1), diag(1,q - 1));
  Phi.mat = D2mat%*%D1mat%*%V%*%t(D1mat)%*%t(D2mat);
  wt = (yt[2:qp1] - yt[1:q])/ss;
  ut = wt[2:q] - wt[1:q - 1];
  ofn = (q/2)*log(2*pi) + (0.5*log(det(Phi.mat)))+
    (0.5*(ut%*%ginv(Phi.mat)%*%ut));
  return(ofn);
}