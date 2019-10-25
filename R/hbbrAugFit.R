#' @title hbbrAug.Fit (Fits processed response data to the augmented hbbr  model)
#'
#' @description Fits processed benefit-risk survey data from an appropriately
#'     designed discrete choice experiment to the augmented hbbr  model that
#'     includes patients' baseline characteristics. For details see article by Mukhopadhyay, S.,
#'     Dilley, K., Oladipo, A., & Jokinen, J. (2019). Hierarchical Bayesian
#'     Benefit–Risk Modeling and Assessment Using Choice Based Conjoint.
#'     Statistics in Biopharmaceutical Research, 11(1), 52-60.
#'
#' @author Saurabh Mukhopadhyay
#'
#' @param brdta processed and coded survey response data to be fitted to the hbbr model.
#'      It is a data frame in which 1st two columns indicate subject id and
#'      subject response (y = 0 or 1), and remaining columns contain information
#'      on design matrix (X). See Details below for more information.
#'
#' @param Z matrix of observed baseline characteristics of the patients. If there are
#'      N patients responded to the survey and we have included g
#'      characteristics of each patient then Z is a matrix of (g+1) x N, with all
#'      elements of the first column equal to 1.
#'      Note that when g=0, the model reduces to regular hbbr model.
#'
#' @param design design information of the experiment:
#'      design = list(b, r, bl, rl, blbls, rlbls) where, b is number of benefit attributes,
#'      r is number of risk attributes, bl and rl are vectors of integers of length b, and r
#'      indicating number of levels in j-th benefit attribute and k-th risk attribute,
#'      respectively. blbls, rlbls consists of labels of benefit and risk attributes.
#'      When blbls is NULL, it uses "B1", "B2", ... and similarly for rlbls.
#'
#' @param tune.param a list of tuning hyper-parameters to be used;
#'      default tune.param=list(tau=0.01, eta=NULL). See Details below for more
#'      information.
#'
#' @param mcmc a list of mcmc parameters to be used in the Gibbs sampler to obtain
#'      posterior samples of the parameters of interests; default:
#'      mcmc=list(burnin=1000, iter=10000, nc=4, thin=10). See Details below for
#'      more information.
#'      
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate output 
#'      summary which might be helpful to see convergence results.     
#'
#' @return returns a list of useful output of interest and input specifications:
#'     (del.mcmc, del.means, del.sds, summary, logL, design, model, brdata, other.inputs).
#'
#' @details brdta is a processed and coded survey response data to be fitted to the
#'      hbbr model. It is a data frame in which 1st column contains ID of respondent,
#'      2nd column contains response (y = 0 or 1) - each value corresponds to each
#'      choice-pair card evaluated by the respondent: y =1 if the 1st choice of the
#'      pair was preferred; 0 otherwise, 3rd column onwards contain information on
#'      design matrix (X). Each row of X is a vector of indicator variables taking
#'      values 0, 1, or -1; a value of 0 is used to denote absence of an attribute
#'      level; a value of 1 or -1 is used to indicate presence of an attribute
#'      level in the 1st choice, or in the 2nd choice, respectively in the choice-pair
#'      presented to the respondent.
#'      Note that column corresponding to the 1st level for each attribute would not be
#'      included as the part-worth parameter (beta) for the 1st level of each attribute
#'      is assumed to be 0 without loss of generality. So, if there are b benefit attributes
#'      and r risk attributes, and then have bl_j and rl_k levels (j=1,...,b; k=1,...,r)
#'      then total number of columns brdta is Sum_over_j(bl_j-1) + Sum_over_k(rl_k-1).
#'      If there are B respondents each responding to k choice-pairs then brdta will
#'      have B*k rows.
#'
#' @details tune.param is a list of tuning hyper-parameters (tau, eta) for the hbbr model.
#'      Specifically, in the hbbr model beta.h ~ MVN(beta.bar, V.beta) where the hyper-prior
#'      of beta.bar is assumed to be MVN (beta0, B) with B = 1/tau*I;  and
#'      hyper-prior of V.beta is assumed to follow inverse Wishart IW(nue, V) with V = 1/eta*I.
#'      When eta is NULL then eta will take the default value of m+3 which is the DF
#'      for the Wishart distribution. If we think the respondents have very similar
#'      part-worth vectors, then use eta=1.
#'
#' @details mcmc is a list of MCMC specification parameters to be used for rjags package:
#'      (a) burnin - contains the number of burn-in values to be generated,
#'      (b) iter - is the total number of iterations of each chain beyond burn-in,
#'      (c) nc - is the number of independent chains, and
#'      (d) thin = posterior samples to be saved for every 'thin' values of the MCMC
#'      samples in each of the 'nc' chains. For more details see rjags package help files.
#' @examples ## Sample calls:
#' # fits simulated response data included with this package to augmented hbbr model
#' # and then plots the estimated part-worth utilities.
#'
#' \donttest{
#'  data("simAugData")
#'  hbA = hbbrAug.Fit(brdta= simAugData$brdtaAug, Z=simAugData$Z,
#'                design=simAugData$design,
#'                tune.param=list(tau=0.01, eta=NULL, df.add=2),
#'                mcmc=list(burnin=500, iter=10000, nc=2, thin=10))
#'
#'  # define an appropriate function to plot the part-worth values...
#'  partworth.plot = function(attr.lvl, beta.mns, nb=3, new=TRUE, pnt =15, cl=clrs)
#'  {
#'  #check dimension
#'    k = length(attr.lvl) # no of attributes
#'    bk = length(unlist(attr.lvl)) # no of levels acrosss attributes
#'    if (bk - k != length(beta.mns)) stop("error 1")
#'    mns = rep(0, length(unlist(attr.lvl)))
#'    cntr = 0
#'    for (j in 1:k)
#'    {
#'      for (i in 1:length(attr.lvl[[j]])){
#'        cntr = cntr +1
#'        if (i > 1) mns[cntr]= beta.mns[cntr-1-(j-1)]
#'      }
#'    }
#'    indx = list()
#'    j0=1
#'    for (j in 1:k) {
#'      j1 = (j0+length(attr.lvl[[j]])-1)
#'      indx[[j]]= j0:j1
#'      j0=j1+1
#'   }
#'    if (new) {
#'      plot(c(1,bk), c(floor(min(beta.mns)*1.2),ceiling(max(beta.mns)*1.2)),
#'                                         type="n", axes=FALSE, xlab="",ylab="")
#'      axis(2, at=0:ceiling(max(beta.mns)*1.2), las=1, cex.axis=.7)
#'      axis(4, at=floor(min(beta.mns)*1.2):0, las=1, cex.axis=.7)
#'    }
#'    vl=c()
#'    for (j in 1:k)
#'    {
#'      points(indx[[j]], mns[indx[[j]]], type="b", pch =pnt, col=cl[j])
#'      vl=c(vl, max(indx[[j]])+.5)
#'    }
#'    abline(v=vl,col="gray", h=0)
#'    box()
#'  }
#'
#'  # Plotting estimated betas (part-worth) for some selected baseline characteristics:
#'
#'  augattr.lvl = list(b1=paste("B1",1:3,sep=""),b2=paste("B2",1:3,sep=""),
#'           r1=paste("R1",1:3,sep=""),r2=paste("R2",1:3,sep=""))
#'  clrs = c("blue", "green4","orange4", "red3")
#'
#'  mns = hbA$del.means
#'  # est. part-worth values
#'  betmn1 = mns %*% matrix(c(1, 0, 1), ncol=1)  # at mean age with disease staus=1
#'  betmn2 = mns %*% matrix(c(1, 0, -1), ncol=1) # at mean age with disease staus=-1
#'  betmn3 = mns %*% matrix(c(1, 1, -1), ncol=1) # at age = mean+1*SD, disease staus=-1
#'
#'  partworth.plot(attr.lvl = augattr.lvl, beta.mns = betmn1)
#'  partworth.plot(attr.lvl = augattr.lvl, beta.mns = betmn2, new=FALSE, pnt=17)
#'  partworth.plot(attr.lvl = augattr.lvl, beta.mns = betmn3, new=FALSE, pnt=16)
#'
#'  # Plotting true betas at those baseline characteristics
#'  Del = simAugData$Del
#'  clrs = rep("darkgrey", 4)
#'  # true part-worth values
#'  bmn1 = Del %*% matrix(c(1, 0, 1), ncol=1)  # at mean age with disease staus=1
#'  bmn2 = Del %*% matrix(c(1, 0, -1), ncol=1) # at mean age with disease staus=-1
#'  bmn3 = Del %*% matrix(c(1, 1, -1), ncol=1) # at age = mean+1*SD, disease staus=-1
#'  partworth.plot(attr.lvl = augattr.lvl, beta.mns = bmn1)
#'  partworth.plot(attr.lvl = augattr.lvl, beta.mns = bmn2, new=FALSE, pnt=17)
#'  partworth.plot(attr.lvl = augattr.lvl, beta.mns = bmn3, new=FALSE, pnt=16)
#' }
#'
#' @import R2jags
#' @export
#'
hbbrAug.Fit <- function(brdta, Z, design,
                         tune.param=list(tau=0.01, eta=NULL, df.add=2),
                         mcmc=list(burnin=5000, iter=100000, nc=2, thin=20), verbose=TRUE) {
  if (verbose) cat("Hello from hbbrAugFit! \n")
  #--------- preparing the data for the model fitting -------
  subj=levels(factor(brdta[,1]))   #
  B=length(subj)                   # no of respondents
  k = dim(brdta)[1]/B              # no of choice pairs per respondent
  m = dim(brdta)[2]-2              # no of components in beta
  g = dim(Z)[2]-1                  # no of characteristics measured for each respondent

  if (!(g>0) | (dim(Z)[1]!=B)   ) {
    stop("Z matrix must have at least 2 columns and appropriate number of rows")
  }

  if (g>0) {
    #--- check that design parameters are consistent to brdta -------
    if ((sum(design$bl)-design$b+sum(design$rl)-design$r) != m)
      stop("column dimention of data does not match with design")
    if (is.null(design$blbls)) design$blbls = paste("B", 1:design$b,sep="")
    if (is.null(design$rlbls)) design$rlbls = paste("R", 1:design$r,sep="")

    X=c()
    y=c()                            # y will arrange all responses according to X matrics
    for (h in 1:B) {
      y=c(y, brdta[brdta[,1]==subj[h],2])
      X=rbind(X, as.matrix(brdta[brdta[,1]==subj[h],c(3:(m+2))]))
    }

    if (tune.param$df.add<2) stop ("df must be m+2 or more")
    df = m+tune.param$df.add
    if (is.null(tune.param$eta))  tune.param$eta = df # default is DF

    S = tune.param$tau*diag(m)
    Omega=tune.param$eta*diag(m)
    zero=rep(0,m)

    data.hbbr.aug = list(y=y, X = X, Z=Z, g=g, B=B, k=k, m=m, df=df,
                         zero=zero, S=S, Omega=Omega)


    #--- The Model -------------------
    # B = number of respondent
    # k = number of choice pair for each respondent
    # m = number of part worth (excluding the 1st levels) (m <100)
    # r = number of characteristics being measured for each patient at baseline

    hbbr.aug.model = function(){
      for (h in 1:B){
        for (i in 1:k) {
          y[k*(h-1)+ i] ~ dbern(p[k*(h-1)+ i])
          p[k*(h-1)+ i] <- ilogit(X[k*(h-1)+ i,] %*% beta[h,])
        }
        beta[h,1:m] <- Del[,] %*% Z[h,] + eps[h,]
        eps[h,1:m] ~ dmnorm(zero[], V[,])
      }
      for (i in 1:(g+1)){
        Del[1:m,i] ~ dmnorm(zero[], S[,])
      }
      V[1:m, 1:m] ~ dwish(Omega[,], df)
    }

    parms = c("Del")
    set.seed(1234)  #jags.seed only works with jags.parallel - so using this for reproducible MCMCs

    start = proc.time()[3]

    jags.out = jags(data.hbbr.aug,  parameters.to.save=parms,
                    model.file=hbbr.aug.model, n.chains=mcmc$nc, n.iter=mcmc$iter, n.burnin=mcmc$burnin,
                    n.thin=mcmc$thin, jags.seed = 123, digits=3,
                    refresh = mcmc$iter/50, progress.bar = "text" )
    end = proc.time()[3]
    if (verbose){
      cat(" Total Time Elapsed: ", round((end - start)/60, 2), "Minutes", fill = TRUE)
  
      #fit.smry = data.frame(jags.out$BUGSoutput$summary)
      cat("\n\n|**************************************************|\n")
      cat(    "|         summary of augmented hbbr output         |\n")
      cat(    "|**************************************************|\n")
      print(round(jags.out$BUGSoutput$summary, 4), digits=3)
    }

    out = jags.out
    out=list(del.mcmc=jags.out$BUGSoutput$sims.list$Del,
             del.means = jags.out$BUGSoutput$mean$Del,
             del.sds = jags.out$BUGSoutput$sd$Del,
             summary = jags.out$BUGSoutput$summary,
             logL=-jags.out$BUGSoutput$sims.list$deviance/2, # note deviance = -2*logL
             design = design,
             model = jags.out$model,
             brdata = brdta,
             other.inputs=data.frame(nc=mcmc$nc, thin=mcmc$thin,
                                     iter=mcmc$iter, tau=tune.param$tau,
                                     eta=tune.param$eta, df=df)
    )
  }
  out
}



