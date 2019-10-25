#' @title hbbr.Fit (Fits processed response data to hbbr model)
#'
#' @description Fits processed benefit-risk survey data from an appropriately
#'     designed discrete choice experiment to the hbbr (Hierarchical Bayesian
#'     Benefit-Risk) model. For details see article by Mukhopadhyay, S.,
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
#' @param design design information of the experiment:
#'      design = list(b, r, bl, rl, blbls, rlbls) where, b is number of benefit attributes,
#'      r is number of risk attributes, bl and rl are vectors of integers of length b and r
#'      indicating number of levels in j-th benefit attribute and k-th risk attribute,
#'      respectively. blbls, rlbls consists of labels for benefit and risk attributes.
#'      When blbls is NULL, it uses "B1", "B2", ... and similarly for rlbls.
#'
#' @param tune.param a list of tuning hyper-parameters to be used;
#'      default tune.param=list(tau=0.01, eta=NULL). See Details below for more
#'      information.
#'
#' @param mcmc a list of mcmc parameters to be used in the Gibbs sampler to obtain
#'      posterior samples of the paramaters of interests; default:
#'      mcmc=list(burnin=5000, iter=100000, nc=2, thin=20). See Details below for
#'      more information.
#'      
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate output 
#'      summary which might be helpful to see convergence results.     
#'
#' @return returns a list of useful output of interest and input specifications:
#'     (bbar.mcmc, bbar.means, bbar.sds, summary, logL, design, model, brdata, other.inputs).
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
#'      included in X as the part-worth parameter (beta) for the 1st level of each attribute
#'      is assumed to be 0 without loss of generality. So, if there are b benefit attributes
#'      and r risk attributes, and then have bl_j and rl_k levels (j=1,...,b; k=1,...,r)
#'      then total number of columns brdta is Sum_over_j(bl_j-1) + Sum_over_k(rl_k-1).
#'      If there are B respondents, each responding to k choice-pairs, then brdta will
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
#' @details mcmc is a list of MCMC specification parameters:
#'      (a) burnin - contains the number of burn-in values to be generated,
#'      (b) iter - is the total number of iterations of each chain beyond burn-in,
#'      (c) nc - is the number of independent chains, and
#'      (d) thin = posterior samples to be saved for every 'thin' values of the MCMC
#'      samples in each of the 'nc' chains. For more details see R2jags package help files.
#'
#' @examples ## Sample calls: fits pilot response data included with the package
#' \donttest{
#'
#'   data(hbbrPilotResp)
#'   hbfit = hbbr.Fit(brdta=hbbrPilotResp$brdta, design=hbbrPilotResp$design,
#'                    mcmc=list(burnin=500, iter=10000, nc=2, thin=10))
#'   hb = hbfit$bbar.mcmc
#'   dgn = hbfit$design
#'   mns = hbfit$bbar.means
#'   sds = hbfit$bbar.sd # same as apply(hbfit$bbar.mcmc, 2, sd)
#'
#'   ## Plots of MCMC draws ---------------------------------------
#'   op=par(mfrow=c(1,2), mar = c(4,2,3,1),oma=c(.1,.1,2,.1))
#'   matplot(hb,type="l",xlab="Iterations",ylab="",
#'           main=paste("Average Part-Worths (beta-bars)"),
#'           cex.main=.8, cex.lab=0.8, axes=FALSE)
#'   axis(1, at=seq(0,dim(hb)[1],length.out = 6),
#'           labels= paste(seq(0,5,1)*dim(hb)[1]/5 *hbfit$other.inputs$thin, sep=""),
#'           cex.axis=0.8)
#'   axis(2, cex.axis=0.8,las=1)
#'
#'   plot(hbfit$logL, type="l",main="Log Likelihood", axes=FALSE,xlab="Iterations",ylab="",
#'       cex.main=.8,cex.lab=0.8)
#'   axis(1, at=seq(0,dim(hb)[1],length.out = 6),
#'       labels= paste(seq(0,5,1)*dim(hb)[1]/5 *hbfit$other.inputs$thin, sep=""),
#'        cex.axis=0.8)
#'   axis(2, cex.axis=0.8,las=1)
#'   title(outer=TRUE, main = paste("MCMC draws plotted at every ",
#'        hbfit$other.inputs$thin,"-th Iteration",sep=""),cex.main=.9)
#'   par(op)
#'   
#'   ## Plots for mean estimated part-worth utilities ------------------
#'   require(ggplot2)
#'   require(gridExtra)
#'
#'   b.mns = c()
#'   b.sds = c()
#'   b.atr = c()
#'   b.lvl = c()
#'   j.now=1
#'   for (j in 1:dgn$b) {
#'     b.mns = c(b.mns,0, mns[j.now:(j.now-1+dgn$bl[j]-1)])
#'     b.sds = c(b.sds,0, sds[j.now:(j.now-1+dgn$bl[j]-1)])
#'     b.atr = c(b.atr, rep(dgn$blbls[j], dgn$bl[j]))
#'     b.lvl = c(b.lvl, paste("E", 1:dgn$bl[j],sep=""))
#'     j.now = j.now-1+dgn$bl[j]
#'   }
#'
#'   r.mns = c()
#'   r.sds = c()
#'   r.atr = c()
#'   r.lvl = c()
#'   k.now=j.now
#'   for (k in 1:dgn$r) {
#'     r.mns = c(r.mns,0,mns[k.now:(k.now-1+dgn$rl[k]-1)])
#'     r.sds = c(r.sds,0, sds[k.now:(k.now-1+dgn$rl[k]-1)])
#'     r.atr = c(r.atr, rep(dgn$rlbls[k], dgn$rl[k]))
#'     r.lvl = c(r.lvl, paste("H", 1:dgn$rl[k],sep=""))
#'     k.now = k.now-1+dgn$rl[k]
#'   }
#'
#'   d0.b = data.frame(Attributes =b.atr, lvl=b.lvl, util = b.mns, se = b.sds)
#'   d0.r = data.frame(Attributes =r.atr, lvl=r.lvl, util = r.mns, se = r.sds)
#'   y.max = max(abs(mns) + max(sds))
#'   pd <- position_dodge(0.2) # move them .2 to the left and right
#'
#'   pb = ggplot(data = d0.b, aes(x=lvl, y=util, group=Attributes,color=Attributes)) +
#'     ylim(0, y.max) +
#'     geom_hline(yintercept = 0) +
#'     geom_line(size=1.5, position=pd) +
#'     geom_point(size=4, shape=22, fill="green",color="darkgreen", position=pd) +
#'     geom_errorbar(aes(ymin=util-se, ymax=util+se), width=0.2, position=pd) +
#'     xlab("Benefit-Attribute Levels") + ylab("Estimated Utility") +
#'     ggtitle("Estimated Partworth Utilities of Benefits") +
#'     scale_color_manual(values=c("deepskyblue3" , "#9999CC", "cyan3" )) +
#'     theme(legend.position="bottom",plot.title = element_text(size = 10))
#'
#'   pr = ggplot( data = d0.r, aes(x=lvl, y=util, group=Attributes,color=Attributes)) +
#'     ylim(-y.max,0)+
#'     geom_hline(yintercept = 0) +
#'     geom_line(size=1.5, position=pd) +
#'     geom_point(size=4, shape=22, fill="pink",color="darkred", position=pd) +
#'     geom_errorbar(aes(ymin=util-se, ymax=util+se), width=0.2, position=pd) +
#'     xlab("Risk-Attribute Levels") + ylab("Estimated Utility") +
#'     ggtitle("Estimated Partworth Utilities of Risks") +
#'     scale_color_manual(values=c("orange" , "maroon" )) +
#'     theme(legend.position="bottom",plot.title = element_text(size = 10))
#'
#'   grid.arrange(pb, pr, nrow = 1)
#'
#' ##------------------------------------------------------------------
#' }
#'
#' @import R2jags
#' @export
hbbr.Fit <- function(brdta, design,
                      tune.param=list(tau=0.01, eta=NULL, df.add=2),
                      mcmc=list(burnin=5000, iter=100000, nc=2, thin=20), verbose=TRUE) {
  if (verbose) cat("Hello from hbbrFit!  \n")
  #--------- preparing the data for the model fitting -------
  subj=levels(factor(brdta[,1]))   #
  B=length(subj)                   # no of respondents
  k = dim(brdta)[1]/B              # no of choice pairs per respondent
  m = dim(brdta)[2]-2              # no of components in beta

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

  data.hbbr = list(y=y, X = X, B=B, k=k, m=m, df=df, zero=zero, S=S, Omega=Omega)

  #--- Fitting the DCE data with the basic HB model using JAGS --------------

  #--- The Model -------------------
  # B = number of respondent
  # k = number of choice pair for each respondent
  # m = number of part worth (excluding the 1st levels) (m <100)

  hbbr.model = function() {
    logL <- sum(log(y*p + (1-y)*(1-p)))
    for (h in 1:B){
      for (i in 1:k) {
        y[k*(h-1)+ i] ~ dbern(p[k*(h-1)+ i])
        p[k*(h-1)+ i] <- ilogit(X[k*(h-1)+ i,] %*% beta[h,])
      }
      beta[h,1:m] ~ dmnorm(bbar[], V[,]) # note V is precision matrix
    }
    bbar[1:m]~ dmnorm(zero[], S[,])
    V[1:m, 1:m] ~ dwish(Omega[,], df) #note E[V] = Inv(Omega)*(df)
  }

  parms = c('bbar') # note that for hbbr - we do not really need individual beta's
  set.seed(1234)  #jags.seed only works with jags.parallel - so using this for reproducible MCMCs

  start = proc.time()[3]

  jags.out = jags(data.hbbr,  parameters.to.save=parms,
                  model.file=hbbr.model, n.chains=mcmc$nc, n.iter=mcmc$iter, n.burnin=mcmc$burnin,
                  n.thin=mcmc$thin, jags.seed = 123, digits=3,
                  refresh = mcmc$iter/50, progress.bar = "text" )
  end = proc.time()[3]
  if (verbose) {
    cat(" Total Time Elapsed: ", round((end - start), 0), "Seconds", fill = TRUE)
  #fit.smry = data.frame(jags.out$BUGSoutput$summary)
    cat("\n\n|**************************************************|\n")
    cat("             summary of hbbr output                 \n")
    cat("|**************************************************|\n")
    print(round(jags.out$BUGSoutput$summary, 5), digits=3)
  }

  para= data.frame(jags.out$BUGSoutput$sims.list[[1]])
  out=list(bbar.mcmc=jags.out$BUGSoutput$sims.list$bbar,
           bbar.means = jags.out$BUGSoutput$mean$bbar,
           bbar.sds = jags.out$BUGSoutput$sd$bbar,
           summary = jags.out$BUGSoutput$summary,
           logL=-jags.out$BUGSoutput$sims.list$deviance/2, # note deviance = -2*logL
           design = design,
           model = jags.out$model,
           brdata = brdta,
           other.inputs=data.frame(nc=mcmc$nc, thin=mcmc$thin,
                                   iter=mcmc$iter, tau=tune.param$tau,
                                   eta=tune.param$eta, df=df)
           )
  out
}


