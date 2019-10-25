test_that("running test for hbbrFit", {
  data(hbbrPilotResp)
  hbfit = hbbr.Fit(brdta=hbbrPilotResp$brdta, design=hbbrPilotResp$design,
                   mcmc=list(burnin=500, iter=10000, nc=2, thin=10))

  expect_that( length(hbfit$bbar.means), equals(13) ) 
  expect_lt( abs(round(hbfit$bbar.means[1],1)-1.8) , .2) 
})


test_that("running test for hbbrAugFit", {
  data("simAugData")
  hbA = hbbrAug.Fit(brdta= simAugData$brdtaAug, Z=simAugData$Z,
                  design=simAugData$design,
                  tune.param=list(tau=0.01, eta=NULL, df.add=2),
                  mcmc=list(burnin=500, iter=3000, nc=2, thin=4))
  
  expect_that( length(hbA$del.means), equals(24) ) 
  expect_lt( abs(round(hbA$del.means[1,1],1)-4.7) ,.2) 
})