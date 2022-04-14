library(ape)
library(bayou)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (!is.numeric(as.numeric(args[1]))){
  stop("Argument must be a number", call.=FALSE)
}

#Read data files
all_data <- read.csv("/n/holyscratch01/informatics/sshakya/avonet/AVONET1_BirdLife.csv", header = T)
all_tree <- read.nexus("/n/holyscratch01/informatics/sshakya/avonet/HackettStage1_0001_1000_MCCTreeTargetHeights.nex")
output_folder <- "/n/holyscratch01/informatics/sshakya/avonet/runs/"
runs <- 5000000

#Some data cleanup
all_data$species <- gsub(" ", "_", all_data$Species1)
all_data_in_tree <- subset(all_data, species %in% all_tree$tip.label)
all_data_in_tree <- all_data_in_tree[,c("Group", "species", "Mass", "Tarsus.Length"),]
all_data_in_tree$ln_tarsus <- log(all_data_in_tree$Tarsus.Length)
all_data_in_tree$ln_mass <- log(all_data_in_tree$Mass)

#gather group subset
all_groups <- unique(all_data$Group)
group <- grep(paste("^",as.numeric(args[1]),"_",sep=""), all_groups, value = T)
group_init <- substr(group, 1, 8)

#assign output iteration number
if (dir.exists(output_folder)){
  already_files <- list.files(path = output_folder, pattern = group_init)
  iter <- round(length(already_files)/6,digits = 0) + 1  #each run creates 6 files including the start tree png
}

#prepare data for bayou input
all_data.1 <- subset(all_data_in_tree, Group %in% group)
all_tree.1 <- keep.tip(all_tree, all_data.1$species)
all_data.1$species <- factor(all_data.1$species, levels = all_tree.1$tip.label)
all_data.1 <- all_data.1[order(all_data.1$species),]
#note: order of tips and data should be same

dat <- all_data.1$ln_tarsus
names(dat) <- all_data.1$species

pred <- data.frame(lnMass = all_data.1$ln_mass)
rownames(pred) <- all_data.1$species

priors <- make.prior(all_tree.1, plot.prior = F,
                     dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                dsb="dsb", dk="cdpois", dtheta="dnorm", dloc="dunif"),
                     param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                dbeta_lnMass=list(mean=2, sd=0.5),
                                dk=list(lambda=if(nrow(all_data.1)*0.02 < 1){1}else{floor(nrow(all_data.1)*0.02)}, 
                                        kmax=floor(nrow(all_data.1)*0.2)),
                                dsb=list(bmax=1,prob=1), dtheta=list(mean=0.3, sd=0.5)))
#slope is theta and intercept is beta
#make sure lambda of dk is not 0 otherwise the distribution is centered at 0 only

DN1 <- list(list(alpha=4, sig2=3.6, beta_lnMass=0.1, k=c(1,1), theta=1, slide=1),
            list(alpha=7, sig2=1.3, beta_lnMass=0.2, k=c(1,1), theta=2.5, slide=1),
            list(alpha=6, sig2=2, beta_lnMass=0.2, k=c(1,1), theta=1, slide=1),
            list(alpha=9, sig2=1.2, beta_lnMass=0.2, k=c(1,1), theta=2, slide=1),
            list(alpha=7, sig2=1.3, beta_lnMass=0.4, k=c(1,1), theta=2.5, slide=1),
            list(alpha=6, sig2=2.3, beta_lnMass=0.2, k=c(1,1), theta=1, slide=1),
            list(alpha=8, sig2=1, beta_lnMass=0.3, k=c(1,1), theta=3, slide=1),
            list(alpha=8, sig2=1.2, beta_lnMass=0.2, k=c(1,1), theta=2, slide=1),
            list(alpha=9, sig2=1.2, beta_lnMass=0.2, k=c(1,1), theta=2.5, slide=1),
            list(alpha=6, sig2=2.6, beta_lnMass=0.3, k=c(1,1), theta=2, slide=1),
            list(alpha=7, sig2=1.8, beta_lnMass=0.2, k=c(1,1), theta=1.5, slide=1),
            list(alpha=5, sig2=1, beta_lnMass=0.5, k=c(1,1), theta=3.5, slide=1),
            list(alpha=7, sig2=1.8, beta_lnMass=0.1, k=c(1,1), theta=0.5, slide=1),
            list(alpha=7, sig2=1.2, beta_lnMass=0.3, k=c(1,1), theta=2.5, slide=1),
            list(alpha=7, sig2=1, beta_lnMass=0.3, k=c(1,1), theta=3, slide=1),
            list(alpha=6, sig2=1, beta_lnMass=0.3, k=c(1,1), theta=3, slide=1),
            list(alpha=7, sig2=3, beta_lnMass=0.2, k=c(1,1), theta=0.5, slide=1),
            list(alpha=6, sig2=1, beta_lnMass=0.3, k=c(1,1), theta=3, slide=1))
#run short runs and adjust values based on acceptance ratio (in output as .alpha for alpha ...)
#manual says good acceptance is 0.2-0.4 so if acceptance ratio is higher than 0.4 increase value of parameter
#if acceptance ratio is lower than 0.2 decrease value of parameter

png(paste(output_folder, group_init, ".png", sep=""))
model.N1 <- makeBayouModel(dat ~ lnMass, rjpars = c("theta", "lnMass"),  
                           tree=all_tree.1, dat=dat, pred=pred, SE=0.1, prior=priors, D=DN1[[args[1]]])
dev.off()

mcmc.N1 <- bayou.makeMCMC(all_tree.1, dat, pred=pred, SE=0.1, model=model.N1$model, prior=priors, startpar=model.N1$startpar, 
                          new.dir=output_folder, outname=paste(group_init, iter, sep="_"), plot.freq=NULL)
mcmc.N1$run(runs)

pdf(paste(output_folder, group_init, ".pdf", sep=""))
chain.N1 <- set.burnin(mcmc.N1$load(), 0.5)
save.image(paste(output_folder, group_init, ".RData", sep=""))
#plot(chain.N1, auto.layout=FALSE)
par(mfrow=c(1,1))
plotSimmap.mcmc(chain.N1, burnin = 0.3, pp.cutoff = 0.3)
shiftsumsN1 <- shiftSummaries(chain.N1, mcmc.N1, pp.cutoff=0.3)
save.image(paste(output_folder, group_init, ".RData", sep=""))
plotShiftSummaries(shiftsumsN1, lwd=2, single.plot=TRUE, label.pts=FALSE)
dev.off()
