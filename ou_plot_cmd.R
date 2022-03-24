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
runs <- 10000000

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

dat <- all_data.1$ln_tarsus
names(dat) <- all_data.1$species

pred <- data.frame(lnMass = all_data.1$ln_mass)
rownames(pred) <- all_data.1$species

priors <- make.prior(all_tree.1, plot.prior = T,
                     dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm",
                                dsb="dsb", dk="cdpois", dtheta="dnorm"),
                     param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                dbeta_lnMass=list(mean=mean(all_data.1$ln_mass), sd=sd(all_data.1$ln_mass)),
                                dk=list(lambda=nrow(all_data.1)*0.02, kmax=nrow(all_data.1)*0.05),
                                dsb=list(bmax=1,prob=1), dtheta=list(mean=0.33, sd=0.5))
) #slope is theta and intercept is beta
DN1 <- list(alpha=2, sig2=2, beta_lnMass=0.1, k=1, theta=0.1, slide=1)

png(paste(output_folder, group_init, ".png", sep=""))
model.N1 <- makeBayouModel(dat ~ lnMass, rjpars = c("theta"),  
                           tree=all_tree.1, dat=dat, pred=pred, SE=0.1, prior=priors, D=DN1)
dev.off()

mcmc.N1 <- bayou.makeMCMC(all_tree.1, dat, pred=pred, SE=0.1, model=model.N1$model, prior=priors, startpar=model.N1$startpar, 
                          new.dir=output_folder, outname=paste(group_init, iter, sep="_"), plot.freq=NULL)
mcmc.N1$run(runs)

#chain.N1 <- set.burnin(mcmc.N1$load(), 0.3)
#plot(chain.N1, auto.layout=FALSE)
#par(mfrow=c(1,1))
#plotSimmap.mcmc(chain.N1, burnin = 0.3, pp.cutoff = 0.1)
#shiftsumsN1 <- shiftSummaries(chain.N1, mcmc.N1, pp.cutoff=0.3)
#plotShiftSummaries(shiftsumsN1, lwd=2, single.plot=TRUE, label.pts=FALSE)

