### generate simulated data to test new ComBat approach ###

# run harmonization to calculate site variance parameters as reference #
testdata  <- read.csv("/Users/adamgabriellang/Desktop/CSFBiomarkerProject/OutputFiles/tidybiomarkerdataPREHARM.csv")
testdata  <- subset(testdata, baseline == TRUE)
csf.train <- FindIndices(testdata, "DIAGNOSIS", 1)
csf.feat  <- testdata[ ,c("ABETA", "PTAU", "TAU")]
csf.cov   <- testdata[ ,c("AGE", "ApoEInd", "GENDER", "STUDY")]
csf.cov$ApoEInd <-as.factor(csf.cov$ApoEInd)
csf.cov$GENDER  <- as.factor(csf.cov$GENDER)
csf.harm  <- ComGamHarm(feature.data = csf.feat,
                        covar.data   = csf.cov,
                        training.indicies = csf.train,
                        smooth.terms      = c("AGE"),
                        k.val   = c(4),
                        verbose = TRUE,
                        model.diagnostics = TRUE)

rownames(csf.cov) <- 1:nrow(csf.cov)
rownames(csf.feat) <- 1:nrow(csf.feat)
adni1rows <- which(csf.cov$STUDY=="ADNI1" , arr.ind = TRUE)


data.std     <- as.data.frame(t(csf.harm$stan.dict$std.data))[adni1rows,]
std.mean     <- as.data.frame(t(csf.harm$stan.dict$stand.mean))[adni1rows,]
var.pooled   <- as.data.frame(t(csf.harm$stan.dict$var.pooled))[adni1rows,]
gamma.effect <- csf.harm$shift.scale.params$gamma.hat
delta.effect <- csf.harm$shift.scale.params$delta.hat


#simulate covariates based on distribution from data set
length(which(csf.cov$GENDER==1)) / 2033

sim.age  <- round(rnorm(n=1000, mean = 72, sd=9))
sim.apoe <- sample(x = c(0,1,2),
                   size = 1000,
                   replace = TRUE,
                   prob = c(0.5435317,  0.3566158, 0.09985243))
sim.apoe <- as.factor(sim.apoe)
sim.gender <- sample(x=c(1,2),
                     size = 1000,
                     replace = TRUE,
                     prob = c(0.53, 0.47))

sim.gender <-as.factor(sim.gender)

newapo<- rep(0, 1000)
sim.study <- rep("ADNI1", 1000)

sim.dframe <- data.frame("AGE" = sim.age,
                         "ApoEInd" = sim.apoe,
                         "GENDER" = sim.gender,
                         "STUDY" = sim.study)
sim.abeta.gam <- csf.harm$models.list$ABETA
sim.ptau.gam  <- csf.harm$models.list$PTAU
sim.tau.gam  <- csf.harm$models.list$TAU

sim.abeta <- predict.gam(sim.abeta.gam, newdata = sim.dframe, type = "response")
sim.ptau <- predict.gam(sim.ptau.gam, newdata = sim.dframe, type = "response")
sim.tau <- predict.gam(sim.tau.gam, newdata = sim.dframe, type = "response")

abeta.add <- rnorm(n=1000, mean = 0, sd=60)
ptau.add <- rnorm(n=1000, mean = 0, sd=5)
tau.add <- rnorm(n=1000, mean = 0, sd=38)

sim.dframe$abeta <- sim.abeta
sim.dframe$ptau  <- sim.ptau
sim.dframe$tau   <- sim.tau 

site.b.abeta <- c(rnorm(n=250, mean = 100, sd=100), rnorm(n=250, mean = 250, sd=100), rnorm(n=250, mean = 350, sd=100), rnorm(n=250, mean = 450, sd=100))
site.b.ptau <- c(rnorm(n=250, mean = 5, sd=10), rnorm(n=250, mean = 10, sd=10), rnorm(n=250, mean = 25, sd=20), rnorm(n=250, mean = 35, sd=20))
site.b.tau <- c(rnorm(n=250, mean = 50, sd=20), rnorm(n=250, mean = 100, sd=30), rnorm(n=250, mean = 150, sd=40), rnorm(n=250, mean = 200, sd=50))

sim.dframe.a <- sim.dframe.b <- sim.dframe[order(sim.dframe$abeta), ]
sim.dframe.b$order <- c(1:1000)
sim.dframe.b$abeta <- sim.dframe.b$abeta + site.b.abeta

sim.dframe.b <- sim.dframe.b[order(sim.dframe.b$ptau), ]
sim.dframe.b$ptau <- sim.dframe.b$ptau + site.b.ptau

sim.dframe.b <- sim.dframe.b[order(sim.dframe.b$tau), ]
sim.dframe.b$tau <- sim.dframe.b$tau + site.b.tau

sim.dframe.b <- sim.dframe.b[order(sim.dframe.b$order),]
sim.dframe.b$order <- NULL

sim.dframe.a$site <- rep("site-a", 1000)
sim.dframe.a$abeta <- sim.dframe.a$abeta + abeta.add
sim.dframe.a$ptau <- sim.dframe.a$ptau + ptau.add
sim.dframe.a$tau <- sim.dframe.a$tau + tau.add

sim.dframe.b$site <- rep("site-b", 1000)
testframe <- rbind(sim.dframe.a, sim.dframe.b)


