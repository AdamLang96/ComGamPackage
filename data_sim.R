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
site.b.ptau  <- c(rnorm(n=250, mean = 5, sd=10), rnorm(n=250, mean = 10, sd=10), rnorm(n=250, mean = 25, sd=20), rnorm(n=250, mean = 35, sd=20))
site.b.tau   <- c(rnorm(n=250, mean = 50, sd=20), rnorm(n=250, mean = 100, sd=30), rnorm(n=250, mean = 150, sd=40), rnorm(n=250, mean = 200, sd=50))

sim.dframe.a       <- sim.dframe.b <- sim.dframe[order(sim.dframe$abeta), ]
sim.dframe.b$order <- c(1:1000)
sim.dframe.b$abeta <- sim.dframe.b$abeta + site.b.abeta

sim.dframe.b      <- sim.dframe.b[order(sim.dframe.b$ptau), ]
sim.dframe.b$ptau <- sim.dframe.b$ptau + site.b.ptau

sim.dframe.b     <- sim.dframe.b[order(sim.dframe.b$tau), ]
sim.dframe.b$tau <- sim.dframe.b$tau + site.b.tau

sim.dframe.b       <- sim.dframe.b[order(sim.dframe.b$order),]
sim.dframe.b$order <- NULL

sim.dframe.a$site  <- rep("site-a", 1000)
sim.dframe.a$abeta <- sim.dframe.a$abeta + abeta.add
sim.dframe.a$ptau  <- sim.dframe.a$ptau + ptau.add
sim.dframe.a$tau   <- sim.dframe.a$tau + tau.add

sim.dframe.b$site <- rep("site-b", 1000)
full.sim <- rbind(sim.dframe.a, sim.dframe.b)
full.sim$STUDY <-NULL
full.sim$STUDY <- full.sim$site
full.sim$site  <- NULL

#significant site variance#
check.var <- gam(ptau ~  s(AGE, k=4) + ApoEInd + GENDER + STUDY, data=full.sim, method = "REML")
check.var <- gam(ptau ~  s(AGE, k=4) + ApoEInd + GENDER , data=full.sim, method = "REML")
round(runif(250, min = 0, max = 2000))

####### see how well it improves on standard harmonization
sim.feat <- full.sim[, c("ptau", "tau")]
sim.cov  <- full.sim[, c("AGE", "ApoEInd", "GENDER", "STUDY")]
sim.cov$STUDY <- as.factor(sim.cov$STUDY)
sim.harm <- ComGamHarm(feature.data = sim.feat,
                       covar.data = sim.cov,
                       training.indicies = round(runif(250, min = 0, max = 2000)),
                       smooth.terms = c("AGE"),
                       k.val = c(4),
                       model.diagnostics = TRUE,
                       verbose = TRUE)

feat.harm <- as.data.frame(t(sim.harm$harm.results))
feat.harm <- cbind(feat.harm, sim.cov)


ggplot(full.sim, aes(x=AGE, y=ptau, colour=STUDY)) + geom_point() + geom_smooth(method = "gam", formula = y~s(x, k=4)) + labs(title = "preharm ptau") + ylim(0,100)
ggplot(full.sim, aes(x=AGE, y=tau, colour=STUDY)) + geom_point() + geom_smooth(method = "gam", formula = y~s(x, k=4))  + labs(title = "preharm tau")

ggplot(feat.harm, aes(x=AGE, y=ptau, colour=STUDY)) + geom_point() + geom_smooth(method = "gam", formula = y~s(x, k=4))  + labs(title = "postharm ptau") + ylim(0, 100)
ggplot(feat.harm, aes(x=AGE, y=tau, colour=STUDY)) + geom_point() + geom_smooth(method = "gam", formula = y~s(x, k=4))  + labs(title = "postharm tau")

summary(gam(ptau ~  s(AGE, k=4) + ApoEInd + GENDER + STUDY, data = feat.harm, method = "REML"))
summary(gam(ptau ~  s(AGE, k=4) + ApoEInd + GENDER , data  = feat.harm, method = "REML"))

view.std <- as.data.frame(t(sim.harm$stan.dict$std.data))
view.std <- cbind(view.std, sim.cov)



ggplot(view.std, aes(x=AGE, y=ptau, colour=STUDY)) + geom_point() +geom_smooth(method = "lm")


## manually run nl adj ##
split.std <- split(view.std, view.std$STUDY)
std.adj.df <- data.frame("ptau.site.a" = split.std$`site-a`$ptau,
                         "tau.site.a" = split.std$`site-a`$tau,
                         "ptau.site.b" = split.std$`site-b`$ptau,
                         "tau.site.b" = split.std$`site-b`$tau)

nl.adj.gam <- gam(ptau.site.a ~ s(ptau.site.b, k=4), data = std.adj.df)
summary(nl.adj.gam)

site.a.compare <- predict.gam(nl.adj.gam, type = "response")

nl.adj.df <- data.frame("ptau.site.a" = std.adj.df$ptau.site.a,
                        "ptau.b.mapped" = site.a.compare)
nl.adj.df <- melt(nl.adj.df, measure.vars = c("ptau.site.a", "ptau.b.mapped"))
nl.adj.df <- cbind(nl.adj.df, sim.cov)

ggplot(nl.adj.df, aes(x=AGE, y=value, colour=variable)) + geom_point()
View(t(sim.harm$stan.dict$var.pooled))

std.mn   <-t(sim.harm$stan.dict$stand.mean)[,1]
std.var <- t(sim.harm$stan.dict$var.pooled)[,1]

nl.adj.df$value <- nl.adj.df$value * std.var
nl.adj.df$value <- nl.adj.df$value + std.mn

nl.adj.df <- cbind(nl.adj.df, sim.cov)
ggplot(nl.adj.df, aes(x=AGE, y=value, colour=variable)) + geom_point() + geom_smooth(method = "gam", formula = y~s(x, k=4)) + labs(title = "postharm ptau nl adjusted") + ylim(0, 100)

## looks good ##
summary(gam(data = nl.adj.df, formula = value ~ s(AGE, k=4) + ApoEInd + GENDER + variable, method = "REML"))

