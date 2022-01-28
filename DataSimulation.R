library(simstudy)

set.seed(123) #for reproducibility

#Creating simulated covariates
site <- defData(      varname = "Age",       dist = "normal",      formula = 75  , variance = 25 ,    id = "id")
site <- defData(site, varname = "Sex",       dist = "categorical", formula = ".5 ;.5",                id = "id")
site <- defData(site, varname = "Education", dist = "normal",      formula = 16,   variance = 3,      id = "id")
site <- defData(site, varname = "ApoE",      dist = "categorical", formula = ".5;.5",                 id = "id")
site <- genData(300, site)

#Adding random effects
site <- addCorData(site, idname   = "id", mu  = c(0), sigma   = c(10), rho = .5, corstr = "cs", cnames = c("a0"))
site <- addPeriods(site, nPeriods = 4, idvars = "id", timevarName = "Age")
site$Age <- site$Age + site$period
site1 <- site2 <- site3 <- site

#Create simulated outcome based on covariates

# Site1 
ROI1.sitea <- defDataAdd(varname  = "ROI1", 
                   formula  = "(10000 + a0) + -250*Sex + 20*Education + 350*ApoE  + (-.85)*(Age^2)", 
                   variance = 600000)

ROI2.sitea <- defDataAdd(varname  = "ROI2", 
                         formula  = "(1700 + a0) + -29.45*Sex + 6*Education + 15*ApoE  + (-.095)*(Age^2)", 
                         variance = 18000)

ROI3.sitea <- defDataAdd(varname  = "ROI3", 
                         formula  = "(50000 + a0) + -350.45*Sex + 45*Education + 150*ApoE  + (-2.95)*(Age^2)", 
                         variance = 6000000)

ROI4.sitea <- defDataAdd(varname  = "ROI4", 
                         formula  = "(20000 + a0) + -150.45*Sex + 15*Education + 75*ApoE  + (-1.15)*(Age^2)", 
                         variance = 1000000)

site1  <- addColumns(ROI1.sitea, site1)
site1  <- addColumns(ROI2.sitea, site1)
site1  <- addColumns(ROI3.sitea, site1)
site1  <- addColumns(ROI4.sitea, site1)

# Site2
ROI1.siteb <- defDataAdd(varname  = "ROI1", 
                         formula  = "(10000 + a0) + -250*Sex + 20*Education + 350*ApoE  + (-.85)*(Age^2) + 450", 
                         variance = 750000)

ROI2.siteb <- defDataAdd(varname  = "ROI2", 
                         formula  = "(1700 + a0) + -29.45*Sex + 6*Education + 15*ApoE  + (-.095)*(Age^2) + 150 ", 
                         variance = 12500)

ROI3.siteb <- defDataAdd(varname  = "ROI3", 
                         formula  = "(50000 + a0) + -350.45*Sex + 45*Education + 150*ApoE  + (-2.95)*(Age^2) + 1500", 
                         variance = 5000000)

ROI4.siteb <- defDataAdd(varname  = "ROI4", 
                         formula  = "(20000 + a0) + -150.45*Sex + 15*Education + 75*ApoE  + (-1.15)*(Age^2) + 2500", 
                         variance = 840000)

site2  <- addColumns(ROI1.siteb, site2)
site2  <- addColumns(ROI2.siteb, site2)
site2  <- addColumns(ROI3.siteb, site2)
site2  <- addColumns(ROI4.siteb, site2)


# Site3
ROI1.sitec <- defDataAdd(varname  = "ROI1", 
                         formula  = "(10000 + a0) + -250*Sex + 20*Education + 350*ApoE  + (-.85)*(Age^2) - 675", 
                         variance = 30000)

ROI2.sitec <- defDataAdd(varname  = "ROI2", 
                         formula  = "(1700 + a0) + -29.45*Sex + 6*Education + 15*ApoE  + (-.095)*(Age^2) - 120 ", 
                         variance = 10005)

ROI3.sitec <- defDataAdd(varname  = "ROI3", 
                         formula  = "(50000 + a0) + -350.45*Sex + 45*Education + 150*ApoE  + (-2.95)*(Age^2) - 1750", 
                         variance = 7000000)

ROI4.sitec <- defDataAdd(varname  = "ROI4", 
                         formula  = "(20000 + a0) + -150.45*Sex + 15*Education + 75*ApoE  + (-1.15)*(Age^2) - 1500", 
                         variance = 900000)

site3  <- addColumns(ROI1.sitec, site3)
site3  <- addColumns(ROI2.sitec, site3)
site3  <- addColumns(ROI3.sitec, site3)
site3  <- addColumns(ROI4.sitec, site3)

site2$id <- site2$id + 300
site3$id <- site3$id + 600

site1$STUDY <- "A"
site2$STUDY <- "B"
site3$STUDY <- "C"

full.simulated.data <- rbind(site1, site2, site3)
full.simulated.data <- full.simulated.data[,c("ROI1", "ROI2", "ROI3", "ROI4", "Age", "Sex", "ApoE", "Education", "STUDY", "id")]
full.simulated.data$Sex <- factor(full.simulated.data$Sex)
full.simulated.data$ApoE <- factor(full.simulated.data$ApoE)
full.simulated.data$STUDY <- factor(full.simulated.data$STUDY)
full.simulated.data.cs <- full.simulated.data[!duplicated(full.simulated.data$id),]

#Feature Data
feature.data <- full.simulated.data.cs[,c("ROI1", "ROI2", "ROI3", "ROI4")]
covariate.data <- full.simulated.data.cs[,c("Age","Sex", "ApoE", "Education", "STUDY")]

feature.data.long <- full.simulated.data[,c("ROI1", "ROI2", "ROI3", "ROI4")]
covariate.data.long <- full.simulated.data[,c("Age","Sex", "ApoE", "Education", "STUDY")]


if(FALSE) {

#Plots Preharm
ggplot(full.simulated.data.cs, aes(x=Age, y=ROI1, colour=STUDY)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))
ggplot(full.simulated.data.cs, aes(x=Age, y=ROI2, colour=STUDY)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))
ggplot(full.simulated.data.cs, aes(x=Age, y=ROI3, colour=STUDY)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))
ggplot(full.simulated.data.cs, aes(x=Age, y=ROI4, colour=STUDY)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))



#Harmonize
harmfeats <- ComGamHarm(feature.data = feature.data,
                        covar.data   = covariate.data,
                        eb           = TRUE,
                        parametric   = TRUE,
                        smooth.terms = c("Age"),
                        k.val        = 5)

harmonized.features <- as.data.frame(t(harmfeats$harm.results))

harm.data <- cbind(harmonized.features, covariate.data)
ggplot(harm.data, aes(x=Age, y=ROI4, colour=STUDY)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x), se=FALSE)

}

