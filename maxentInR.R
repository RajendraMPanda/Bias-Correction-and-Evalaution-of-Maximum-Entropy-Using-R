library(sp)
library(raster)
library(dismo)
library(rgdal)
library(rJava)
library(ENMeval)
setwd(system.file(package="dismo"))

#import list of environmental variables (.img files)
files <- list.files(paste(system.file(package="dismo"), "/ASCIFiles", sep=""), pattern='asc', full.names=TRUE)
pred <- stack(files)
projection(pred) <- CRS('+proj=longlat')

#import presence only data
pres <- read.csv(paste(system.file(package="dismo"), "/Species.csv", sep=""))
coordinates(pres) <- ~LON+LAT
projection(pres) <- CRS('+proj=longlat')

#fix random background points
set.seed(0)
bgpt<-randomPoints(pred, 10000)
nrpt <- nrow(bgpt)

#partition background points into test and train sets
spt<- sample(nrpt, 0.25 * nrpt)
back_train <- bgpt[-spt, ]
back_test <- bgpt[spt, ]

#k-fold partitioning of presence-only dataset
fold <- kfold(pres, k=4) #4-fold partition
pres_test <- pres[fold == 1, ] #25% held_up data for testing
pres_train <- pres[fold != 1, ] #75% data for training 

#spatial sorting bias (Hijmans, 2012)
sbpt <- ssb(pres_test, back_test, pres_train)
sbpt[,1] / sbpt[,2]

#Removing SSB (sample sorting bias)
i <- pwdSample(pres_test, back_test, pres_train, n=1, tr=0.33)
pres_test_pwd <- pres_test[!is.na(i[,1]), ]
back_test_pwd <- back_test[na.omit(as.vector(i)), ]
sb2 <- ssb(pres_test_pwd, back_test_pwd, pres_train)
sb2[1]/ sb2[2]

# improves memory
options(java.parameters = "-mex2g") 
jarfile <- paste(system.file(package="dismo"), "/java/maxent", sep="") #import javafile for maxent run

#maxent run and tune parameters
mxt <- maxent(x=pred, p=pres_train, remove.duplicates=T, outputfiletype = asc, appendtoresultsfile = TRUE, bias=ct_bias_asc2.asc, writeplotdata = TRUE, m = 5000, pictures = TRUE, 
              replicates = 15, replicatetype = Subsample, writeplotdata = TRUE, args=c('outputformat=Logistic','jackknife=true', "responsecurves=true" ),
             path = paste(system.file(package="dismo"), '/Output_LC', sep=''))

# predict to entire dataset
t2<-predict(mxt,pred)
plot(t2)

#write raster
writeRaster(t2, filename = "t2.grd")
MyRaster <- raster("t2.grd")

#write raster in .asc format
t2<-writeRaster(MyRaster, file="t2.asc", format="ascii")

#project raster
projection(t2<- CRS('+proj=longlat'))

#evaluate maxent prediction
e1 <- evaluate(mxt, p=pres_test, a=back_test, x=pred)
e1

#bias corrected evaluation
e11 <- evaluate(mxt, p=pres_test_pwd, a=back_test_pwd, x=pred)
e11
tr<-threshold(e1)
tr
tr1<-threshold(e11)
tr1

# null geographic model to get cAUC for maxent
distm <- geoDist(pres_train, lonlat=TRUE)
colnames(back_test) <- c("LON", "LAT")
e<-evaluate(p=pres_test, a=back_test, model=mxt, x=pred)
e_pwd<-evaluate(p=pres_test_pwd, a=back_test_pwd, model=mxt, x=pred)
e_null <- evaluate(distm, p=pres_test, a=back_test)
colnames(back_test_pwd) <- c("LON", "LAT")
e_null_pwd <- evaluate(distm, p=pres_test_pwd, a=back_test_pwd)

e <- list(e)
e_pwd <- list(e_pwd)
e_null <- list(e_null)
e_null_pwd <- list(e_null_pwd)
#auc value calculation
(auc <- sapply(e, function(x) {slot(x, 'auc' ) } ) )
mean(auc)
auc_adj <- sapply( e_pwd, function(x) {slot(x, 'auc' ) } ) 
mean(auc_adj)
(null_auc <- sapply( e_null, function(x) {slot(x, 'auc' ) } ) )
mean(null_auc)
(null_auc_adj <- sapply( e_null_pwd, function(x) {slot(x, 'auc' ) } ) )
mean(null_auc_adj)

#import shape file
myshp<-shapefile(file.choose("/india_wgs"),stringsAsFactors = FALSE,verbose=FALSE)
#projection
projection(myshp) <- CRS('+proj=longlat')
#add boundary
plot(myshp, add=TRUE, border="black")
#add points
points(pres, col='blue', pch=20, cex=0.25)
# presence points and define colour
points(pres, pch=20, cex=0.5, col="darkgreen")
# add pseudo-absence points and define colour
points(pres, cex=0.5, col="darkorange3")
#calculate AICc for Maxent (Warren and Seifert, 2011)
get.params(mxt)
calc.aicc(127, pres,t2)


