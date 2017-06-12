library(data.table)
library(alphahull)
library(spatstat)
library(rgeos)
#library(splancs)

source("utils.r")

data = fread("WashingtonDC_Incidents_2006-2013_Raw_Data_2.csv")
data = data[,3:6,with=FALSE]
setnames(data,c("datetime","Type","long","lat"))
data = data[lat > 37]
head(data)

data$date = as.Date(data$datetime, format="%m/%d/%y")
data$timestamp = as.POSIXct(data$datetime, format="%m/%d/%y %H:%M")
data$hour = format(data$timestamp,"%H") 
data$dow = format(data$timestamp,"%u") 

data = data[complete.cases(data$timestamp),]
n.orig = nrow(data)
data = unique(data)
n.duplicates = n.orig - nrow(data)

Time = as.matrix(as.numeric(data$timestamp))
Time = (Time - min(Time)) / 60 / 60  # hours
X = as.matrix(lat_long_to_xy(data$lat, data$long))

X[,1] =(X[,1] - min(X[,1])) / 1000  # kilometers
X[,2] =(X[,2] - min(X[,2])) / 1000 # kilometers
plot(X)
data = data.frame(X=X[,1], Y=X[,2], T=Time, date=data$date)
data = unique(data)

xy = data[,1:2]
xy = unique(xy)
dd = apply(as.matrix(dist(xy)) + diag(nrow(xy)),1,min)

plot(xy, pch=20, col=rgb(0,0,0,.4))
points(xy[dd>=.3,],col="red", pch=20)
xy = xy[dd<.3,]

x.ah <- ahull(xy[,1], xy[,2], alpha=.5) 
plot(xy)
plot(x.ah,col='purple',add=T,asp=1,cex=0, pch=4)

my.win = ah2sp(x.ah)
plot(my.win, col="lightgrey", pbg="white") 
my.win = gBuffer(my.win,width=.2)
plot(my.win)
points(xy)

plot(my.win)
points(data[,1:2], pch=20, col=rgb(0,0,0,.4))
my.owin = as.owin(my.win)

inside.window <- inside.owin(x=data[,1],y=data[,2],w=my.owin)
xyt = data[inside.window,1:3]
data = data[inside.window,]
xyt[,1] = xyt[,1] + rnorm(nrow(xyt)) * 1e-6 # add jitter
xy = ppp(xyt[,1],xyt[,2], window=my.owin)
plot(xy,pch=20,col=rgb(0,0,0,.4))

doy = format(data$date,"%m-%d")
# code holidays as: 
# July 2, 3, 4, 5, 6, December 29, 30, 31, January 1, 2

holidays = c("07-01","07-02", "07-03", "07-04", "07-05", "07-06", "12-29", "12-30", "12-31", "01-01", "01-02")
data$holiday = doy %in% holidays

# make sure it's sorted, since the Hawkes model code assumes it is
ii = order(data$T)
data = data[ii,]
xyt = xyt[ii,]
xy = xy[ii,]

save(data,xyt,xy,my.win,my.owin,file="WashingtonDC-agls.rdata")
