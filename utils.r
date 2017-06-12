library(pdist)
# source: http://r-sig-geo.2731867.n2.nabble.com/alpha-hull-ahull-to-polygon-shapefile-td7342734.html

ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){ 
  require(alphahull) 
  require(maptools) 
  if (class(x) != "ahull"){ 
    stop("x needs to be an ahull class object") 
  } 
  # Extract the edges from the ahull object as a dataframe 
  xdf <- as.data.frame(x$arcs) 
  # Remove all cases where the coordinates are all the same       
  xdf <- subset(xdf,xdf$r > 0) 
  res <- NULL 
  if (nrow(xdf) > 0){ 
    # Convert each arc to a line segment 
    linesj <- list() 
    prevx<-NULL 
    prevy<-NULL 
    j<-1 
    for(i in 1:nrow(xdf)){ 
      rowi <- xdf[i,] 
      v <- c(rowi$v.x, rowi$v.y) 
      theta <- rowi$theta 
      r <- rowi$r 
      cc <- c(rowi$c1, rowi$c2) 
      # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment. 
      ipoints <- 2 + round(increment * (rowi$theta / 2),0) 
      # Calculate coordinates from arc() description for ipoints along the arc. 
      angles <- anglesArc(v, theta) 
      seqang <- seq(angles[1], angles[2], length = ipoints) 
      x <- round(cc[1] + r * cos(seqang),rnd) 
      y <- round(cc[2] + r * sin(seqang),rnd) 
      # Check for line segments that should be joined up and combine their coordinates 
      if (is.null(prevx)){ 
        prevx<-x 
        prevy<-y 
      } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){ 
        if (i == nrow(xdf)){ 
          #We have got to the end of the dataset 
          prevx<-append(prevx,x[2:ipoints]) 
          prevy<-append(prevy,y[2:ipoints]) 
          prevx[length(prevx)]<-prevx[1] 
          prevy[length(prevy)]<-prevy[1] 
          coordsj<-cbind(prevx,prevy) 
          colnames(coordsj)<-NULL 
          # Build as Line and then Lines class 
          linej <- Line(coordsj) 
          linesj[[j]] <- Lines(linej, ID = as.character(j)) 
        } else { 
          prevx<-append(prevx,x[2:ipoints]) 
          prevy<-append(prevy,y[2:ipoints]) 
        } 
      } else { 
        # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset. 
        prevx[length(prevx)]<-prevx[1] 
        prevy[length(prevy)]<-prevy[1] 
        coordsj<-cbind(prevx,prevy) 
        colnames(coordsj)<-NULL 
        # Build as Line and then Lines class 
        linej <- Line(coordsj) 
        linesj[[j]] <- Lines(linej, ID = as.character(j)) 
        j<-j+1 
        prevx<-NULL 
        prevy<-NULL 
      } 
    } 
    # Promote to SpatialLines 
    lspl <- SpatialLines(linesj) 
    # Convert lines to polygons 
    # Pull out Lines slot and check which lines have start and end points that are the same 
    lns <- slot(lspl, "lines") 
    polys <- sapply(lns, function(x) { 
      crds <- slot(slot(x, "Lines")[[1]], "coords") 
      identical(crds[1, ], crds[nrow(crds), ]) 
    }) 
    # Select those that do and convert to SpatialPolygons 
    polyssl <- lspl[polys] 
    list_of_Lines <- slot(polyssl, "lines") 
    sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string) 
    # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame 
    hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID")) 
    areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area")) 
    df <- data.frame(hid,areas) 
    names(df) <- c("HID","Area") 
    rownames(df) <- df$HID 
    res <- SpatialPolygonsDataFrame(sppolys, data=df)
    print(res@data$Area)
    res <- res[which(res@data$Area > 0),] 
  }   
  return(res) 
} 

# source: https://stat.ethz.ch/pipermail/r-sig-geo/2009-May/005781.html

owin2Polygons <- function(x, id="1") {
  stopifnot(is.owin(x))
  x <- as.polygonal(x)
  closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
  pieces <- lapply(x$bdry,
                   function(p) {
                     Polygon(coords=closering(cbind(p$x,p$y)),
                             hole=is.hole.xypolygon(p))  })
  z <- Polygons(pieces, id)
  return(z)
}

tess2SP <- function(x) {
  stopifnot(is.tess(x))
  y <- tiles(x)
  nam <- names(y)
  z <- list()
  for(i in seq(y))
    z[[i]] <- owin2Polygons(y[[i]], nam[i])
  return(SpatialPolygons(z))
}

owin2SP <- function(x) {
  stopifnot(is.owin(x))
  y <- owin2Polygons(x)
  z <- SpatialPolygons(list(y))
  return(z)
}

# this is for Maryland
lat_long_to_xy = function(lat,long) {
  library(rgdal)
  library(sp)
  data = data.frame(long=long, lat=lat)
  coordinates(data) <- ~ long+lat
  proj4string(data) <- CRS("+init=epsg:4326")
  xy = data.frame(spTransform(data, CRS("+init=epsg:2804")))
  setnames(xy,c("x","y"))
  return(xy[,c("x","y")])
}

max_v = function(x,x0) vapply(x,function(a) max(a,x0), x0)

epanech = function(x, lengthscale) {
  return(max_v((0.75 * (1 - (x/3.1415926/lengthscale)^2)/3.1415926)/lengthscale,0.0));
}

quartic = function(x,lengthscale) {
  return(max_v((0.9375 * (1 - (x/lengthscale)^2)^2)/lengthscale,0.0))
}
epanech2 = function(x, lengthscale) {
  return(max_v((0.75 * (1 - (x/lengthscale)^2))/lengthscale,0.0));
}

daniell = function(x, lengthscale) {
  return(max_v(0.2387324*(1-abs(x)/3.1415926/lengthscale)/lengthscale,0.0));
}
gaussian2 = function(x, lengthscale) {
  return(exp(-(x^2/lengthscale^2)/2)/sqrt(pi*2)/lengthscale)
  
}
gaussian = function(x, lengthscale) {
  return(0.3989423*exp(-.5 * x^2/lengthscale^2)/lengthscale);
}

intensity.estimate.separable = function(xyt, lengthscale.s, lengthscale.t) {
  return(apply(xyt, 1, function(x) sum(epanech2(pdist(x[1:2],xyt[,1:2])@dist, lengthscale.s)) *
                 sum(epanech2(pdist(matrix(x[3]),matrix(xyt[,3]))@dist, lengthscale.t))/nrow(xyt) ))
}

intensity.estimate = function(xyt, lengthscale.s, lengthscale.t, kernel=epanech2) {
  return(apply(xyt, 1, function(x) sum(kernel(pdist(x[1:2],xyt[,1:2])@dist, lengthscale.s) *
                                         (kernel(pdist(matrix(x[3]),matrix(xyt[,3]))@dist, lengthscale.t)))))
}

