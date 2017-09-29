O.moult <- function(nev="baseline",ss=list(a=0:7,e=0:2,f2=-6:10,f1=-6:10,r=0:12,
                                  t=0:51)){
  f.nev <- paste(nev,".bin",sep="");
  zz <- file(f.nev,"rb");
  n <- readBin(zz,what="int",n=1);
  seek(zz,where=2*8*n,origin="current");
  n <- readBin(zz,what="int",n=1);
  seek(zz,where=8*n+4*n,origin="current");
  m1 <- readBin(zz,what="int",n=n);
  m2 <- readBin(zz,what="int",n=n);
  close(zz);
  states <- gen.states(ss);
  states[["m1"]] <- m1;
  states[["m2"]] <- m2;
  invisible(states);
}


O.pol <- function(nev="baseline",ss=list(a=0:7,e=0:2,f2=-6:10,f1=-6:10,r=0:12,
                                  t=0:51)){
  f.nev <- paste(nev,".bin",sep="");
  zz <- file(f.nev,"rb");
  n <- readBin(zz,what="int",n=1);
  seek(zz,where=2*8*n,origin="current");
  n <- readBin(zz,what="int",n=1);
  u <- readBin(zz,what="double",n=n);
  act <- readBin(zz,what="int",n=n);
  m1 <- readBin(zz,what="int",n=n);
  m2 <- readBin(zz,what="int",n=n);
  close(zz);
  states <- gen.states(ss);
  states[["u"]] <- u;
  states[["act"]] <- act;
  states[["m1"]] <- m1;
  states[["m2"]] <- m2;
  invisible(states);
}


O.sim <- function(y=NULL,nev="moult") {
                                        #quick read in dat file
  f.nev <- paste(nev,".dat",sep="");
  gzipped <- FALSE;
  if (!file.exists(f.nev)) {
    f.nev <- paste(nev,".dat.gz",sep="");
    if(file.exists(f.nev)) {
      system(paste("gunzip ",f.nev));
      gzipped <- TRUE;
    }
  }
  ds <- NA;
  if(file.exists(f.nev)) {
    if(is.null(y)) {
      zz2 <- pipe(paste("awk '{if ($1!=\"year\") print $0}' ",f.nev, sep=""));
    } else {
      zz2 <- pipe(paste("awk '{if ($1==",y,") print $0}' ",f.nev, sep=""));
    }
    ds <- scan(zz2,what=list(year=integer(0),week=integer(0),bird=integer(0),
                     event=character(0),moult1=character(0),moult2=character(0),
                     action=character(0),res=integer(0),fea1=integer(0),
                     ml1=integer(0),fea2=integer(0),ml2=integer(0),
                     brood=integer(0),exper=integer(0),uval=double(0)));
    close(zz2);
  }
  invisible(ds);
}

O.yrv <- function(nev="baseline"){
  f.nev <- paste(nev,".yrv",sep="");
  r <- read.table(f.nev);
  names(r) <- c("rf","FOOD","n.fora","G0","theta0","N.for0","G1","theta1",
                "N.for1","G2","theta2","N.for2"); 
  invisible(r);
}


O.pop <- function(nev="moult",ee=3,ll=1) {
  f.nev <- paste(gsub("\\.C$","",nev),".pop",sep="");
  if(file.exists(f.nev)) {
    pop <- read.table(f.nev,header=FALSE);
    v.nevek <- c("alive","breed","moultP","moultS");
    szamok <- paste(rep(1:ee,ll),rep(1:ll,rep(ee,ll)),sep=".");
    nevek <- paste(v.nevek,rep(szamok,rep(length(v.nevek),length(szamok))),
                   sep="");
    colnames(pop) <- c("week","total",nevek);
    invisible(pop);
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
  }
}

O.rv <- function(nev="baseline",ss=list(a=0:7,e=0:2,f2=-6:10,f1=-6:10,r=0:12,
                                  t=0:52)){
  f.nev <- paste(nev,".bin",sep="");
  zz <- file(f.nev,"rb");
  n <- readBin(zz,what="int",n=1);
  rv <- readBin(zz,what="double",n=n);
  close(zz);
  states <- gen.states(ss);
  states[["rv"]] <- rv;
  invisible(states);
}

gen.states <- function (...) 
{
    nargs <- length(args <- list(...))
    if (nargs == 1 && is.list(a1 <- args[[1]])) 
        nargs <- length(args <- a1)
    cargs <- args
    nmc <- paste("Var", 1:nargs, sep = "")
    nm <- names(args)
    if (is.null(nm)) 
        nm <- nmc
    if (any(ng0 <- nchar(nm) > 0)) 
        nmc[ng0] <- nm[ng0]
    names(cargs) <- nmc
    rep.fac <- 1
    orep <- final.len <- prod(sapply(args, length))
    for (i in 1:nargs) {
        x <- args[[i]]
        nx <- length(x)
        orep <- orep/nx
        x <- rep(rep(x, rep(rep.fac, nx)), orep)
        if (!is.factor(x) && is.character(x)) 
            x <- factor(x, levels = unique(x))
        cargs[[i]] <- x
        rep.fac <- rep.fac * nx
    }
    invisible(cargs);
}
