

mort.data <- function(nev="baseline",y=2,ext=".dat") {
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  r <- NULL;
  nevek <- c("surv","pred","starv","y.surv","y.pred","y.starv",
             "ye.surv","ye.pred","ye.starv");
  if(file.exists(f.nev)) {
    command <- paste(.path.package("Rmoult"),"/scripts/mort.awk -v y=",y," ",
                     f.nev,sep="");
    zz1 <- pipe(command);
    r <- scan(zz1,quiet=TRUE);
    close(zz1);
    names(r) <- nevek;
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
  }
  return(r);
}

mort.ser <- function(y=2,m.path=".",patt="\\.dat$",f.list=NULL) {
  if(is.null(f.list)) {
    f.list <- list.files(path=m.path,pattern=patt,full.names=TRUE);
    f.list <- gsub("\\.dat$","",f.list);
  }
  r <- NULL;
  for(i in f.list){
    r <- rbind(r,mort.data(nev=i,y=y,ext=NULL));
  }
  rownames(r) <- f.list;
  r <- data.frame(r);
  r;
}

breed.eff.data <- function(nev="baseline",y=2,ext=".dat") {
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  r <- NULL;
  nevek <- c("n","breed.eff","breed.att","prod");
  if(file.exists(f.nev)) {
    command <- paste(.path.package("Rmoult"),"/scripts/breed_eff.awk -v y=",
                     y," ", f.nev,sep="");
    zz1 <- pipe(command);
    r <- scan(zz1,quiet=TRUE);
    close(zz1);
    names(r) <- nevek;
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
  }
  return(r);
}

breed.eff.ser <- function(y=2,m.path=".",patt="\\.dat$",f.list=NULL) {
  if(is.null(f.list)) {
    f.list <- list.files(path=m.path,pattern=patt,full.names=TRUE);
    f.list <- gsub("\\.dat$","",f.list);
  }
  r <- NULL;
  for(i in f.list){
    r <- rbind(r,breed.eff.data(nev=i,y=y));
  }
  rownames(r) <- f.list;
  r <- data.frame(r);
  r;
}
