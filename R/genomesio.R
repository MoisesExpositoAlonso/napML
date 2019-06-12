####************************************************************************####
#### SIMULATE ####
#' Simulate a big genome matrix with filebacked system
#'
#' @param backingpath
#' @param filename
#' @param n
#' @param m
#' @param force
#'
#' @return
#' @export
#'
#' @examples
XBMsimulate<-function(backingpath="../databig",filename="example",n,m,force=F){
  # Create a filebacked big matrix in R, and fill it
  # with 2s (homozygous for minor) at a MAF simulated from runif 0-0.5
  # To read matrix, attach.big.matrix("descriptionfile.desc")
  if(force){
    system(paste("rm",paste0(backingpath,"/",filename,".bk")))
  }
  x <- filebacked.big.matrix(n, m, type='double', init= 0,
                             backingpath = backingpath,
                             backingfile=paste0(filename,".bk"),
                             descriptorfile=paste0(filename,".desc")
                            )
  x
  BMsimulate(x@address)
  return(x)
}

## Need to also produce a Map and BIM files
MAPwrite<-function(x,path="../databig/example"){
  #### Write bim map
  dmap<-data.frame(CHR=1,
                   SNP=paste(sep="_",1,1:ncol(x)),
                   X=0,
                   POS=1:ncol(x),
                   R="C",A="A")
  write.table(row.names = F,col.names = F, quote = F,
              file=paste0(path, ".map"),
              dmap[,1:4] # the last columns are for a .bim file, but gets reconstructed based on map and ped using plink later
              )
  return(TRUE)
}

#### write 012 PED
RAW012write<-function(x,path="../databig/example"){
  BMwrite012(x@address,paste0(path,".012")) # only needed for my own purposes
  return(TRUE)
}

PEDwrite<-function(x,path="../databig/example"){
  tmp1<-tempfile()
  tmp2<-tempfile()
  BMwritePED(x@address,tmp1)
  write.table(row.names = F,col.names = F, quote = F,
              file=tmp2,
              data.frame(FID=1:nrow(x),IID=1:nrow(x),PAT=0,MAT=0,SEX=0,PHENO=-9)
              )
  system(paste("paste ", tmp2, tmp1," > ", paste0(path,".ped")))
  return(TRUE)
}
#
# PEDwrite<-function(x,path="../databig/example"){
#   BMwritePED(x@address,paste0(path,".ped.tmp")) # to convert to bed
#   dum<-data.frame(FID=1:nrow(x),IID=1:nrow(x),PAT=0,MAT=0,SEX=0,PHENO=-9)
#   write.table(row.names = F,col.names = F, quote = F,
#               file="tmp.fam",
#               dum
#               )
#   system(paste("cut -d ' ' -f-6 tmp.fam > ", paste0(path,".ped.tmp2")) )
#   system("paste ", paste0(path,".ped.tmp2"), paste0(path,".ped.tmp")," > ", paste0(path,".ped"))
#   return(TRUE)
# }

BEDmake<-function(path="../databig/example"){
  system(paste("plink --noweb --file ",path," --make-bed --out ",path) )
}

####************************************************************************####



readplink<-function(bfile = "databig/toy",force=T){
    require(bigmemory)
    require(data.table)

    famfile<-paste0(bfile,".fam")
    bedfile<-paste0(bfile,".bed")
    bimfile<-paste0(bfile,".bim")
    mapfile<-paste0(bfile,".map")
    rdafile <- paste0(bfile, ".rda")
    bkfile<-paste0(bfile,".bk")
    descfile<-paste0(bfile,".desc")

    if(file.exists(rdafile) & !force){
      message("Genome already read, loading from RDA file")
      thedata<-readRDS(rdafile)
      thedata$G<-attach.big.matrix(descfile)
    }else{
    # if (!all(c(basename(g012file), basename(mapfile), basename(famfile)) %in%
    #     list.files(dirname(g012file)))) {
    #     stop("Not all bim bed and fam files were found!")
    # }
    NAMES.BIM <- c("chromosome", "marker.ID", "genetic.dist",
        "physical.pos", "allele1", "allele2")
    NAMES.FAM <- c("family.ID", "sample.ID", "paternal.ID", "maternal.ID",
        "sex", "affection")
    message("Reading the FAM file...")
    fam <- data.table::fread(famfile)
    N<-nrow(fam)
    colnames(fam) <- NAMES.FAM
    message("Reading the BIM/MAP file...")
    if(file.exists(bimfile)){
      bim <- data.table::fread(bimfile)
    }else{
      bim <- data.table::fread(mapfile)
    }
    p<-nrow(bim)
    colnames(bim) <- NAMES.BIM
    message("Reading the BED file...")
    if(force & any(file.exists(c(bkfile,descfile))) ){
      system(paste("rm",bkfile))
      system(paste("rm",descfile))
    }
    x <- filebacked.big.matrix(N, p, type='double', init= 0,
                             backingpath = dirname(bfile),
                             backingfile=basename(bkfile),
                             descriptorfile=basename(descfile)
                            )
    BMreadbed(bedfile,x@address,N,p)
    # save structure
    thedata <- structure(list(genotypes = describe(x),
                              desc=descfile,
                               fam = fam,
                               map = bim,
                               savedIn = rdafile))
    saveRDS(object = thedata, file = rdafile)
    message("The genome file is stored in ", rdafile)
    thedata$G<-attach.big.matrix(descfile)
    }

    return(thedata)
}


readplink012<-function (g012file = "data-raw/515g.012", backingfile = "genomes",
    backingpath = "databig/")
{
    require(bigmemory)
    require(data.table)
    # bimfile <- sub("\\.012$", ".bim", g012file)
    famfile <- sub("\\.012$", ".fam", g012file)
    # mapfile <- sub("\\.012$", ".map", g012file)
    mapfile <- sub("\\.012$", ".bim", g012file)

    if (!all(c(basename(g012file), basename(mapfile), basename(famfile)) %in%
        list.files(dirname(g012file)))) {
        stop("Not all bim bed and fam files were found!")
    }
    NAMES.MAP <- c("chromosome", "marker.ID", "genetic.dist",
        "physical.pos", "allele1", "allele2")
    NAMES.FAM <- c("family.ID", "sample.ID", "paternal.ID", "maternal.ID",
        "sex", "affection")
    message("Reading the FAM file...")
    fam <- data.table::fread(famfile)
    colnames(fam) <- NAMES.FAM
    message("Reading the BIM/MAP file...")
    bim <- data.table::fread(mapfile)
    colnames(bim) <- NAMES.MAP
    message("Reading the 012 file...")

    geno <- bigmemory::read.big.matrix(filename = g012file, sep = " ",
        # type = "integer", // will dramatically affect implementation
        type = "double",
        backingfile = paste0(backingfile, ".bk"),
        backingpath = backingpath,
        descriptorfile = paste0(backingfile,".desc"),
        skip = 1)

    rda <- paste0(file.path(backingpath, backingfile), ".rda")
    snp_list <- structure(list(genotypes = describe(geno), fam = fam,
        map = bim, savedIn = rda))

    saveRDS(object = snp_list, file = rda)
    message("The genome file is stored in ", rda)
    return(snp_list)
}


read_n_subset_genome<-function(genomefile='databig/genome.rda',selectedsnps)
{
    genomes<-readRDS(genomefile)
    map<-genomes$map
    map$SNP<-paste0(map$chromosome, "_",map$physical.pos)

    Go<-bigmemory::attach.big.matrix(
      paste0(tools::file_path_sans_ext(genomefile),'.desc'))

    sg<- data.frame(Go[,which(map$SNP %in% selectedsnps)])
    sg<-sg / 2
    head(sg)
    colnames(sg) <- selectedsnps

    return(sg)
}

####************************************************************************####

get_genomematrix<-function (size = 1000, positions = NULL, type = c("random", "window",
    "top"), start = 1)
{
    genomes <- readRDS("databig/genomes.rda")
    G <- attachgenomes(genomes)
    Map <- genomes$map
    if (!is.null(positions)) {
        message("subsetting the provided positions")
        stopifnot(is.numeric(positions))
        X = inputena.mat(G[, positions])
    }
    else if (type == "random") {
        message("subsetting ", size, " random positions ")
        X = inputena.mat(G[, sample(1:ncol(G), size = size)])
    }
    else if (type == "window") {
        message("subsetting a window of ", size, " base pairs starting at ",
            start)
        X = inputena.mat(G[, start:(start + size)])
    }
    else if (type == "top") {
        X <- get_topSNPs()
    }
    else {
        stop("None of the possibilities were given")
    }
    return(X)
}

readattachG<-function(genomesfile="databig/genome.desc", doattach=F)
{
  genomes<-readRDS(genomesfile)
  if(doattach==T){
    genomes<-attachgenomes(genomes)
  }
  return(genomes)
}

attachgenomes<-function (genomes)
{
    genomes$g <- bigmemory::attach.big.matrix(genomes$genotypes)
    return(genomes)
}

# read_plink<-function (root, snps = NULL, impute = c("none", "avg", "random"),
#     verbose = FALSE)
# {
#     proot <- path.expand(root)
#     impute <- match.arg(impute)
#     impute_int <- switch(impute, none = 0L, avg = 1L, random = 2L)
#     bedfile <- paste(proot, ".bed", sep = "")
#     famfile <- paste(proot, ".fam", sep = "")
#     bimfile <- paste(proot, ".bim", sep = "")
#     bim <- read.table(bimfile, header = FALSE, sep = "", stringsAsFactors = FALSE)
#     fam <- read.table(famfile, header = FALSE, sep = "", stringsAsFactors = FALSE)
#     geno <- .Call("read_plink", PACKAGE = "gws", bedfile, famfile,
#         impute_int, verbose)
#     colnames(geno) <- bim[, 2]
#     rownames(geno) <- paste(fam[, 1], fam[, 2], sep = ":")
#     list(bed = geno, fam = fam, bim = bim)
# }

