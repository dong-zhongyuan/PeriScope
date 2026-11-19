# SuSiE on protein and imaging, coloc per signal pair, MR-signal rule. One result file per pair; finished pairs are skipped.
suppressPackageStartupMessages({library(coloc);library(susieR);library(parallel);library(data.table)})
a <- commandArgs(TRUE); WORK <- a[1]; RES <- a[2]; cores <- as.integer(a[3]); glist <- if(length(a)>3 && nzchar(a[4])) readLines(a[4]) else NULL
plist <- if(length(a)>4 && nzchar(a[5])) readLines(a[5]) else NULL
Lv <- as.integer(Sys.getenv('SUSIE_L','5')); MAXIT <- as.integer(Sys.getenv('SUSIE_MAXIT','200')); PILOT <- Sys.getenv('SUSIE_PILOT_N','0')=='1'; COV <- .95; PUR <- .5
GD <- file.path(WORK,'groups'); dir.create(RES,recursive=TRUE,showWarnings=FALSE)
clean <- function(e) substr(gsub('[\r\n,]+',' ',conditionMessage(e)),1,200)
fit <- function(z,s,b,R,N,k,pos){ ds <- list(beta=setNames(b,k),varbeta=setNames(s^2,k),snp=k,position=pos,type='quant',N=N,sdY=1,LD=R)
  tryCatch(suppressWarnings(suppressMessages(runsusie(ds,repeat_until_convergence=FALSE,maxit=MAXIT,L=Lv,coverage=COV,min_abs_corr=PUR))),error=function(e) e) }
one_pair <- function(pr,gdir,S0,Rf,Gm,inst,gi,prot,sig,cache){
  out <- list(pair_id=pr$pair_id,status=NA_character_,reason=NA_character_,n_snps=NA_integer_,n_protein_cs=NA_integer_,n_imaging_cs=NA_integer_,n_protein_signals_low_overlap=NA_integer_,conv_protein=NA,conv_imaging=NA,lambda_protein=NA_real_,lambda_imaging=NA_real_,min_p_protein=NA_real_,min_p_imaging=NA_real_,n_af_discordant_protein=NA_integer_,n_af_discordant_imaging=NA_integer_,instrument_in_panel=isTRUE(inst$in_panel),instrument_in_snpset=NA,n_protein_cs_qualifying=0L,instrument_in_signal=FALSE,instrument_basis='none',instrument_max_r2_to_cs=NA_real_,H4_susie_mr_signal=NA_real_,H3_susie_mr_signal=NA_real_,H4_susie_mr_signal_p12_1e6=NA_real_,H3_susie_mr_signal_p12_1e6=NA_real_,mr_signal_protein_cs=NA_integer_,mr_signal_imaging_cs=NA_integer_,mr_signal_hit_protein=NA_character_,mr_signal_hit_imaging=NA_character_,H4_susie_any=NA_real_,H3_susie_any=NA_real_,H4_susie_any_p12_1e6=NA_real_,H3_susie_any_p12_1e6=NA_real_,any_protein_cs=NA_integer_,any_imaging_cs=NA_integer_,any_hit_protein=NA_character_,any_hit_imaging=NA_character_,seconds=NA_real_)
  long <- NULL; t0 <- Sys.time()
  fin <- function(st,why=NA_character_){ out$status <<- st; out$reason <<- why; out$seconds <<- round(as.numeric(difftime(Sys.time(),t0,units='secs')),2); list(out=as.data.frame(out,stringsAsFactors=FALSE),long=long) }
  if(pr$prep_status!='ok') return(fin('not_evaluable',pr$prep_reason))
  y <- fread(file.path(gdir,paste0('img_',pr$pair_id,'.csv')),showProgress=FALSE); S <- intersect(S0$k,y$k); isA <- pr$platform=='AGES'
  if(isA){ S <- intersect(S,prot$k) }
  out$n_snps <- length(S); out$min_p_imaging <- if(length(S)) min(y$p[match(S,y$k)],na.rm=TRUE) else NA_real_
  if(length(S)<100) return(fin('not_evaluable',paste0('fewer than 100 panel-matched SNPs (',length(S),')')))
  ii <- match(S,S0$k); y <- y[match(S,y$k),]; pos <- S0$pos[ii]; afp <- S0$af_a1[ii]
  # cached LD for this SNP set
  hit <- 0; for(j in seq_along(cache$e$list)) if(identical(cache$e$list[[j]]$S,S)){hit <- j;break}
  if(hit==0){ R <- Rf[ii,ii,drop=FALSE]; dimnames(R) <- list(S,S); eg <- eigen(R,symmetric=TRUE); R2 <- R; attr(R2,'eigen') <- eg
    cache$e$list <- c(cache$e$list,list(list(S=S,R=R,R2=R2,fit=NULL,lam=NULL))); hit <- length(cache$e$list); if(hit>4){cache$e$list[[1]] <- NULL;hit <- hit-1} }
  ent <- cache$e$list[[hit]]; R <- ent$R; R2 <- ent$R2
  Ny <- if(PILOT) 21000 else as.numeric(pr$n_imaging); Nx <- if(PILOT) 5362 else as.numeric(pr$n_protein)
  out$lambda_imaging <- tryCatch(suppressWarnings(estimate_s_rss(y$z,R2,n=Ny)),error=function(e) NA_real_)
  out$n_af_discordant_imaging <- if(all(is.na(y$eaf_a1))) NA_integer_ else sum(abs(y$eaf_a1-afp)>.2,na.rm=TRUE)
  if(isA){ x <- prot[match(S,prot$k),]; out$min_p_protein <- min(x$p,na.rm=TRUE); out$n_af_discordant_protein <- sum(abs(x$eaf_a1-afp)>.2,na.rm=TRUE)
    if(is.null(ent$fit)){ ent$fit <- fit(x$z,x$s,x$b,R,Nx,S,pos); ent$lam <- tryCatch(suppressWarnings(estimate_s_rss(x$z,R2,n=Nx)),error=function(e) NA_real_); cache$e$list[[hit]] <- ent }
    Sx <- ent$fit; out$lambda_protein <- ent$lam
    if(inherits(Sx,'error')){ m <- clean(Sx); return(fin(if(grepl('converge',m)) 'not_converged' else 'not_evaluable',paste0('protein SuSiE: ',m))) }
    out$conv_protein <- isTRUE(Sx$converged); idx1 <- Sx$sets$cs_index; cs1 <- if(length(idx1)) lapply(Sx$sets$cs,function(v) S[v]) else list()
    out$n_protein_cs <- length(idx1); out$n_protein_signals_low_overlap <- 0L
    if(!length(idx1)) return(fin('no_protein_cs','SuSiE returned no 95% credible set passing purity'))
    bf1 <- Sx$lbf_variable[idx1,,drop=FALSE]; plab <- paste0('L',idx1)
  } else { # UKB: published ST16 signals
    sigs <- split(sig,sig$signal); bf1 <- NULL; cs1 <- list(); plab <- character(0); low <- 0L
    for(sg in sigs){ g <- sg[nzchar(sg$k)&sg$k %in% S,]; mass <- sum(g$cs_pip); if(mass<.5){low <- low+1L;next}
      l <- rep(-50,length(S)); names(l) <- S; pl <- max(sg$lead_pip[1],1e-6); ln <- sg$log10bf[1]*log(10); l[g$k] <- ln+log(pmax(g$cs_pip,1e-6)/pl); bf1 <- rbind(bf1,l); cs1[[length(cs1)+1]] <- sg$k[nzchar(sg$k)]; plab <- c(plab,paste0('S',sg$signal[1])) }
    out$n_protein_signals_low_overlap <- low; out$n_protein_cs <- length(cs1)
    if(!length(cs1)) return(fin('no_protein_cs',if(length(sigs)) 'no ST16 cis signal has >=50% of its credible-set probability among panel-matched window SNPs' else 'no ST16 cis signal'))
  }
  # imaging
  if(!is.finite(out$min_p_imaging) || out$min_p_imaging>=1e-5) return(fin('no_imaging_signal',NA))
  Sy <- fit(y$z,y$s,y$b,R,Ny,S,pos)
  if(inherits(Sy,'error')){ m <- clean(Sy); return(fin(if(grepl('converge',m)) 'not_converged' else 'not_evaluable',paste0('imaging SuSiE: ',m))) }
  out$conv_imaging <- isTRUE(Sy$converged); idx2 <- Sy$sets$cs_index; out$n_imaging_cs <- length(idx2)
  if(!length(idx2)) return(fin('no_imaging_cs','SuSiE returned no 95% credible set passing purity'))
  bf2 <- Sy$lbf_variable[idx2,,drop=FALSE]; rownames(bf1) <- NULL; rownames(bf2) <- NULL
  run <- function(p12){ r <- suppressWarnings(suppressMessages(coloc.bf_bf(bf1,bf2,p1=1e-4,p2=1e-4,p12=p12))); s <- r$summary; if(is.null(s)) NULL else as.data.frame(s) }
  s5 <- run(1e-5); s6 <- run(1e-6); if(is.null(s5)) return(fin('not_evaluable','coloc.bf_bf returned no summary'))
  m6 <- match(paste(s5$idx1,s5$idx2),paste(s6$idx1,s6$idx2))
  # MR-signal rule
  ik <- inst$key; ig <- if(!is.null(gi)) { g <- as.numeric(gi); g[g==255] <- NA; g[is.na(g)] <- mean(g,na.rm=TRUE); g } else NULL
  out$instrument_in_snpset <- ik %in% S; q <- rep(FALSE,length(cs1)); basis <- rep('none',length(cs1)); r2m <- rep(NA_real_,length(cs1))
  for(j in seq_along(cs1)){ cs <- cs1[[j]]; if(ik %in% cs){q[j] <- TRUE;basis[j] <- 'in_credible_set'}
    ix <- match(intersect(cs,S0$k),S0$k); if(length(ix)&&!is.null(ig)&&stats::sd(ig)>0){ rr <- suppressWarnings(cor(ig,Gm[,ix,drop=FALSE]))^2; r2m[j] <- max(rr,na.rm=TRUE); if(!q[j]&&is.finite(r2m[j])&&r2m[j]>=.8){q[j] <- TRUE;basis[j] <- 'r2_ge_0.8_with_credible_set_variant'} } }
  out$n_protein_cs_qualifying <- sum(q); out$instrument_in_signal <- any(q); if(any(is.finite(r2m))) out$instrument_max_r2_to_cs <- max(r2m,na.rm=TRUE); if(any(q)) out$instrument_basis <- if(any(basis=='in_credible_set')) 'in_credible_set' else 'r2_ge_0.8_with_credible_set_variant'
  long <- data.frame(pair_id=pr$pair_id,protein_cs=plab[s5$idx1],imaging_cs=paste0('L',idx2[s5$idx2]),hit_protein=s5$hit1,hit_imaging=s5$hit2,nsnps=s5$nsnps,PP.H0=s5$PP.H0.abf,PP.H1=s5$PP.H1.abf,PP.H2=s5$PP.H2.abf,PP.H3=s5$PP.H3.abf,PP.H4=s5$PP.H4.abf,PP.H3_p12_1e6=s6$PP.H3.abf[m6],PP.H4_p12_1e6=s6$PP.H4.abf[m6],mr_signal=q[s5$idx1],mr_basis=basis[s5$idx1],instrument_r2_to_cs=r2m[s5$idx1],stringsAsFactors=FALSE)
  pick <- function(d){ d <- d[!is.na(d$PP.H4),,drop=FALSE]; if(!nrow(d)) return(NULL); d[which.max(d$PP.H4),,drop=FALSE] }
  a <- pick(long); if(!is.null(a)){ out$H4_susie_any <- a$PP.H4; out$H3_susie_any <- a$PP.H3; out$H4_susie_any_p12_1e6 <- a$PP.H4_p12_1e6; out$H3_susie_any_p12_1e6 <- a$PP.H3_p12_1e6; out$any_protein_cs <- as.integer(sub('^[LS]','',a$protein_cs)); out$any_imaging_cs <- as.integer(sub('^L','',a$imaging_cs)); out$any_hit_protein <- a$hit_protein; out$any_hit_imaging <- a$hit_imaging }
  b <- pick(long[long$mr_signal,,drop=FALSE]); if(!is.null(b)){ out$H4_susie_mr_signal <- b$PP.H4; out$H3_susie_mr_signal <- b$PP.H3; out$H4_susie_mr_signal_p12_1e6 <- b$PP.H4_p12_1e6; out$H3_susie_mr_signal_p12_1e6 <- b$PP.H3_p12_1e6; out$mr_signal_protein_cs <- as.integer(sub('^[LS]','',b$protein_cs)); out$mr_signal_imaging_cs <- as.integer(sub('^L','',b$imaging_cs)); out$mr_signal_hit_protein <- b$hit_protein; out$mr_signal_hit_imaging <- b$hit_imaging }
  fin('ran')
}
do_group <- function(gid){
  gdir <- file.path(GD,gid); rd <- file.path(RES,gid); dir.create(rd,showWarnings=FALSE,recursive=TRUE)
  P <- as.data.frame(fread(file.path(gdir,'pairs.csv'),showProgress=FALSE,colClasses=list(character='pair_id'))); if(!is.null(plist)) P <- P[P$pair_id %in% plist,,drop=FALSE]
  todo <- P[!file.exists(file.path(rd,paste0(P$pair_id,'.csv'))),,drop=FALSE]; if(!nrow(todo)) return(0L)
  S0 <- as.data.frame(fread(file.path(gdir,'snps.csv'),showProgress=FALSE)); m <- nrow(S0); inst <- jsonlite_read(file.path(gdir,'instrument.json'))
  Gm <- matrix(as.integer(readBin(file.path(gdir,'G.bin'),'raw',n=m*503)),nrow=503,ncol=m); Gm[Gm==255L] <- NA; storage.mode(Gm) <- 'double'; mu <- colMeans(Gm,na.rm=TRUE); w <- which(is.na(Gm),arr.ind=TRUE); if(nrow(w)) Gm[w] <- mu[w[,2]]
  Rf <- suppressWarnings(cor(Gm)); Rf[!is.finite(Rf)] <- 0; diag(Rf) <- 1
  gi <- if(file.exists(file.path(gdir,'inst.bin'))) as.integer(readBin(file.path(gdir,'inst.bin'),'raw',n=503)) else NULL
  plat <- todo$platform[1]; prot <- if(plat=='AGES') as.data.frame(fread(file.path(gdir,'protein.csv'),showProgress=FALSE)) else NULL; sig <- if(plat=='UKB') as.data.frame(fread(file.path(gdir,'signals.csv'),showProgress=FALSE)) else NULL
  cache <- new.env(); cache$e <- new.env(); cache$e$list <- list()
  for(i in seq_len(nrow(todo))){ pr <- as.list(todo[i,]); r <- tryCatch(one_pair(pr,gdir,S0,Rf,Gm,inst,gi,prot,sig,cache),error=function(e) list(out=data.frame(pair_id=pr$pair_id,status='not_evaluable',reason=paste0('worker error: ',clean(e)),stringsAsFactors=FALSE),long=NULL))
    base <- file.path(rd,pr$pair_id); if(!is.null(r$long)&&nrow(r$long)){ tl <- paste0(base,'.sig.tmp',Sys.getpid()); fwrite(r$long,tl,showProgress=FALSE); file.rename(tl,paste0(base,'.sig.csv')) }
    tf <- paste0(base,'.tmp',Sys.getpid()); fwrite(r$out,tf,showProgress=FALSE); file.rename(tf,paste0(base,'.csv')) }
  nrow(todo)
}
jsonlite_read <- function(f){ tx <- paste(readLines(f,warn=FALSE),collapse=''); gs <- function(k){ m <- regmatches(tx,regexec(paste0('"',k,'": *"([^"]*)"'),tx))[[1]]; if(length(m)==2) m[2] else NA_character_ }; list(key=gs('key'),in_panel=grepl('"in_panel": *true',tx)) }
gs <- if(is.null(glist)) list.dirs(GD,full.names=FALSE,recursive=FALSE) else glist; gs <- gs[file.exists(file.path(GD,gs,'DONE'))]
set.seed(7); gs <- sample(gs); cat('groups',length(gs),'L',Lv,'maxit',MAXIT,'pilotN',PILOT,'\n')
res <- mclapply(gs,function(g) tryCatch(do_group(g),error=function(e){ writeLines(clean(e),file.path(RES,paste0(g,'.error'))); -1L }),mc.cores=cores,mc.preschedule=FALSE)
cat('DONE groups',length(gs),'errors',sum(vapply(res,function(x) !is.numeric(x)||x<0,TRUE)),'\n')
