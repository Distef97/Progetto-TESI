library(MTS)
library(fable)
library(tsibble)
library(dplyr)
library(tidyr)
#################
#generazione matrici
#################
genera_4_matrici_3x3_rango_ridotto <- function(rango = 2) {
  stopifnot(rango %in% 1:2)
  
  # Base ortonormale V del sottospazio colonna comune
  V <- qr.Q(qr(matrix(rnorm(3 * rango), nrow = 3)))
  
  # Base completa per generare matrici stabili
  Q <- qr.Q(qr(matrix(rnorm(3 * 3), nrow = 3)))
  
  generate_stable_matrix <- function() {
    eigvals <- c(runif(rango, min = -0.9, max = 0.9), rep(0, 3 - rango))
    A_tmp <- Q %*% diag(eigvals) %*% t(Q)       # stabile
    A_proj <- V %*% t(V) %*% A_tmp              # proiettata nel sottospazio
    return(A_proj)
  }
  
  A1 <- generate_stable_matrix()
  A2 <- generate_stable_matrix()
  A3 <- generate_stable_matrix()
  A4 <- generate_stable_matrix()
  
  C <- cbind(A1, A2, A3, A4)
  
  return(list(A1 = A1, A2 = A2, A3 = A3, A4 = A4, concatenata = C))
}

genera_2matrici <- function(grado = 3) {
  
  stopifnot(grado >= 2)
  
  
  
  # La radice che sarà ripetuta almeno due volte
  
  radice_dup <- runif(1, -0.9, 0.9)
  
  
  
  # Due radici uguali (le restanti random)
  
  if (grado > 2) {
    
    autovalori_AR <- c(radice_dup, radice_dup, runif(grado - 2, -0.9, 0.9))
    
  } else {
    
    autovalori_AR <- c(radice_dup, radice_dup)
    
  }
  
  
  
  # Almeno una radice di B uguale a quella duplicata di A
  
  if (grado > 1) {
    
    autovalori_B <- c(radice_dup, runif(grado - 1, -0.9, 0.9))
    
  } else {
    
    autovalori_B <- radice_dup
    
  }
  
  
  
  # Matrici ortogonali casuali per diagonalizzazione
  
  Q1 <- qr.Q(qr(matrix(rnorm(grado^2), nrow = grado)))
  
  Q2 <- qr.Q(qr(matrix(rnorm(grado^2), nrow = grado)))
  
  
  
  # Costruzione delle matrici
  
  AR <- Q1 %*% diag(autovalori_AR) %*% t(Q1)
  
  MA  <- Q2 %*% diag(autovalori_B) %*% t(Q2)
  
  
  
  return(list(A1 = AR, A2 = MA, autovalori_AR = autovalori_AR, autovalori_MA = autovalori_B))
  
}
###########
#PRima funz###
#############
gen_varma_noid<- function(seed=123, maxl=7, part="varma2"){
  set.seed(seed)
  
  if(part=="varma2"){res<-genera_2matrici()
  ar<-as.vector(res$A1)
  ma<-as.vector(res$A2)}
  
  if(part=="varma3"){
    res<-genera_4_matrici_3x3_rango_ridotto()
    
    ar<-as.vector(cbind(res$A1,res$A2))
    ma<-as.vector(cbind(res$A3,res$A4))}
  if(part=="var"){
    phi1 <- matrix(runif(9), 3, 3)
    phi1[upper.tri(phi1, TRUE)] <- 0
    res<-c()
    res$A1<-phi1
    phi2<-phi1%*%phi1
    ar<-as.vector(phi1)
    ma<-as.vector(cbind(phi1,phi2))
  }
  
  w<-nrow(res$A1)
  s<-length(ar)/(w*w)
  l<-length(ma)/(w*w)
  sig<-matrix(0.2,w,w)
  diag(sig)<-1
  mod1<- structure(
    list(
      ar  = array(ar,c(w,w,s)),
      ma  = array(-ma, c(w,w, l)),
      cov = matrix(sig, w,w)
    ),
    class = "varma"
  )
  
  
  #controllo stazionarietà
  roots <- inv_roots(mod1, plot = F)
  stationary <- all(Mod(roots$ar) < 1)
  invertible <- all(Mod(roots$ma) < 1)
  #c(stationary&invertible)# T,T
  #controllo non identificatbilità
  #!is_identified(mod1)#T
  
  while(!(stationary&invertible&!is_identified(mod1))){
    res<-genera_4_matrici_3x3_rango_ridotto()
    
    ar<-as.vector(cbind(res$A1,res$A2))
    ma<-as.vector(cbind(res$A3,res$A4))
    w<-nrow(res$A1)
    s<-length(ar)/(w*w)
    sig<-matrix(0.2,w,w)
    diag(sig)<-1
    mod1<- structure(
      list(
        ar  = array(ar,c(w,w,s)),
        ma  = array(-ma, c(w,w, s)),
        cov = matrix(sig, w,w)
      ),
      class = "varma"
    )
    
    
    #controllo stazionarietà
    roots <- inv_roots(mod1, plot = F)
    stationary <- all(Mod(roots$ar) < 1)
    invertible <- all(Mod(roots$ma) < 1)
    c(stationary&invertible)# T,T
    #controllo non identificatbilità
    !is_identified(mod1)#T
  }
  
  roots <- inv_roots(mod1, plot = T)
  
  irf1 <- irf(
    mod1,
    maxlag = maxl,
    orth = "none",
    plot = FALSE
  )
  results<-list(
    mod1=mod1,
    irf1=irf1,
    roots=roots,
    identified=is_identified(mod1),
    w=w,
    s=s,
    l=l,
    maxl=maxl
  )
  cat("To continue, save: mod1, maxl, w, s,l, irf1")
  return(results)
}


############
#SImulazione IH
##############
simulvarmanoidih<-function(NV=1000, simulazioni=5){
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  risultati_varma <- vector("list", simulazioni) # Pre-allocazione per efficienza
  for (o in 1:simulazioni) {
    
    
    Y <- sim_varma(mod1, n = NV)
    #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
    
    #Stima del VARMA(secondo hannan o altri processi)
    ih_fit <- fit_varma_ihr(Y, s, l, TRUE)
    #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
    #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
    # simte con: VARMACpp, ss_varma
    roots <- inv_roots(ih_fit, plot = F)
    while (!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))) {
      Y <- sim_varma(mod1, n = NV)
      ih_fit <- fit_varma_ihr(Y, s, s, TRUE)
      roots <- inv_roots(ih_fit, plot = F)}
    
    risultati_varma[[o]]<-ih_fit
    
    irf2 <- irf(
      ih_fit,
      maxlag = maxl,
      orth = "none",
      plot = FALSE
    )
    
    
    
    psi_err<-c()
    for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
    
    diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
    
    irf_err[o,,1]<-as.numeric(psi_err)
    irf2<-as.vector(irf2)
    irf2_save<-cbind(irf2_save,irf2)
  }
  
  #nrow(irf_err[,,1]) #nrow-> sono ogni simulazione
  means1<-colMeans(irf_err[,,1])#errori della stima delle
  RMSE1<-sqrt(colMeans((irf_err[,,1])^2))
  #Irf con stima del modello di Hannan
  #length(means1)
  valori<-c()
  for (i in 1:(maxl+1)) {
    valori<-cbind(valori, mean(means1[indices[i,]]))
  }
  #valori
  RMS<-c()
  for (i in 1:(maxl+1)) {
    RMS<-cbind(RMS, mean(RMSE1[indices[i,]]))
  }
  #varianze errori
  var1<-apply(irf_err[,,1], 2, var)
  
  values<-c()
  for (i in 1:(maxl+1)) {
    values<-cbind(values, mean(var1[indices[i,]]))
  }
  #values
  #coverage
  #da rifareeee
  irf_min<-c()
  irf_max<-c()
  long<-length(irf2_save[,1])
  for (i in 1:long) {
    irf_min<-cbind(irf_min,min(irf2_save[i,]))
    irf_max<-cbind(irf_max,max(irf2_save[i,]))
  }
  irf1_val<-as.vector(irf1)
  irf_min<-as.vector(irf_min)
  irf_max<-as.vector(irf_max)
  coverage<-mean(irf1_val<irf_max&irf1_val>irf_min)
  
  results<-list(
    NV=NV,
    simulazioni=simulazioni,
    valori=valori,
    RMSE=RMS,
    values=values,
    coverage=coverage,
    differenze_irf=diff_irf,
    diff_irf=mean(diff_irf),
    irf2=irf2_save,
    varianza=mean(apply(irf2_save, 1,var)),
    risultati_varma=risultati_varma
    
    
    
    
    
  )
  return(results)
}
#################
#Simulazione .net
#################

simulvarmanoidnet<-function(NV=1000, simulazioni=5){
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  risultati_varma <- vector("list", simulazioni) # Pre-allocazione per efficienza
  for (o in 1:simulazioni) {
    
    
    Y <- sim_varma(mod1, n = NV)
    #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
    
    #Stima del VARMA(secondo hannan o altri processi)
    ih_fit <- fit_varma_net(Y, s, l, TRUE)
    #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
    #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
    # simte con: VARMACpp, ss_varma
    roots <- inv_roots(ih_fit, plot = F)
    while (!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))) {
      Y <- sim_varma(mod1, n = NV)
      ih_fit <- fit_varma_net(Y, s, s, TRUE)
      roots <- inv_roots(ih_fit, plot = F)}
    
    risultati_varma[[o]]<-ih_fit
    
    irf2 <- irf(
      ih_fit,
      maxlag = maxl,
      orth = "none",
      plot = FALSE
    )
    
    
    
    psi_err<-c()
    for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
    
    diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
    
    irf_err[o,,1]<-as.numeric(psi_err)
    irf2<-as.vector(irf2)
    irf2_save<-cbind(irf2_save,irf2)
  }
  
  #nrow(irf_err[,,1]) #nrow-> sono ogni simulazione
  means1<-colMeans(irf_err[,,1])#errori della stima delle
  RMSE1<-sqrt(colMeans((irf_err[,,1])^2))
  #Irf con stima del modello di Hannan
  #length(means1)
  valori<-c()
  for (i in 1:(maxl+1)) {
    valori<-cbind(valori, mean(means1[indices[i,]]))
  }
  #valori
  RMS<-c()
  for (i in 1:(maxl+1)) {
    RMS<-cbind(RMS, mean(RMSE1[indices[i,]]))
  }
  #varianze errori
  var1<-apply(irf_err[,,1], 2, var)
  
  values<-c()
  for (i in 1:(maxl+1)) {
    values<-cbind(values, mean(var1[indices[i,]]))
  }
  #values
  #coverage
  #da rifareeee
  irf_min<-c()
  irf_max<-c()
  long<-length(irf2_save[,1])
  for (i in 1:long) {
    irf_min<-cbind(irf_min,min(irf2_save[i,]))
    irf_max<-cbind(irf_max,max(irf2_save[i,]))
  }
  irf1_val<-as.vector(irf1)
  irf_min<-as.vector(irf_min)
  irf_max<-as.vector(irf_max)
  coverage<-mean(irf1_val<irf_max&irf1_val>irf_min)
  
  results<-list(
    NV=NV,
    simulazioni=simulazioni,
    valori=valori,
    RMSE=RMS,
    values=values,
    coverage=coverage,
    differenze_irf=diff_irf,
    diff_irf=mean(diff_irf),
    irf2=irf2_save,
    varianza=mean(apply(irf2_save, 1,var)),
    risultati_varma=risultati_varma
    
    
    
    
    
  )
  return(results)
}
#################
#Simulazione DMA
#################

simulvarmanoiddma<-function(NV=1000, simulazioni=5){
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  risultati_varma <- vector("list", simulazioni) # Pre-allocazione per efficienza
  for (o in 1:simulazioni) {
    
    
    Y <- sim_varma(mod1, n = NV)
    #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
    
    #Stima del VARMA(secondo hannan o altri processi)
    ih_fit <- fit_varma_dma(Y, s, l, TRUE)
    #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
    #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
    # simte con: VARMACpp, ss_varma
    roots <- inv_roots(ih_fit, plot = F)
    while (!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))) {
      Y <- sim_varma(mod1, n = NV)
      ih_fit <- fit_varma_dma(Y, s, s, TRUE)
      roots <- inv_roots(ih_fit, plot = F)}
    
    risultati_varma[[o]]<-ih_fit
    
    irf2 <- irf(
      ih_fit,
      maxlag = maxl,
      orth = "none",
      plot = FALSE
    )
    
    
    
    psi_err<-c()
    for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
    
    diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
    
    irf_err[o,,1]<-as.numeric(psi_err)
    irf2<-as.vector(irf2)
    irf2_save<-cbind(irf2_save,irf2)
  }
  
  #nrow(irf_err[,,1]) #nrow-> sono ogni simulazione
  means1<-colMeans(irf_err[,,1])#errori della stima delle
  RMSE1<-sqrt(colMeans((irf_err[,,1])^2))
  #Irf con stima del modello di Hannan
  #length(means1)
  valori<-c()
  for (i in 1:(maxl+1)) {
    valori<-cbind(valori, mean(means1[indices[i,]]))
  }
  #valori
  RMS<-c()
  for (i in 1:(maxl+1)) {
    RMS<-cbind(RMS, mean(RMSE1[indices[i,]]))
  }
  #varianze errori
  var1<-apply(irf_err[,,1], 2, var)
  
  values<-c()
  for (i in 1:(maxl+1)) {
    values<-cbind(values, mean(var1[indices[i,]]))
  }
  #values
  #coverage
  #da rifareeee
  irf_min<-c()
  irf_max<-c()
  long<-length(irf2_save[,1])
  for (i in 1:long) {
    irf_min<-cbind(irf_min,min(irf2_save[i,]))
    irf_max<-cbind(irf_max,max(irf2_save[i,]))
  }
  irf1_val<-as.vector(irf1)
  irf_min<-as.vector(irf_min)
  irf_max<-as.vector(irf_max)
  coverage<-mean(irf1_val<irf_max&irf1_val>irf_min)
  
  results<-list(
    NV=NV,
    simulazioni=simulazioni,
    valori=valori,
    RMSE=RMS,
    values=values,
    coverage=coverage,
    differenze_irf=diff_irf,
    diff_irf=mean(diff_irf),
    irf2=irf2_save,
    varianza=mean(apply(irf2_save, 1,var)),
    risultati_varma=risultati_varma
    
    
    
    
    
  )
  return(results)
}

#################
#Simulazione FKF
#################

simulvarmanoidfkf<-function(NV=1000, simulazioni=5){
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  
  for (o in 1:simulazioni) {
    tryCatch({
      print(o)
      
      Y <- sim_varma(mod1, n = NV)
      #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
      
      #Stima del VARMA(secondo hannan o altri processi)
      ih_fit <- tryCatch({fit_varma_fkf(Y, s, l, TRUE)
      }, error = function(e) {
        message(paste("Errore nella simulazione", o, "- Tentativo :", e$message))
        return(NULL) # Ritorna NULL se fallisce
      })
      #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
      #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
      # simte con: VARMACpp, ss_varma
      roots <- inv_roots(ih_fit, plot = F)
      while (!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))) {
        Y <- sim_varma(mod1, n = NV)
        ih_fit <- tryCatch({fit_varma_fkf(Y, s, l, TRUE)
        }, error = function(e) {
          message(paste("Errore nella simulazione", o, "- Tentativo :", e$message))
          return(NULL) # Ritorna NULL se fallisce
        })
        while (is.null(ih_fit)) {
          Y <- sim_varma(mod1, n = NV)
          ih_fit <- tryCatch({fit_varma_fkf(Y, s, l, TRUE)
          }, error = function(e) {
            message(paste("Errore nella simulazione", o, "- Tentativo :", e$message))
            return(NULL) # Ritorna NULL se fallisce
          })
        }
        roots <- inv_roots(ih_fit, plot = F)}
      
      risultati_varma[[o]]<-ih_fit
      
      irf2 <- irf(
        ih_fit,
        maxlag = maxl,
        orth = "none",
        plot = FALSE
      )
      
      
      
      psi_err<-c()
      for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
      
      diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
      
      irf_err[o,,1]<-as.numeric(psi_err)
      irf2<-as.vector(irf2)
      irf2_save<-cbind(irf2_save,irf2)
    }, error = function(e) {
      
      # --- BLOCCO "ERROR": COSA FARE SE SI VERIFICA UN ERRORE ---
      
      # 1. (Opzionale) Stampa un messaggio per sapere cosa è successo
      cat("AVVISO: Errore in un'iterazione della simulazione. Messaggio:", e$message, "\n")
      
      # 2. Restituisci NULL per segnalare il fallimento
      
      
    })}
  
  
  #nrow(irf_err[,,1]) #nrow-> sono ogni simulazione
  means1<-colMeans(irf_err[,,1])#errori della stima delle
  RMSE1<-sqrt(colMeans((irf_err[,,1])^2))
  #Irf con stima del modello di Hannan
  #length(means1)
  valori<-c()
  for (i in 1:(maxl+1)) {
    valori<-cbind(valori, mean(means1[indices[i,]]))
  }
  #valori
  RMS<-c()
  for (i in 1:(maxl+1)) {
    RMS<-cbind(RMS, mean(RMSE1[indices[i,]]))
  }
  #varianze errori
  var1<-apply(irf_err[,,1], 2, var)
  
  values<-c()
  for (i in 1:(maxl+1)) {
    values<-cbind(values, mean(var1[indices[i,]]))
  }
  #values
  #coverage
  #da rifareeee
  irf_min<-c()
  irf_max<-c()
  long<-length(irf2_save[,1])
  for (i in 1:long) {
    irf_min<-cbind(irf_min,min(irf2_save[i,]))
    irf_max<-cbind(irf_max,max(irf2_save[i,]))
  }
  irf1_val<-as.vector(irf1)
  irf_min<-as.vector(irf_min)
  irf_max<-as.vector(irf_max)
  coverage<-mean(irf1_val<irf_max&irf1_val>irf_min)
  
  results<-list(
    NV=NV,
    simulazioni=simulazioni,
    valori=valori,
    RMSE=RMS,
    values=values,
    coverage=coverage,
    differenze_irf=diff_irf,
    diff_irf=mean(diff_irf),
    irf2=irf2_save,
    varianza=mean(apply(irf2_save, 1,var)),
    risultati_varma=risultati_varma
    
    
    
    
    
  )
  return(results)
}

#################
#Simulazione KFAS
#################


simulvarmanoidkfas<-function(NV=1000, simulazioni=5){
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  for (o in 1:simulazioni) {
    print(o)
    
    
    
    Y <- sim_varma(mod1, n = NV)
    #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
    
    #Stima del VARMA(secondo hannan o altri processi)
    ih_fit <- tryCatch({fit_varma_kfas(Y, s, l, TRUE)
    }, error = function(e) {
      message(paste("Errore nella simulazione", o, "- Tentativo :", e$message))
      return(NULL) # Ritorna NULL se fallisce
    })
    #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
    #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
    # simte con: VARMACpp, ss_varma
    roots <- inv_roots(ih_fit, plot = F)
    while (!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))) {
      Y <- sim_varma(mod1, n = NV)
      ih_fit <- tryCatch({fit_varma_kfas(Y, s, l, TRUE)
      }, error = function(e) {
        message(paste("Errore nella simulazione", o, "- Tentativo :", e$message))
        return(NULL) # Ritorna NULL se fallisce
      })
      while (is.null(ih_fit)) {
        Y <- sim_varma(mod1, n = NV)
        ih_fit <- tryCatch({fit_varma_kfas(Y, s, l, TRUE)
        }, error = function(e) {
          message(paste("Errore nella simulazione", o, "- Tentativo :", e$message))
          return(NULL) # Ritorna NULL se fallisce
        })
      }
      roots <- inv_roots(ih_fit, plot = F)}
    
    risultati_varma[[o]]<-ih_fit
    
    irf2 <- irf(
      ih_fit,
      maxlag = maxl,
      orth = "none",
      plot = FALSE
    )
    
    
    
    psi_err<-c()
    for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
    
    diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
    
    irf_err[o,,1]<-as.numeric(psi_err)
    irf2<-as.vector(irf2)
    irf2_save<-cbind(irf2_save,irf2)
  }
  
  #nrow(irf_err[,,1]) #nrow-> sono ogni simulazione
  means1<-colMeans(irf_err[,,1])#errori della stima delle
  RMSE1<-sqrt(colMeans((irf_err[,,1])^2))
  #Irf con stima del modello di Hannan
  #length(means1)
  valori<-c()
  for (i in 1:(maxl+1)) {
    valori<-cbind(valori, mean(means1[indices[i,]]))
  }
  #valori
  RMS<-c()
  for (i in 1:(maxl+1)) {
    RMS<-cbind(RMS, mean(RMSE1[indices[i,]]))
  }
  #varianze errori
  var1<-apply(irf_err[,,1], 2, var)
  
  values<-c()
  for (i in 1:(maxl+1)) {
    values<-cbind(values, mean(var1[indices[i,]]))
  }
  #values
  #coverage
  #da rifareeee
  irf_min<-c()
  irf_max<-c()
  long<-length(irf2_save[,1])
  for (i in 1:long) {
    irf_min<-cbind(irf_min,min(irf2_save[i,]))
    irf_max<-cbind(irf_max,max(irf2_save[i,]))
  }
  irf1_val<-as.vector(irf1)
  irf_min<-as.vector(irf_min)
  irf_max<-as.vector(irf_max)
  coverage<-mean(irf1_val<irf_max&irf1_val>irf_min)
  
  results<-list(
    NV=NV,
    simulazioni=simulazioni,
    valori=valori,
    RMSE=RMS,
    values=values,
    coverage=coverage,
    differenze_irf=diff_irf,
    diff_irf=mean(diff_irf),
    irf2=irf2_save,
    varianza=mean(apply(irf2_save, 1,var)),
    risultati_varma=risultati_varma
    
    
    
    
    
  )
  return(results)
}

simulvarmanoidcpp<-function(NV=1000, simulazioni=5){
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  for (o in 1:simulazioni) {
    print(o)
    
    Y <- sim_varma(mod1, n = NV)
    #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
    
    #Stima del VARMA(secondo hannan o altri processi)
    ih_fit <- fit_varma_cpp(Y, s, l, TRUE)
    #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
    #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
    # simte con: VARMACpp, ss_varma
    roots <- inv_roots(ih_fit, plot = F)
    while (!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))) {
      Y <- sim_varma(mod1, n = NV)
      ih_fit <- fit_varma_cpp(Y, s, s, TRUE)
      roots <- inv_roots(ih_fit, plot = F)}
    
    risultati_varma[[o]]<-ih_fit
    
    irf2 <- irf(
      ih_fit,
      maxlag = maxl,
      orth = "none",
      plot = FALSE
    )
    
    
    
    psi_err<-c()
    for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
    
    diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
    
    irf_err[o,,1]<-as.numeric(psi_err)
    irf2<-as.vector(irf2)
    irf2_save<-cbind(irf2_save,irf2)
  }
  
  #nrow(irf_err[,,1]) #nrow-> sono ogni simulazione
  means1<-colMeans(irf_err[,,1])#errori della stima delle
  RMSE1<-sqrt(colMeans((irf_err[,,1])^2))
  #Irf con stima del modello di Hannan
  #length(means1)
  valori<-c()
  for (i in 1:(maxl+1)) {
    valori<-cbind(valori, mean(means1[indices[i,]]))
  }
  #valori
  RMS<-c()
  for (i in 1:(maxl+1)) {
    RMS<-cbind(RMS, mean(RMSE1[indices[i,]]))
  }
  #varianze errori
  var1<-apply(irf_err[,,1], 2, var)
  
  values<-c()
  for (i in 1:(maxl+1)) {
    values<-cbind(values, mean(var1[indices[i,]]))
  }
  #values
  #coverage
  #da rifareeee
  irf_min<-c()
  irf_max<-c()
  long<-length(irf2_save[,1])
  for (i in 1:long) {
    irf_min<-cbind(irf_min,min(irf2_save[i,]))
    irf_max<-cbind(irf_max,max(irf2_save[i,]))
  }
  irf1_val<-as.vector(irf1)
  irf_min<-as.vector(irf_min)
  irf_max<-as.vector(irf_max)
  coverage<-mean(irf1_val<irf_max&irf1_val>irf_min)
  
  results<-list(
    NV=NV,
    simulazioni=simulazioni,
    valori=valori,
    RMSE=RMS,
    values=values,
    coverage=coverage,
    differenze_irf=diff_irf,
    diff_irf=mean(diff_irf),
    irf2=irf2_save,
    varianza=mean(apply(irf2_save, 1,var)),
    risultati_varma=risultati_varma
    
    
    
    
    
  )
  return(results)
}

simulvarmanoidMTS<-function(NV=1000, simulazioni=5){
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  for (o in 1:simulazioni) {
    print(o)
    
    Y <- sim_varma(mod1, n = NV)
    #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
    
    #Stima del VARMA(secondo hannan o altri processi)
    ih_fit <- convert_mts_varma(MTS::VARMA(Y, s, l, TRUE))
    #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
    #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
    # simte con: VARMACpp, ss_varma
    roots <- inv_roots(ih_fit, plot = F)
    while (!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))) {
      Y <- sim_varma(mod1, n = NV)
      ih_fit <- convert_mts_varma(MTS::VARMA(Y, s, s, TRUE))
      roots <- inv_roots(ih_fit, plot = F)}
    
    risultati_varma[[o]]<-ih_fit
    
    irf2 <- irf(
      ih_fit,
      maxlag = maxl,
      orth = "none",
      plot = FALSE
    )
    
    
    
    psi_err<-c()
    for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
    
    diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
    
    irf_err[o,,1]<-as.numeric(psi_err)
    irf2<-as.vector(irf2)
    irf2_save<-cbind(irf2_save,irf2)
  }
  
  #nrow(irf_err[,,1]) #nrow-> sono ogni simulazione
  means1<-colMeans(irf_err[,,1])#errori della stima delle
  RMSE1<-sqrt(colMeans((irf_err[,,1])^2))
  #Irf con stima del modello di Hannan
  #length(means1)
  valori<-c()
  for (i in 1:(maxl+1)) {
    valori<-cbind(valori, mean(means1[indices[i,]]))
  }
  #valori
  RMS<-c()
  for (i in 1:(maxl+1)) {
    RMS<-cbind(RMS, mean(RMSE1[indices[i,]]))
  }
  #varianze errori
  var1<-apply(irf_err[,,1], 2, var)
  
  values<-c()
  for (i in 1:(maxl+1)) {
    values<-cbind(values, mean(var1[indices[i,]]))
  }
  #values
  #coverage
  #da rifareeee
  irf_min<-c()
  irf_max<-c()
  long<-length(irf2_save[,1])
  for (i in 1:long) {
    irf_min<-cbind(irf_min,min(irf2_save[i,]))
    irf_max<-cbind(irf_max,max(irf2_save[i,]))
  }
  irf1_val<-as.vector(irf1)
  irf_min<-as.vector(irf_min)
  irf_max<-as.vector(irf_max)
  coverage<-mean(irf1_val<irf_max&irf1_val>irf_min)
  
  results<-list(
    NV=NV,
    simulazioni=simulazioni,
    valori=valori,
    RMSE=RMS,
    values=values,
    coverage=coverage,
    differenze_irf=diff_irf,
    diff_irf=mean(diff_irf),
    irf2=irf2_save,
    varianza=mean(apply(irf2_save, 1,var)),
    risultati_varma=risultati_varma
    
    
    
    
    
  )
  return(results)
}

simulvarmanoidvar<-function(NV=1000, simulazioni=5){
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  for (o in 1:simulazioni) {
    print(o)
    
    Y <- sim_varma(mod1, n = NV)
    #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
    
    #Stima del VARMA(secondo hannan o altri processi)
    ih_fit <- varma_from_var(Y, s, l, TRUE)
    #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
    #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
    # simte con: VARMACpp, ss_varma
    roots <- inv_roots(ih_fit, plot = F)
    while (!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))) {
      Y <- sim_varma(mod1, n = NV)
      ih_fit <- varma_from_var(Y, s, s, TRUE)
      roots <- inv_roots(ih_fit, plot = F)}
    
    risultati_varma[[o]]<-ih_fit
    
    irf2 <- irf(
      ih_fit,
      maxlag = maxl,
      orth = "none",
      plot = FALSE
    )
    
    
    
    psi_err<-c()
    for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
    
    diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
    
    irf_err[o,,1]<-as.numeric(psi_err)
    irf2<-as.vector(irf2)
    irf2_save<-cbind(irf2_save,irf2)
  }
  
  #nrow(irf_err[,,1]) #nrow-> sono ogni simulazione
  means1<-colMeans(irf_err[,,1])#errori della stima delle
  RMSE1<-sqrt(colMeans((irf_err[,,1])^2))
  #Irf con stima del modello di Hannan
  #length(means1)
  valori<-c()
  for (i in 1:(maxl+1)) {
    valori<-cbind(valori, mean(means1[indices[i,]]))
  }
  #valori
  RMS<-c()
  for (i in 1:(maxl+1)) {
    RMS<-cbind(RMS, mean(RMSE1[indices[i,]]))
  }
  #varianze errori
  var1<-apply(irf_err[,,1], 2, var)
  
  values<-c()
  for (i in 1:(maxl+1)) {
    values<-cbind(values, mean(var1[indices[i,]]))
  }
  #values
  #coverage
  #da rifareeee
  irf_min<-c()
  irf_max<-c()
  long<-length(irf2_save[,1])
  for (i in 1:long) {
    irf_min<-cbind(irf_min,min(irf2_save[i,]))
    irf_max<-cbind(irf_max,max(irf2_save[i,]))
  }
  irf1_val<-as.vector(irf1)
  irf_min<-as.vector(irf_min)
  irf_max<-as.vector(irf_max)
  coverage<-mean(irf1_val<irf_max&irf1_val>irf_min)
  
  results<-list(
    NV=NV,
    simulazioni=simulazioni,
    valori=valori,
    RMSE=RMS,
    values=values,
    coverage=coverage,
    differenze_irf=diff_irf,
    diff_irf=mean(diff_irf),
    irf2=irf2_save,
    varianza=mean(apply(irf2_save, 1,var)),
    risultati_varma=risultati_varma
    
    
    
    
    
  )
  return(results)
}

################################################
##########################FABLE
#############################################################
simulvarmanoidfable<-function(NV=1000, simulazioni=5){
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  risultati_varma <- vector("list", simulazioni)
  
  for (o in 1:simulazioni) {
    print(o)
    
    Y <- sim_varma(mod1, n = 1000)
    #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
    dati_series <- tibble(
      date = seq(as.Date("2000-01-01"), by = "day", length.out = 1000),
      Y1 = Y[,1],
      Y2 = Y[,2],
      Y3 = Y[,3]
    )
    dati_ts <- as_tsibble(dati_series, index = date)
    scan_gaps(dati_ts)
    dati_ts <- dati_ts %>% 
      fill_gaps()
    
    #dati_ts <- dati_ts %>%
    #  fill(Y1, Y2, .direction = "down")
    # 2. Stima del modello VARMA
    # Esempio: VARMA(1,1) -> AR=1, MA=1. d=0 indica nessuna integrazione (stazionario)
    ih_fit <-  tryCatch({dati_ts %>%
        model(
          mio_varma = VARIMA(vars(Y1, Y2, Y3) ~ pdq(2, 0, 2), identification = "kronecker_indices")
        )  %>%
        convert_fable_varma()
    }, error = function(e) {
      message(paste("Errore nella simulazione", o, "- Tentativo :", e$message))
      return(NULL) # Ritorna NULL se fallisce
    })
    
    #coefficienti<-tidy(fit)
    
    #Stima del VARMA(secondo hannan o altri processi)
    #ih_fit <- varma_from_var(Y, s, l, TRUE)
    #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
    #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
    # simte con: VARMACpp, ss_varma
    roots <- inv_roots(ih_fit, plot = F)
    while (tryCatch({!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))}, error = function(e) {
      
      return(TRUE) # Ritorna NULL se fallisce
    })) {
      Y <- sim_varma(mod1, n = 1000)
      dati_series <- tibble(
        date = seq(as.Date("2000-01-01"), by = "day", length.out = 1000),
        Y1 = Y[,1],
        Y2 = Y[,2],
        Y3 = Y[,3]
      )
      dati_ts <- as_tsibble(dati_series, index = date)
      scan_gaps(dati_ts)
      dati_ts <- dati_ts %>% 
        fill_gaps()
      
      #dati_ts <- dati_ts %>%
      #  fill(Y1, Y2, .direction = "down")
      # 2. Stima del modello VARMA
      # Esempio: VARMA(1,1) -> AR=1, MA=1. d=0 indica nessuna integrazione (stazionario)
      ih_fit <- tryCatch({dati_ts %>%
          model(
            mio_varma = VARIMA(vars(Y1, Y2, Y3) ~ pdq(2, 0, 2), identification = "kronecker_indices")
          )  %>%
          convert_fable_varma()
      }, error = function(e) {
        message(paste("Errore nella simulazione", o, "- Tentativo :", e$message))
        return(NULL) # Ritorna NULL se fallisce
      })
      roots <- inv_roots(ih_fit, plot = F)}
    
    risultati_varma[[o]]<-ih_fit
    
    irf2 <- irf(
      ih_fit,
      maxlag = maxl,
      orth = "none",
      plot = FALSE
    )
    
    
    
    psi_err<-c()
    for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
    
    diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
    
    irf_err[o,,1]<-as.numeric(psi_err)
    irf2<-as.vector(irf2)
    irf2_save<-cbind(irf2_save,irf2)
  }
  
  #nrow(irf_err[,,1]) #nrow-> sono ogni simulazione
  means1<-colMeans(irf_err[,,1])#errori della stima delle
  RMSE1<-sqrt(colMeans((irf_err[,,1])^2))
  #Irf con stima del modello di Hannan
  #length(means1)
  valori<-c()
  for (i in 1:(maxl+1)) {
    valori<-cbind(valori, mean(means1[indices[i,]]))
  }
  #valori
  RMS<-c()
  for (i in 1:(maxl+1)) {
    RMS<-cbind(RMS, mean(RMSE1[indices[i,]]))
  }
  #varianze errori
  var1<-apply(irf_err[,,1], 2, var)
  
  values<-c()
  for (i in 1:(maxl+1)) {
    values<-cbind(values, mean(var1[indices[i,]]))
  }
  #values
  #coverage
  #da rifareeee
  irf_min<-c()
  irf_max<-c()
  long<-length(irf2_save[,1])
  for (i in 1:long) {
    irf_min<-cbind(irf_min,min(irf2_save[i,]))
    irf_max<-cbind(irf_max,max(irf2_save[i,]))
  }
  irf1_val<-as.vector(irf1)
  irf_min<-as.vector(irf_min)
  irf_max<-as.vector(irf_max)
  coverage<-mean(irf1_val<irf_max&irf1_val>irf_min)
  
  results<-list(
    NV=NV,
    simulazioni=simulazioni,
    valori=valori,
    RMSE=RMS,
    values=values,
    coverage=coverage,
    differenze_irf=diff_irf,
    diff_irf=mean(diff_irf),
    irf2=irf2_save,
    varianza=mean(apply(irf2_save, 1,var)),
    risultati_varma=risultati_varma
    
    
    
    
    
  )
  return(results)
}

gen_varma_noid_simulation<- function(seed=123, maxl=7, part="varma2", NV=50, simulazioni=1000){
  if (!part %in% c("var","varma2", "varma3")) {
    
    # Se non è presente, ferma l'esecuzione e mostra un errore personalizzato
    stop(paste0("Errore: il valore '", part, "' non è valido. I valori ammessi sono: ", 
                paste(c("var","varma2", "varma3"), collapse = ", ")))
  }
  set.seed(seed)
  
  if(part=="varma2"){res<-genera_2matrici()
  ar<-as.vector(res$A1)
  ma<-as.vector(res$A2)}
  
  if(part=="varma3"){
    res<-genera_4_matrici_3x3_rango_ridotto()
    
    ar<-as.vector(cbind(res$A1,res$A2))
    ma<-as.vector(cbind(res$A3,res$A4))}
  if(part=="var"){
    phi1 <- matrix(runif(9), 3, 3)
    phi1[upper.tri(phi1, TRUE)] <- 0
    res<-c()
    res$A1<-phi1
    phi2<-phi1%*%phi1
    ar<-as.vector(phi1)
    ma<-as.vector(cbind(phi1,phi2))
  }
  #Genero le matrici dei ritardi
  
  w<-nrow(res$A1)
  s<-length(ar)/(w*w)
  l<-length(ma)/(w*w)
  sig<-matrix(0.2,w,w)
  diag(sig)<-1
  mod1<- structure(
    list(
      ar  = array(ar,c(w,w,s)),
      ma  = array(-ma, c(w,w, l)),
      cov = matrix(sig, w,w)
    ),
    class = "varma"
  )
  ##mi salvo gli iperparemtri
  
  #controllo stazionarietà
  roots <- inv_roots(mod1, plot = F)
  stationary <- all(Mod(roots$ar) < 1)
  invertible <- all(Mod(roots$ma) < 1)
  #c(stationary&invertible)# T,T
  #controllo non identificatbilità
  #!is_identified(mod1)#T
  
  while(!(stationary&invertible&!is_identified(mod1))){
    res<-genera_4_matrici_3x3_rango_ridotto()
    
    ar<-as.vector(cbind(res$A1,res$A2))
    ma<-as.vector(cbind(res$A3,res$A4))
    w<-nrow(res$A1)
    s<-length(ar)/(w*w)
    sig<-matrix(0.2,w,w)
    diag(sig)<-1
    mod1<- structure(
      list(
        ar  = array(ar,c(w,w,s)),
        ma  = array(-ma, c(w,w, s)),
        cov = matrix(sig, w,w)
      ),
      class = "varma"
    )
    #controllo la stazionarietà
    
    #controllo stazionarietà
    roots <- inv_roots(mod1, plot = F)
    stationary <- all(Mod(roots$ar) < 1)
    invertible <- all(Mod(roots$ma) < 1)
    c(stationary&invertible)# T,T
    #controllo non identificatbilità
    !is_identified(mod1)#T
  }
  
  roots <- inv_roots(mod1, plot = T)
  
  irf1 <- irf(
    mod1,
    maxlag = maxl,
    orth = "none",
    plot = FALSE
  )
  
  ####################
  irf_err<-array(0,c(simulazioni,(w*w)*(maxl+1),2))
  indices<-matrix(c(1:((maxl+1)*w*w)),ncol=w*w, byrow=TRUE)
  irf2_save<-c()
  diff_irf<-c()
  
  
  for (o in 1:simulazioni) {
    
    
    Y <- sim_varma(mod1, n = NV)
    #matplot(Y, type = "l")#potrei farlo per tutte e tre le serie
    
    #Stima del VARMA(secondo hannan o altri processi)
    ih_fit <- ih_varma1(Y, s, l, TRUE)
    #ss_fit1 <- ss_varma_kfas(Y, s, s, TRUE)#troppa potenza di calcolo?
    #tsay <- VARMACpp(Y, s, s, TRUE, details = FALSE) #mi genera Phi e Theta ma non in diverse matrici
    # simte con: VARMACpp, ss_varma
    roots <- inv_roots(ih_fit, plot = F)
    while (!(all(Mod(roots$ar) < 1) & all(Mod(roots$ma) < 1))) {
      Y <- sim_varma(mod1, n = NV)
      ih_fit <- ih_varma1(Y, s, l, TRUE)
      roots <- inv_roots(ih_fit, plot = F)}
    
    irf2 <- irf(
      ih_fit,
      maxlag = maxl,
      orth = "none",
      plot = FALSE
    )
    
    
    
    psi_err<-c()
    for(i in 1:(maxl+1)){psi_err<-cbind(psi_err, irf1[,,i]-irf2[,,i])}
    
    diff_irf<-c(diff_irf, irf_distance(irf1,irf2))
    
    irf_err[o,,1]<-as.numeric(psi_err)
    irf2<-as.vector(irf2)
    irf2_save<-cbind(irf2_save,irf2)
  }
  
  
  ####################
  results<-list(
    mod1=mod1,
    irf1=irf1,
    roots=roots,
    identified=is_identified(mod1),
    w=w,
    s=s,
    l=l,
    maxl=maxl,
    NV=NV,
    simulazioni=simulazioni,
    differenze_irf=mean(diff_irf),
    varianza=mean(apply(irf2_save, 1,var))
  )
  
  return(results)
}

simulate_varma_performances <- function(n, varma, irf_lags = 5, nsim = 1000, verbose = TRUE) {
  dv <- dim(varma)
  # Arrays to store estimates
  mat_irf <- array(NA_real_, c(nsim, prod(c(dv[1], dv[1], irf_lags)), 3))
  mat_ar  <- array(NA_real_, c(nsim, prod(dv[c(1, 1, 2)]), 3))
  mat_ma  <- array(NA_real_, c(nsim, prod(dv[c(1, 1, 3)]), 3))
  mat_cov <- array(NA_real_, c(nsim, prod(dv[1]*(dv[1]+1)/2), 3))
  # New matrix to store running times (nsim x 3)
  mat_time <- matrix(NA_real_, nrow = nsim, ncol = 3,
                     dimnames = list(NULL, c("MTS", "SSM", "HAN")))
  
  # Simulation iterations
  for (i in 1:nsim) {
    if (verbose && ((i %% 100) == 0) || i == 1) cat("Iteration n.", i, "\n")
    Y <- sim_varma(varma, n)
    
    # Produce estimates and measure running time
    # MTS::VARMA
    time_mts <- system.time({
      sink(nullfile())
      est_mts <- tryCatch(
        convert_mts_varma(MTS::VARMA(Y, dv[2], dv[3])),
        error = function(e) return(list(ar = NA, ma = NA, cov = NA))
      )
      sink()
    })["elapsed"]
    mat_time[i, 1] <- time_mts
    
    # ss_varma_fkf
    time_ssm <- system.time({
      est_ssm <- tryCatch(
        ss_varma_fkf(Y, dv[2], dv[3]),
        error = function(e) return(list(ar = NA, ma = NA, cov = NA))
      )
    })["elapsed"]
    mat_time[i, 2] <- time_ssm
    
    # ih_varma
    time_han <- system.time({
      est_han <- tryCatch(
        ih_varma(Y, dv[2], dv[3]),
        error = function(e) return(list(ar = NA, ma = NA, cov = NA))
      )
    })["elapsed"]
    mat_time[i, 3] <- time_han
    
    # Store estimation results
    if (dv[2] > 0) {
      mat_ar[i, , 1] <- est_mts$ar
      mat_ar[i, , 2] <- est_ssm$ar
      mat_ar[i, , 3] <- est_han$ar
    }
    if (dv[3] > 0) {
      mat_ma[i, , 1] <- est_mts$ma
      mat_ma[i, , 2] <- est_ssm$ma
      mat_ma[i, , 3] <- est_han$ma
    }
    mat_cov[i, , 1] <- est_mts$cov[lower.tri(est_mts$cov, diag = TRUE)]
    mat_cov[i, , 2] <- est_ssm$cov[lower.tri(est_mts$cov, diag = TRUE)]
    mat_cov[i, , 3] <- est_han$cov[lower.tri(est_mts$cov, diag = TRUE)]
    
    if (!anyNA(est_mts)) mat_irf[i, , 1] <- irf(est_mts, irf_lags, varnames = NULL)[ , , -1]
    if (!anyNA(est_ssm)) mat_irf[i, , 2] <- irf(est_ssm, irf_lags, varnames = NULL)[ , , -1]
    if (!anyNA(est_han)) mat_irf[i, , 3] <- irf(est_han, irf_lags, varnames = NULL)[ , , -1]
  }
  
  # Assign names to arrays
  mat_names <- outer(paste0("y", 1:dv[1]), paste0("y", 1:dv[1]), function(x, y) paste(x, y, sep = "_"))
  vec_names <- as.character(mat_names)
  if (dv[2] > 0) dimnames(mat_ar) <- list(NULL,
                                          paste(vec_names, rep(1:dv[2], each = length(vec_names)), sep = "_"),
                                          c("MTS", "SSM", "HAN"))
  if (dv[3] > 0) dimnames(mat_ma) <- list(NULL,
                                          paste(vec_names, rep(1:dv[3], each = length(vec_names)), sep = "_"),
                                          c("MTS", "SSM", "HAN"))
  dimnames(mat_irf) <- list(NULL,
                            paste(vec_names, rep(1:irf_lags, each = length(vec_names)), sep = "_"),
                            c("MTS", "SSM", "HAN"))
  dimnames(mat_cov) <- list(NULL,
                            mat_names[lower.tri(mat_names, diag = TRUE)],
                            c("MTS", "SSM", "HAN"))
  
  # Return results including timing matrix
  list(
    irf = mat_irf,
    ar = mat_ar,
    ma = mat_ma,
    cov = mat_cov,
    time = mat_time
  )
}

gen_varma_noid_grade<-function(n, maxl){
  res<-genera_2matrici(n)
  ar<-as.vector(res$A1)
  ma<-as.vector(res$A2)
  
  
  w<-nrow(res$A1)
  s<-length(ar)/(w*w)
  l<-length(ma)/(w*w)
  sig<-matrix(0.2,w,w)
  diag(sig)<-1
  mod1<- structure(
    list(
      ar  = array(ar,c(w,w,s)),
      ma  = array(-ma, c(w,w, l)),
      cov = matrix(sig, w,w)
    ),
    class = "varma"
  )
  
  
  #controllo stazionarietà
  roots <- inv_roots(mod1, plot = F)
  stationary <- all(Mod(roots$ar) < 1)
  invertible <- all(Mod(roots$ma) < 1)
  
  
  roots <- inv_roots(mod1, plot = T)
  
  irf1 <- irf(
    mod1,
    maxlag = maxl,
    orth = "none",
    plot = FALSE
  )
  results<-list(
    mod1=mod1,
    irf1=irf1,
    roots=roots,
    identified=is_identified(mod1),
    w=w,
    s=s,
    l=l,
    maxl=maxl
  )}
