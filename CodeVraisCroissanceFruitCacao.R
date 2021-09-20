#title: "Rapport de stage"
# author: "Sinclair Tsana"
# date: "`r Sys.Date()`"

## Import des librairies utilisés
library(pracma) # Package nécéssaire pour l'usage de la fonction integral2() 
## Import des données
dataEtats <- "/home/sinclair/Documents/ESSFAR/Projets et stages ESSFAR/Stage 2021-Mr Takam/Tsana/dataEtats.txt"
DataStadeDevelopement <- "/home/sinclair/Documents/ESSFAR/Projets et stages ESSFAR/Stage 2021-Mr Takam/Tsana/DataStadeDevelopement.txt"
Donnee_Climat_V1 <- "/home/sinclair/Documents/ESSFAR/Projets et stages ESSFAR/Stage 2021-Mr Takam/Tsana/Donnee_Climat_V1.txt"

data1 <- read.table(dataEtats, header = TRUE, sep = "")
data2 <- read.table(DataStadeDevelopement, header = TRUE, sep = "")
data3 <- read.table(Donnee_Climat_V1, header = TRUE, sep = "")

data.e <- data.frame(data1)[,-1]  # Suppression de la première colonne
data.s <- data.frame(data2)[,-1]
data.cl <- data.frame(data3)[,-1]

plv <- data.cl[,4]                # pluviométrie
Ta <- (data.cl[,5]+data.cl[,6])/2 # température ambiante
Tc <- (data.cl[,7]+data.cl[,8])/2 # température sous cacaoyère
plv[plv==0]=0.01
cl = data.frame(plv,Ta,Tc)

## 1.Fonction qui donne les numéros des individus par stade

Ik = function(k, datas){
  res=c()
  for(i in 1:nrow(datas)){
    y = datas[i,]
    if((sum(y==k)!=0) & (sum(y==k)!=1)) res[i]=i
    else res[i]=0
    }
  res
}

## 2.Statistique du nbre d'individus par stade

nb1 = length(Ik(1,data.s)[Ik(1,data.s) !=0])
nb2 = length(Ik(2,data.s)[Ik(2,data.s) !=0])
nb3 = length(Ik(3,data.s)[Ik(3,data.s) !=0])

## 3.Durée maximale discrète par stade

duree.disc = function(k,datas,L=7){
  ik = Ik(k, datas); ik1=ik[ik!=0];res=c()
  for(i in 1:length(ik1)){
    y=datas[ik1[i],]
    res[i] = (sum(y==k))*L + 1
  }
  res
}

## 4.Vecteur des poids

w0k<- function(datas){
  w01 = 1/max(duree.disc(1,datas))
  w02 = 1/max(duree.disc(2,datas))
  w03 = 1/max(duree.disc(3,datas))
  res = c(w01,w02,w03)
  res
}

## 5. Calcul de x_tild

x_tild <- function(data.cl.cur,datas=data.s,u,v,k){#data.cl.cur=Ta,Tc ou plv
  u <- ceiling(u) 
  v <- floor(v)   
  n = v-u+1
  if(u>v){stop("not defined for u>v")}
  if(v==u)
    return(0)
  else
    res <- sum(data.cl.cur[1:(n-1)]) + (1/2)*(data.cl.cur[u]+data.cl.cur[v])
  return(res*w0k(datas)[k])
}

## 6. Fonction de Hazard Cumulée

alpha.1 = rep(1,3) #Initialisation des paramètres 
mu.1 = rep(1,3)

Hk <- function(v,u,climat=cl,datas=data.s,alpha=alpha.1,mu=mu.1,k=1){
  if(v<u){
    stop("not defined for v<u")}
  else if(v==u)
    return(0)
  else{
    plv = climat[,1];Ta=climat[,2];Tc=climat[,3]
    v.cum_clim <- c(x_tild(plv,datas,u,v,k),
                    x_tild(Ta,datas,u,v,k),
                    x_tild(Tc,datas,u,v,k))
    som <- sum(alpha* (v.cum_clim^mu))
    return(som)
  }
}

## 7.Fonction de Survie

Sk <- function(v,u,climat=cl,datas=data.s,alpha=alpha.1,mu=mu.1,k=1){
  if(v>=u){
    return(exp(-Hk(v,u,climat,datas,alpha,mu,k)))
  }else{
    stop("Verifiez que v>=u")
  }
}

## 8. La densité de la durée au stade k

phi_k <- function(v,u,climat=cl,datas=data.s,alpha=alpha.1,mu=mu.1,k=1){
  Ta.tild <- x_tild(Ta,datas,u,v,k)
  Tc.tild <- x_tild(Tc,datas,u,v,k)
  plv.tild <- x_tild(plv,datas,u,v,k)
  v.cum_clim <- c(plv.tild,Ta.tild,Tc.tild)
  v.inter <- c(plv[v],Ta[v], Tc[v])
  betha <- alpha*mu*v.inter
  s <- sum(betha*v.cum_clim^(mu-1))
  return(Sk(v,u,climat,datas,alpha,mu,k)*s)
}

## 9. Statut d'un individu à la sortie du stade k

statut <- function(i,k=1,tp=20,datas=data.s, datae=data.e){
  ik = Ik(k,datas)
  if(ik[i]==0){stop("L'individu n'est pas dans ce stade")}
  else{
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) # Premier jour de i au stade k
    eikp1 <- max((1:length(y))[y==k]) # Dernier jour de i au stade k
    if(eikp1<tp){
      #Cas ou i a changé de stade
      if((datae[i,eikp1]==0) & (datas[i,eikp1] != datas[i,eikp1+1]) & (datae[i,eikp1+1]==0))
        res=1
      if((datae[i,eikp1]==0) & (datas[i,eikp1] != datas[i,eikp1+1]) & (datae[i,eikp1+1]!=0))
        res=2
      # Cas ou i n'a pas changé de stade
      if(datae[i,eikp1]!=0)
        res=3
    }
    if(eikp1==tp){
      #Ici l'individu n'a pas pu changé car il est observé au stade k à tp
      # Cas ou i est observé sain a tp
      if(datae[i,eikp1]==0) res=4
      # Cas ou i est observé attaqué à tp
      if(datae[i,eikp1]!=0) res=5
    }
  }
  res
}

## 10.Statut de tous les individus dans le stade k

statut.data.k = function(k=1,tp=20,datas=data.s, datae=data.e){
  res=c()
  ik = Ik(k,datas);ik=ik[ik!=0]
  for(i in ik)
    res[i] = statut(i,k,tp,datas, datae)
  res
}

####################### Stade 1 ###################
## 11.Calcul de la vraisemblance au stade 1
### 11.1 Calcul de la probabilité d'etre observé sain au stade suivant

fct1 = function(u,i=7,climat=cl,datas=data.s,alpha=alpha.1,mu=mu.1,k=1,L=7){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1; eikp1.c = eikp1*L + 1
    res1 = Sk(eikp1.c-L,u,climat,datas,alpha,mu,k)
    res2 = Sk(eikp1.c,u,climat,datas,alpha,mu,k)
    res = (1/L)*(res1-res2)
    res
  }
}

fct1v = Vectorize(fct1)
integrate(fct1v,1,8,i=7,datas=data.s,climat=cl,alpha=alpha.1,mu=mu.1,k=1,L=7)

### 11.2 Calcul de la probabilité d'etre observé infecté au stade suivant

fct2 = function(m,u,i=7,climat=cl,datas=data.s,alpha=alpha.1,mu=mu.1,k=1,L=7){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) # Premier jour de i au stade k
    eikp1 <- max((1:length(y))[y==k]) # dernier jour de i au stade k
    eik.c = eik*L + 1; eikp1.c = eikp1*L + 1
    res1 = Sk(eikp1.c-L,u,climat,datas,alpha,mu,k)
    res2 = Sk(m,u,climat,datas,alpha,mu,k)
    res = (1/L^2)*(res1-res2)
    res
  }
}

fct2v = Vectorize(fct2)
integral2(fct2v,xmin=15,xmax=22,ymin=1,ymax=8)
  
### 11.3 Calcul de la probabilité d'etre observé infecté au dernier jour d'observation(avant la fin de l'etude)

fct3 = function(m,u,i=7,climat=cl,datas=data.s,alpha=alpha.1,mu=mu.1,k=1){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    res = (1/L^2)*Sk(m,u,climat,datas,alpha,mu,k)
  }
  res
}

fct3v = Vectorize(fct3)
integral2(fct3v, xmin=15,xmax=22,ymin=1,ymax=8)

### 11.4 Calcul de la probabilité d'etre observé sain au dernier jour d'observation de l'étude

fct4 = function(u,i=7,climat=cl,datas=data.s,alpha=alpha.1,mu=mu.1,k=1,L=7,tp=20){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    res = (1/L)*Sk(tp*L+1,u,climat,datas,alpha,mu,k)
  }
  res
}
fct4v = Vectorize(fct4)
integrate(fctv4,lower =1 ,upper = 8)

### 11.5 Calcul de la probabilité d'etre observé attaqué au dernier jour d'observation de l'étude
tp=20;L=7
integral2(fct3v, xmin=tp*L+1-L,xmax= tp*L+1,ymin=1,ymax=8)

### 11.6 Calcul de la vraisemblance au stade 1
# On rappel qu'on ne considerera que les individus ayant été observé au moins deux fois dans le stade

Lvrais1 = function(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20){
    alpha = coefs[1:3]
    mu = coefs[4:6]
    ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
    for(i in ik1){
      y = datas[i,]
      eik <- min((1:length(y))[y==k]) 
      eikp1 <- max((1:length(y))[y==k]) 
      eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
      if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
      if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
      if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
      if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
      if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
    }
    res=res[!is.na(res)]; res[res==0]=0.0001
    resf=-sum(log(res))
    resf
}

### 11.7 Estimation des paramètres au stade 1

coefs.init=rep(1,6)
reg1=nlm(Lvrais1,coefs.init,hessian = T,climat=cl,datas=data.s,k=1,L=7,tp=20)
estimation.stade1=reg1$estimate

### 11.8 Analyse de la variabilité des paramètres du stade 1

boot.stade1 = function(B=100,coefs.init,climat=cl,datas=data.s,k=1,L=7,tp=20){
  res=matrix(0,B,6)
  for(i in 1:B){
    inter1 = sample(1:nrow(datas),replace = T)
    datas.inter=datas[inter1,]
    datae.inter=datae[inter1,]
    res[i,]=nlm(Lvrais1,datas=datas.inter,datae=datae.inter)$estimate
    write.table(res, file="Résultat1.txt")
  }
}
#################### Stade 2 #########################
## 12.Calcul de la vraisemblance au stade 2
### 12.1 Fonction qui calcule l'intégral de phi2
### 12.1.1 Densité de la durée au stade 2
param1=estimation.stade1 # param1 représente les paramètres estimés au stade 1
alpha.2 = param1[1:3] 
mu.2 = param1[4:6]

phi2=function(v,u,climat=cl,datas=data.s,alpha=alpha.2,mu=mu.2,k=2){
  phi_k(v,u,climat,datas,alpha,mu,k)
}

### 12.1.2 Calcul de l'intégral de phi2

int.phi2=function(v,u,i=7,climat=cl,datas=data.s,alpha=alpha.2,mu=mu.2,k=2,L=7,tp=20){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eik.c = eik*L + 1
    f=phi2(tp*L-L,u,climat=cl,datas=data.s,alpha=alpha.2,mu=mu.2,k=2)
    res=integrate(f,lower = eik.c-L, upper = eik.c)
  }
  resf = ifelse(res==0,0.0001,res)
  resf
}

### 12.2 Calcul de la probabilité d'etre observé sain au stade suivant

fcx1 = function(u,i=7,climat=cl,datas=data.s,alpha=alpha.2,mu=mu.2,k=2,L=7){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1; eikp1.c = eikp1*L + 1
    res1 = Sk(eikp1.c-L,u,climat,datas,alpha,mu,k)
    res2 = Sk(eikp1.c,u,climat,datas,alpha,mu,k)
    res3 = phi2(eikp1.c,u,climat,datas,alpha,mu,k=2)
    res4= int.phi2(eikp1.c,u,i,climat,datas,alpha,mu,k,tp)
    res = (1/L)*(res1-res2)*(res3/res4)
    res
  }
}
fcx1v = Vectorize(fcx1)
integrate(fcx1v,lower = 15,upper=22)

## 12.3 Calcul de la probabilité d'etre observé infecté au stade suivant

fcx2 = function(m,u,i=7,climat=cl,datas=data.s,alpha=alpha.2,mu=mu.2,k=2,L=7){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    y = datas[i,]
    eikp1 <- max((1:length(y))[y==k]) 
    eikp1.c = eikp1*L + 1
    res1 = Sk(eikp1.c-L,u,climat,datas,alpha,mu,k)
    res2 = Sk(m,u,climat,datas,alpha,mu,k)
    res3 = phi2(m,u,climat,datas,alpha,mu,k=2)
    res4= int.phi2(m,u,i,climat,datas,alpha,mu,k,tp)
    res = (1/L)*(res1-res2)*(res3/res4)
    res
  }
}

fcx2v = Vectorize(fcx2)
integral2(fcx2v,xmin=15,xmax=22,ymin=1,ymax=8)

### 12.4 Calcul de la probabilité d'etre observé infecté au dernier jour d'observation(avant la fin de l'etude)

fcx3 = function(m,u,i=7,climat=cl,datas=data.s,alpha=alpha.2,mu=mu.2,k=2,tp=20){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    res1 = Sk(m,u,climat,datas,alpha,mu,k)
    res2 = phi2(tp*L+1,u,climat,datas,alpha,mu,k=2)
    res3 = int.phi2(tp*L+1,u,i,climat,datas,alpha,mu,k,tp)
    res = (1/L)*res1*(res2/res3)
  }
  res
}
fcx3v = Vectorize(fcx3)
integral2(fcx3v, xmin=15,xmax=22,ymin=1,ymax=8)

### 12.5 Calcul de la probabilité d'etre observé sain au dernier jour d'observation de l'etude

fcx4 = function(u,i=7,climat=cl,datas=data.s,alpha=alpha.2,mu=mu.2,k=2,L=7,tp=20){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    res1 = phi2(tp*L+1,u,climat,datas,param1,k=2)
    res2= int.phi2(tp*L+1,u,i,climat,datas,alpha,mu,k,tp)
    res = (1/L)*Sk(tp*L+1,u,climat,datas,alpha,mu,k)*(res1/res2)
  }
  res
}
fcx4v = Vectorize(fcx4)
integrate(fcx4v,lower =1 ,upper = 8)

### 12.6 Calcul de la probabilité d'etre observé attaqué au dernier jour d'observation de l'etude
tp=20;L=7
integral2(fcx3v, xmin=tp*L+1-L,xmax= tp*L+1,ymin=eik.c-L,ymax=eik.c)

Lvrais2 = function(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = coefs[1:3]
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fcx1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fcx2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fcx3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fcx4v,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fcx3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

### 12.7 Estimation des paramètres au stade 2
param1=estimation.stade1
reg2=nlm(Lvrais2,param1,hessian = T,climat=cl,datas=data.s,k=2,L=7,tp=20)
estimation.stade2=reg2$estimate

### 12.8 Analyse de la variabilité des paramètres

boot.stade2 = function(B=100,param1,climat=cl,datas=data.s,k=2,L=7,tp=20){
  res=matrix(0,B,6)
  for(i in 1:B){
    inter1 = sample(1:nrow(datas),replace = T)
    datas.inter=datas[inter1,]
    datae.inter=datae[inter1,]
    res[i,]=nlm(Lvrais2,datas=datas.inter,datae=datae.inter)$estimate
    write.table(res, file="Résultat2.txt")
  }
}

#################### Stade 3 #########################
## 13.Calcul de la vraisemblance au stade 3
### 13.1 Fonction qui calcule l'intégral de phi3
### 13.1.1 Densité de la durée au stade 3
param2=estimation.stade2 # param2 représente les paramètres estimés au stade 2
alpha.3 = param2[1:3] 
mu.3 = param2[4:6]

phi3=function(v,u,climat=cl,datas=data.s,alpha=alpha.2,mu=mu.2,k=3){
  phi2(v,u,climat,datas,alpha,mu,k)
}

### 13.1.2 Calcul de l'intégral de phi3
int.phi3=function(v,u,i=7,climat=cl,datas=data.s,alpha=alpha.3,mu=mu.3,k=3,L=7,tp=20){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eik.c = eik*L + 1
    f=phi2(tp*L-L,u,climat=cl,datas=data.s,alpha=alpha.3,mu=mu.3,k)
    res=integrate(f,lower = eik.c-L, upper = eik.c)
  }
  resf = ifelse(res==0,0.0001,res)
  resf
}

### 13.2 Calcul de la probabilité d'etre observé sain au stade suivant

fcz1 = function(u,i=7,climat=cl,datas=data.s,alpha=alpha.3,mu=mu.3,k=3,L=7){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    y = datas[i,]
    eikp1 <- max((1:length(y))[y==k]) 
    eikp1.c = eikp1*L + 1
    res1 = Sk(eikp1.c-L,u,climat,datas,alpha,mu,k)
    res2 = Sk(eikp1.c,u,climat,datas,alpha,mu,k)
    res3 = phi2(eikp1.c,u,climat,datas,alpha,mu,k)
    res4= int.phi2(eikp1.c,u,i,climat,datas,alpha,mu,k,tp)
    res5 = ifelse(res4==0,0.0001,res4)
    res = (1/L)*(res1-res2)*(res3/res4)
    res
  }
}
fcz1v = Vectorize(fcz1)
integrate(fcz1v,lower = 1,upper = 8)

## 13.3 Calcul de la probabilité d'etre observé infecté au stade suivant

fcz2 = function(m,u,i=7,climat=cl,datas=data.s,alpha=alpha.3,mu=mu.3,k=3,L=7){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    y = datas[i,]
    eikp1 <- max((1:length(y))[y==k]) 
    eikp1.c = eikp1*L + 1
    res1 = Sk(eikp1.c-L,u,climat,datas,alpha,mu,k)
    res2 = Sk(m,u,climat,datas,alpha,mu,k)
    res3 = phi2(m,u,climat,datas,alpha,mu,k)
    res4= int.phi2(m,u,i,climat,datas,alpha,mu,k,tp)
    res = (1/L)*(res1-res2)*(res3/res4)
    res
  }
}

fcz2v = Vectorize(fcz2)
integral2(fcz2v,xmin=15,xmax=22,ymin=1,ymax=8)

### 13.4 Calcul de la probabilité d'etre observé infecté au dernier jour d'observation(avant la fin de l'etude)

fcz3 = function(m,u,i=7,climat=cl,datas=data.s,alpha=alpha.3,mu=mu.3,k=3){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    res1 = phi2(tp*L+1,u,climat,datas,alpha,mu,k)
    res2= int.phi2(tp*L+1,u,i,climat,datas,alpha,mu,k,tp)
    res = (1/L)*Sk(m,u,climat,datas,alpha,mu,k)*(res1/res2)
  }
  res
}
fcz3v = Vectorize(fcz3)
integral2(fcz3v, xmin=15,xmax=22,ymin=1,ymax=8)

### 13.5 Calcul de la probabilité d'etre observé sain au dernier jour d'observation de l'étude

fcz4 = function(u,i=7,climat=cl,datas=data.s,alpha=alpha.3,mu=mu.3,k=3,L=7,tp=20){
  ik=Ik(k,datas);ik=ik[ik!=0]
  if(sum(which(ik==i))==0)
    stop("L'individu n'est pas dans le stade concerné")
  else{
    res1 = phi2(tp*L+1,u,climat,datas,param1,k)
    res2= int.phi2(tp*L+1,u,i,climat,datas,alpha,mu,k,tp)
    res = (1/L)*Sk(tp*L+1,u,climat,datas,alpha,mu,k)*(res1/res2)
  }
  res
}
fcz4v = Vectorize(fcz4)
integrate(fcz4v,lower =1 ,upper = 8)

### 13.6 Calcul de la probabilité d'etre observé attaqué au dernier jour d'observation de l'étude
tp=20;L=7
integral2(fcz3v, xmin=tp*L+1-L,xmax= tp*L+1,ymin=1,ymax=8)

Lvrais3 = function(coefs,climat=cl,datas=data.s,k=3,L=7,tp=20){
  alpha = coefs[1:3]
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fcx1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fcx2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fcx3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fcx4v,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fcx3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

### 13.7 Estimation des paramètres au stade 3
param2=estimation.stade2
reg3=nlm(Lvrais3,param2,hessian = T,climat=cl,datas=data.s,k=3,L=7,tp=20)
estimation.stade3=reg3$estimate

### 13.8 Analyse de la variabilité des paramètres

boot.stade3 = function(B=100,param2,climat=cl,datas=data.s,k=3,L=7,tp=20){
  res=matrix(0,B,6)
  for(i in 1:B){
    inter1 = sample(1:nrow(datas),replace = T)
    datas.inter=datas[inter1,]
    datae.inter=datae[inter1,]
    res[i,]=nlm(Lvrais3,datas=datas.inter,datae=datae.inter)$estimate
    write.table(res, file="Résultat3.txt")
  }
}

############ Test du rapport de vraisemblance #######################
## 14. Test de l'effet des variables climatiques 
########## Tests du stade 1 ################
########## Ho: alpha1=0 contre H1: alpha1!=0 ############
Lvrais1.H0.alpha1 = function(coefs=estimation.stade1,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = c(0,coefs[2:3])
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

Lvrais1.H1 = function(coefs=estimation.stade1,climat=cl,datas=data.s,k=1,L=7,tp=20){
  Lvrais1(coefs,climat,datas,k,tp)
}

test1.alpha1 = function(niv = 0.05,coefs=estimation.stade1,climat=cl,datas=data.s,k=1,L=7,tp=20){
  L1 = Lvrais1.H0.alpha1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  L2 = Lvrais1.H1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  lambda = -2*(L1-L2)
  deg.liber = ncol(cl)-1
  dim.param = 3
  valeur.p = 1-pchisq(lambda,deg.liber)
  if(lambda>qchisq(1-niv,dim.param))
    res=paste("On rejette Ho avec un risque : ",niv)
  if(lambda<=qchisq(1-niv,dim.param))
    res=paste("On accepte Ho : ",niv)
  res
}

########## Ho: alpha2=0 contre H1: alpha2!=0 ############
Lvrais1.H0.alpha2 = function(coefs=estimation.stade1,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = c(coefs[1],0,coefs[3])
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

test1.alpha2 = function(niv = 0.05,coefs=estimation.stade1,climat=cl,datas=data.s,k=1,L=7,tp=20){
  L1 = Lvrais1.H0.alpha2(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  L2 = Lvrais1.H1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  lambda = -2*(L1-L2)
  deg.liber = ncol(cl)-1
  dim.param = 3
  valeur.p = 1-pchisq(lambda,deg.liber)
  if(lambda>qchisq(1-niv,dim.param))
    res=paste("On rejette Ho avec un risque : ",niv)
  if(lambda<=qchisq(1-niv,dim.param))
    res=paste("On accepte Ho : ",niv)
  res
}

########## Ho: alpha3=0 contre H1: alpha3!=0 ############
Lvrais1.H0.alpha3 = function(coefs=estimation.stade1,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = c(coefs[1:2],0)
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

test1.alpha3 = function(niv = 0.05,coefs=estimation.stade1,climat=cl,datas=data.s,k=1,L=7,tp=20){
  L1 = Lvrais1.H0.alpha3(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  L2 = Lvrais1.H1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  lambda = -2*(L1-L2)
  deg.liber = ncol(cl)-1
  dim.param = 3
  valeur.p = 1-pchisq(lambda,deg.liber)
  if(lambda>qchisq(1-niv,dim.param))
    res=paste("On rejette Ho avec un risque : ",niv)
  if(lambda<=qchisq(1-niv,dim.param))
    res=paste("On accepte Ho : ",niv)
  res
}

########## Tests du stade 2 ################
########## Ho: alpha1=0 contre H1: alpha1!=0 ############
Lvrais2.H0.alpha1 = function(coefs=estimation.stade2,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = c(0,coefs[2:3])
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

Lvrais2.H1 = function(coefs=estimation.stade2,climat=cl,datas=data.s,k=1,L=7,tp=20){
  Lvrais2(coefs,climat,datas,k,tp)
}

test2.alpha1 = function(niv = 0.05,coefs=estimation.stade2,climat=cl,datas=data.s,k=1,L=7,tp=20){
  L1 = Lvrais2.H0.alpha1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  L2 = Lvrais2.H1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  lambda = -2*(L1-L2)
  deg.liber = ncol(cl)-1
  dim.param = 3
  valeur.p = 1-pchisq(lambda,deg.liber)
  if(lambda>qchisq(1-niv,dim.param))
    res=paste("On rejette Ho avec un risque : ",niv)
  if(lambda<=qchisq(1-niv,dim.param))
    res=paste("On accepte Ho : ",niv)
  res
}

########## Ho: alpha2=0 contre H1: alpha2!=0 ############
Lvrais2.H0.alpha2 = function(coefs=estimation.stade2,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = c(coefs[1],0,coefs[3])
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

test2.alpha2 = function(niv = 0.05,coefs=estimation.stade2,climat=cl,datas=data.s,k=1,L=7,tp=20){
  L1 = Lvrais2.H0.alpha2(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  L2 = Lvrais2.H1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  lambda = -2*(L1-L2)
  deg.liber = ncol(cl)-1
  dim.param = 3
  valeur.p = 1-pchisq(lambda,deg.liber)
  if(lambda>qchisq(1-niv,dim.param))
    res=paste("On rejette Ho avec un risque : ",niv)
  if(lambda<=qchisq(1-niv,dim.param))
    res=paste("On accepte Ho : ",niv)
  res
}

########## Ho: alpha3=0 contre H1: alpha3!=0 ############
Lvrais2.H0.alpha3 = function(coefs=estimation.stade2,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = c(coefs[1:2],0)
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

test2.alpha3 = function(niv = 0.05,coefs=estimation.stade2,climat=cl,datas=data.s,k=1,L=7,tp=20){
  L1 = Lvrais2.H0.alpha3(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  L2 = Lvrais2.H1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  lambda = -2*(L1-L2)
  deg.liber = ncol(cl)-1
  dim.param = 3
  valeur.p = 1-pchisq(lambda,deg.liber)
  if(lambda>qchisq(1-niv,dim.param))
    res=paste("On rejette Ho avec un risque : ",niv)
  if(lambda<=qchisq(1-niv,dim.param))
    res=paste("On accepte Ho : ",niv)
  res
}

########## Tests du stade 3 ################
########## Ho: alpha1=0 contre H1: alpha1!=0 ############
Lvrais3.H0.alpha1 = function(coefs=estimation.stade3,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = c(0,coefs[2:3])
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

Lvrais3.H1 = function(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20){
  Lvrais3(coefs,climat,datas,k,tp)
}

test3.alpha1 = function(niv = 0.05,coefs=estimation.stade3,climat=cl,datas=data.s,k=1,L=7,tp=20){
  L1 = Lvrais3.H0.alpha1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  L2 = Lvrais3.H1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  lambda = -2*(L1-L2)
  deg.liber = ncol(cl)-1
  dim.param = 3
  valeur.p = 1-pchisq(lambda,deg.liber)
  if(lambda>qchisq(1-niv,dim.param))
    res=paste("On rejette Ho avec un risque : ",niv)
  if(lambda<=qchisq(1-niv,dim.param))
    res=paste("On accepte Ho : ",niv)
  res
}

########## Ho: alpha2=0 contre H1: alpha2!=0 ############
Lvrais3.H0.alpha2 = function(coefs=estimation.stade3,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = c(coefs[1],0,coefs[3])
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

test3.alpha2 = function(niv = 0.05,coefs=estimation.stade3,climat=cl,datas=data.s,k=1,L=7,tp=20){
  L1 = Lvrais3.H0.alpha2(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  L2 = Lvrais3.H1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  lambda = -2*(L1-L2)
  deg.liber = ncol(cl)-1
  dim.param = 3
  valeur.p = 1-pchisq(lambda,deg.liber)
  if(lambda>qchisq(1-niv,dim.param))
    res=paste("On rejette Ho avec un risque : ",niv)
  if(lambda<=qchisq(1-niv,dim.param))
    res=paste("On accepte Ho : ",niv)
  res
}

########## Ho: alpha3=0 contre H1: alpha3!=0 ############
Lvrais3.H0.alpha3 = function(coefs=estimation.stade3,climat=cl,datas=data.s,k=1,L=7,tp=20){
  alpha = c(coefs[1:2],0)
  mu = coefs[4:6]
  ik=Ik(k,data.s);ik1=ik[ik!=0];res=c()
  for(i in ik1){
    y = datas[i,]
    eik <- min((1:length(y))[y==k]) 
    eikp1 <- max((1:length(y))[y==k]) 
    eik.c = eik*L + 1;eikp1.c = eikp1*L + 1
    if(statut(i,k)==1) res[i]=integrate(fct1v,eik-L,eik,i)
    if(statut(i,k)==2) res[i]=integral2(fct2v,xmin=eikp1.c-L,xmax=eikp1,ymin=eikp.c-L,ymax=eikp.c)
    if(statut(i,k)==3) res[i]=integral2(fct3v,xmin=eikp1.c-L,xmax=eikp1.c,ymin=eik.c-L,ymax=eik.c)
    if(statut(i,k)==4) res[i]=integrate(fctv4,lower =eik.c-L,upper= eik.c)
    if(statut(i,k)==5) res[i]=integral2(fct3v,xmin=tp*L+1-L,xmax=tp*L+1,ymin=eik.c-L,ymax=eik.c)
  }
  res=res[!is.na(res)]; res[res==0]=0.0001
  resf=-sum(log(res))
  resf
}

test3.alpha3 = function(niv = 0.05,coefs=estimation.stade3,climat=cl,datas=data.s,k=1,L=7,tp=20){
  L1 = Lvrais3.H0.alpha3(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  L2 = Lvrais3.H1(coefs,climat=cl,datas=data.s,k=1,L=7,tp=20)
  lambda = -2*(L1-L2)
  deg.liber = ncol(cl)-1
  dim.param = 3
  valeur.p = 1-pchisq(lambda,deg.liber)
  if(lambda>qchisq(1-niv,dim.param))
    res=paste("On rejette Ho avec un risque : ",niv)
  if(lambda<=qchisq(1-niv,dim.param))
    res=paste("On accepte Ho : ",niv)
  res
}

################ Fin des programmes ##############