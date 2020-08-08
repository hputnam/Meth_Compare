
all.subsets.gam <-
function(y,x.smooth,x.parametric=NULL,family=binomial(link=logit),
	maxp=5,select='all',delta=7,rank=10,...){

options(warn=-1)

#binary function copied from wle library
binary<-function(x,dim){
    if(x==0){
        pos<-1
	    }
    else{
        pos<-floor(log(x,2))+1
	    }
    if(!missing(dim)){
        if(pos<=dim){
           pos<-dim
        	}
        else{
           warning("the value of `dim` is too small")
        	}
    	}
    bin<-rep(0,pos)
    dicotomy<-rep(FALSE,pos)
    for(i in pos:1){
        bin[i]<-floor(x/2^(i-1))
        dicotomy[i]<-bin[i]==1
        x<-x-((2^(i-1))*bin[i])
	    }
return(dicotomy)
}

if(is.null(x.parametric)){
	x<-x.smooth
	}
else{
	x<-cbind(x.smooth,x.parametric)
	}
	
p<-ncol(x)
m<-(2^p)-1
N<-nrow(x)
terms<-NULL
model<-NULL
id<-NULL
AIC<-NULL
K<-NULL
AICc<-NULL
adjR2<-NULL
D2<-NULL

for(i in 1:m){
	bin<-binary(i,p)
	if(sum(bin)<=maxp){
		x.names<-names(x)[bin]

		if(is.null(x.parametric)){
			x.subset<-subset(x,select=x.names)
			formula<-as.formula(paste('y~',paste(paste("s(",names(x.subset),")",sep=""),
				sep="",collapse="+"),sep=""))			
			}

		else{
			sv.smooth<-x.names %in% names(x.smooth)
			x.smooth.names<-x.names[sv.smooth]
			x.smooth.subset<-subset(x,select=x.smooth.names)
			sv.parametric<-x.names %in% names(x.parametric)
			x.parametric.names<-x.names[sv.parametric]
			x.parametric.subset<-subset(x,select=x.parametric.names)
			x.subset<-cbind(x.smooth.subset,x.parametric.subset)
			
			if(sum(sv.smooth)==0){
				formula<-as.formula(paste('y~',paste(x.parametric.names,
					sep='',collapse='+'),sep=''))
				}
			else if(sum(sv.parametric)==0){
				formula<-as.formula(paste('y~',paste(paste('s(',x.smooth.names,')',sep=''),
					sep='',collapse='+'),sep=''))
				}
			else{
				formula<-as.formula(paste('y~',paste(c(x.parametric.names,
					paste('s(',x.smooth.names,')',sep='')),
					sep='',collapse='+'),sep=''))
				}
			}
		
		model.gam<-gam(formula,family=family,data=x.subset,...)

		#save terms
		temp<-as.data.frame(x.names)
		names(temp)<-'Term'
		temp$id<-i
		terms<-rbind(terms,temp)

		#save model aic stats
		model<-rbind(model,paste(x.names,collapse='+'))
		id<-rbind(id,i)
		AIC<-round(model.gam$aic,3)
		K<-model.gam$rank
		AICc<-rbind(AICc,round(AIC+(2*K*(K+1)/(N-K-1)),3))
		
		#save R-sq(adj) and deviance explained
		adjr2<-summary(model.gam)$r.sq
		adjR2<-rbind(adjR2,round(adjr2,3))
		d2<-summary(model.gam)$dev.expl
		D2<-rbind(D2,round(d2,3))
		}
	}

#create model AIC table
AIC.stats<-as.data.frame(cbind(id,model,adjR2,D2,as.data.frame(AICc)))
names(AIC.stats)<-c('id','Model','adjR2','D2','AICc')
AIC.stats<-AIC.stats[order(AIC.stats$AICc),]
AIC.stats$Delta<-AIC.stats$AICc-min(AIC.stats$AICc)
wgt<-exp(-0.5*AIC.stats$Delta)
AIC.stats$Wgt<-round(wgt/sum(wgt),3)
AIC.stats$Rank<-seq(1,nrow(AIC.stats))
if(select=='delta'){
	AIC.stats<-AIC.stats[AIC.stats$Delta<=delta,]
	}
else if(select=='rank'){
	AIC.stats<-AIC.stats[AIC.stats$Rank<=rank,]
	}

#create coefficients table
sv<-NULL
for(i in 1:nrow(terms)){	
	sv<-c(sv,any(terms$id[i]==AIC.stats$id))
	}
terms<-terms[sv,]
terms<-terms[order(terms$Term),]

#create variable importance table
z<-merge(terms,AIC.stats,by='id')
N<-aggregate(is.finite(z$Term),list(z$Term),sum)
names(N)<-c('Term','N')
AICc<-aggregate(z$Wgt,list(z$Term),sum)
names(AICc)<-c('Term','AICc')
varimp<-merge(N,AICc,by='Term')

#final clean up of tables
AIC.stats<-subset(AIC.stats,select=-id)

#save results to list
return(list(modelAIC=AIC.stats,variable.importance=varimp))

}
all.subsets.glm <-
function(y,x1=NULL,x2=NULL,x2.linear=FALSE,x.force=NULL,
  family=binomial(link=logit),offset=NULL,weights=NULL,maxp=5,
  select='AIC',sensitivity=.95,delta.AIC=7,delta.FP=.1,delta.Kappa=.1,
  rank=10,varimp=FALSE,gof='logLik',coef.table=FALSE,model.ave=FALSE,...){

owarn<-options("warn")
on.exit(options(warn=owarn$warn))
options(warn=-1)

#check on valid model parameterization
if(x2.linear==TRUE & is.null(x2)) stop('if x2.linear=TRUE, then x2 cannot be null')

################################################################################
#internal functions
################################################################################

#binary function (copied from wle library)
binary<-function(x,dim){
    if(x==0){
      pos<-1
	    }
    else{
      pos<-floor(log(x,2))+1
	    }
    if(!missing(dim)){
      if(pos<=dim){
        pos<-dim
        }
      else{
        warning("the value of `dim` is too small")
        }
      }
    bin<-rep(0,pos)
    dicotomy<-rep(FALSE,pos)
    for(i in pos:1){
      bin[i]<-floor(x/2^(i-1))
      dicotomy[i]<-bin[i]==1
      x<-x-((2^(i-1))*bin[i])
    	}
return(dicotomy)
}

#Kappa function
cohen.kappa<-function(y){
N<-sum(y)
ccr<-sum(diag(y))/sum(y)
p<-apply(y,1,sum)/N
q<-apply(y,2,sum)/N
num<-ccr-sum(p*q)
den<-1-sum(p*q)
k<-num/den
k[k<0]<-0
return(k)		
}

#Kappa optimization function
kappa.opt<-function(threshold){
obs<-out[,1]>0
if(threshold==0){
	pred<-out[,2]>=threshold
	}
else {
	pred<-out[,2]>threshold
	}
temp<-c(sum(pred&obs),sum(!pred&obs),sum(pred&!obs),sum(!pred&!obs))
temp<-as.table(matrix(temp,nrow=2))
cohen.kappa(temp)
}
	
################################################################################

#add quadratic terms to x2 and define svars, the selection variable vector
if(x2.linear==TRUE){ 
  svars<-c(names(x1),names(x2)) #svar does not include quadratic terms
  svars.sq<-paste(svars,'.sq',sep='') #all quadratic terms
  svars.sq.x2<-paste(names(x2),'.sq',sep='') #x2 quadratic terms
  }

if(!is.null(x2)){
  nvars<-ncol(x2)
  for(i in 1:nvars){
    temp<-x2[,i]^2
    x2<-cbind(x2,temp)
    names(x2)[nvars+i]<-paste(names(x2[i]),'.sq',sep='')
    }
  }

if(x2.linear==FALSE){
  svars<-c(names(x1),names(x2)) #svar includes quadratic terms
  }

p<-length(svars)

#define fvars, the selection variable vector for the forced terms
fvars<-names(x.force)
if(is.null(x.force)) q<-0 
else q<-ncol(x.force)

#calculate the final dimension of the storage objects
ns<-0 #number of subsets we will consider (if maxp == p this is equal to m)
nc<-0 #number of coefficients

f<-ifelse(x2.linear,2,1)
for(i in 1:min(maxp,p)){
	a<-choose(p,i)
	ns<-ns+a
	nc<-nc+a*(i*f+1+q)
  }

#create parameters and pre-allocate storage objects
m<-(2^p)-1  #number of subsets (disregarding maxp)
if(is.null(x1)) N<-nrow(x2) #number of observations
else N<-nrow(x1)

coefficients<-data.frame(character(nc),numeric(nc),numeric(nc),numeric(nc),
  numeric(nc),numeric(nc),stringsAsFactors=FALSE)
names(coefficients)<-c("Term","Estimate","Std. Error","z value","Pr(>|z|)","id")

model<-character(ns)
id<-numeric(ns)
AICc<-numeric(ns)
D2<-numeric(ns)
cutpoint<-numeric(ns)
if(select=='commission') FP<-numeric(ns)
if(select=='Kappa') Kappa<-numeric(ns)

#assemble all the data in one dataframe
if(is.null(x1) & !is.null(x2)) x<-x2
if(!is.null(x1) & is.null(x2)) x<-x1
if(!is.null(x1) & !is.null(x2)) x<-cbind(x1,x2)

if(is.null(x.force)) d<-cbind(data.frame(y=y),x) 
else d<-cbind(data.frame(y=y),x,x.force) 

#loop thru all subsets models
oi<-0 # output index
cr<-1 # the index of the first empty Coefficient table Row 
for(i in 1:m){
	bin<-binary(i,p)
	if(sum(bin)> maxp) next
	oi<-oi+1 
	
	#create a model formula	
	if(x2.linear==TRUE){
    temp<-svars.sq[bin]
    terms.sq<-temp[svars.sq[bin] %in% svars.sq.x2]
		terms<-c(svars[bin],terms.sq,fvars)
	  }
  else{
	 	terms<-c(svars[bin],fvars)
	  }
	formula<-as.formula(paste("y~",paste(terms,collapse="+"),sep=""))

  #fit model
  model.glm<-glm(formula,family=family,offset=offset,weights=weights,data=d)

	#save coefficient stats 
	if(varimp==TRUE | coef.table==TRUE){
		temp<-as.data.frame(summary(model.glm)$coefficients)
		Term<-row.names(temp)
		temp<-cbind(Term,temp,stringsAsFactors=FALSE)
		temp$id<-i
		coefficients[cr:(cr+nrow(temp)-1), ]<-temp[,]
		cr<-cr+nrow(temp)
  	}

	#save model aic stats (always)
	model[oi]<-paste(terms,collapse='+')
	id[oi]<-i
	AIC<-round(model.glm$aic,3)
	K<-model.glm$rank
	AICc[oi]<-round(AIC+(2*K*(K+1)/(N-K-1)),3)
	d2<-(model.glm$null.deviance-model.glm$deviance)/model.glm$null.deviance
	D2[oi]<-round(d2,3)

	#compute errors of commission
	if(select=='commission'){
		if(length(unique(y))>2) stop('y must be binary (0,1) for select=commission')
		out<-as.data.frame(cbind(model.glm$y,model.glm$fitted.values))
		names(out)<-c('observed','fitted')
		cut<-quantile(out$fitted[out$observed==1],prob=1-sensitivity)
		cutpoint[oi]<-cut
		FP[oi]<-round(nrow(out[out$fitted>=cut & out$observed==0,])/
			nrow(out[out$observed==0,]),3)
  	}

	#compute Kappa
	if(select=='Kappa'){
		if(length(unique(y))>2) stop('y must be binary (0,1) response for select=Kappa')
		out<-as.data.frame(cbind(model.glm$y,model.glm$fitted.values))
		names(out)<-c('observed','fitted')
		temp<-optimize(kappa.opt,interval=c(min(out$fitted),max(out$fitted)),maximum=TRUE)
		cutpoint[oi]<-temp$maximum
		Kappa[oi]<-round(temp$objective,3)
    }	
} # end for loop

#create model stats table
if(select=='commission'){
	model.stats<-data.frame(id,model,D2,AICc,cutpoint,FP)
	model.stats$deltaAICc<-round(model.stats$AICc-min(model.stats$AICc),3)
	wgtAICc<-exp(-0.5*model.stats$deltaAICc)
	model.stats$wgtAICc<-round(wgtAICc/sum(wgtAICc),3)
	model.stats$deltaFP<-model.stats$FP-min(model.stats$FP)
	model.stats<-model.stats[order(model.stats$FP),]
	model.stats$rank<-seq(1,nrow(model.stats))
	model.stats<-model.stats[model.stats$deltaFP<=delta.FP & model.stats$rank<=rank,]
	model.stats<-model.stats[,c(1,2,3,4,7,8,5,6,9,10)]
  }
else if(select=='Kappa'){
	model.stats<-data.frame(id,model,D2,AICc,cutpoint,Kappa)
	model.stats$deltaAICc<-round(model.stats$AICc-min(model.stats$AICc),3)
	wgtAICc<-exp(-0.5*model.stats$deltaAICc)
	model.stats$wgtAICc<-round(wgtAICc/sum(wgtAICc),3)
	model.stats$deltaKappa<-max(model.stats$Kappa)-model.stats$Kappa
	model.stats<-model.stats[order(model.stats$Kappa,decreasing=TRUE),]
	model.stats$rank<-seq(1,nrow(model.stats))
	model.stats<-model.stats[model.stats$deltaKappa<=delta.Kappa & model.stats$rank<=rank,]
	model.stats<-model.stats[,c(1,2,3,4,7,8,5,6,9,10)]
  } 
else{
	model.stats<-data.frame(id,model,D2,AICc)
	model.stats$deltaAICc<-model.stats$AICc-min(model.stats$AICc)
	wgtAICc<-exp(-0.5*model.stats$deltaAICc)
	model.stats$wgtAICc<-round(wgtAICc/sum(wgtAICc),3)
	model.stats<-model.stats[order(model.stats$AICc),]
	model.stats$rank<-seq(1,nrow(model.stats))
	model.stats<-model.stats[model.stats$deltaAICc<=delta.AIC & model.stats$rank<=rank,]
	rownames(model.stats) <- 1:nrow(model.stats)
  }

#create coefficients table
if(coef.table==TRUE){
	sv<-coefficients$id %in% model.stats$id
	coefficients<-coefficients[sv,]
	coefficients<-coefficients[order(coefficients$Term),]
	#coefficients<-coefficients[coefficients$Term!='(Intercept)',]
	coefficients<-coefficients[,c(6,1:5)]
  rownames(coefficients)<-1:nrow(coefficients)
	}

#create variable importance table
if(varimp==TRUE){
	z<-merge(coefficients,model.stats,by='id')
  z$Term<-gsub('.sq','',x=z$Term)
  z<-z[!duplicated(z[,1:2]),]
	N<-aggregate(rep(1, nrow(z)),list(z$Term),sum)
	names(N)<-c('Term','N')
	AICc<-aggregate(z$wgtAICc,list(z$Term),sum)
	names(AICc)<-c('Term','AICc')
	varimport.AICc<-merge(N,AICc,by='Term')
	if(p<=12){
		require(hier.part)
		if(!is.null(x.force)) x<-cbind(x,x.force) 
		HP<-round(hier.part(y,x,family=family,gof=gof)$I.perc/100,3)
		HP$Term<-row.names(HP)
		names(HP)<-c('HP','Term')
		varimport.HP<-HP[,c(2,1)]
		}
	}

#create model averaging table
if(model.ave==TRUE){
  require(MuMIn)
  form.list<-list(nrow(model.stats))
  
  for(i in 1:nrow(model.stats)){
    glm.form<-as.formula(paste('y~',model.stats$model[i],sep=''))
    form.list[[i]]<-glm(glm.form,data=d,family='binomial')
    }
  
  temp<-summary(model.avg(form.list))
  subset.model.ave.results<-as.data.frame(cbind(temp$term.names,temp$avg.model))
  names(subset.model.ave.results)[1]<-'Term'
  row.names(subset.model.ave.results)<-NULL
  full.model.ave.results<-as.data.frame(cbind(temp$term.names,temp$coef.shrinkage))
  names(full.model.ave.results)<-c('Term','Estimate')
  row.names(full.model.ave.results)<-NULL
  }

#save results to list
if(coef.table==TRUE & varimp==TRUE & model.ave==TRUE)
	return(list(model.statistics=model.stats,coefficients=coefficients,
	variable.importance.AICc=varimport.AICc,variable.importance.HP=varimport.HP,
  subset.averaged.coefficients=subset.model.ave.results,
  full.averaged.coefficients=full.model.ave.results))
	
if(coef.table==TRUE & varimp==FALSE & model.ave==TRUE)
	return(list(model.statistics=model.stats,coefficients=coefficients,
  subset.averaged.coefficients=subset.model.ave.results,
  full.averaged.coefficients=full.model.ave.results))

if(coef.table==FALSE & varimp==TRUE & model.ave==TRUE)
	return(list(model.statistics=model.stats,
  variable.importance.AICc=varimport.AICc,variable.importance.HP=varimport.HP,
  subset.averaged.coefficients=subset.model.ave.results,
  full.averaged.coefficients=full.model.ave.results))

if(coef.table==FALSE & varimp==FALSE & model.ave==TRUE)
  return(list(model.statistics=model.stats,
  subset.averaged.coefficients=subset.model.ave.results,
  full.averaged.coefficients=full.model.ave.results))

if(coef.table==TRUE & varimp==TRUE & model.ave==FALSE)
  return(list(model.statistics=model.stats,coefficients=coefficients,
  variable.importance.AICc=varimport.AICc,variable.importance.HP=varimport.HP))

if(coef.table==TRUE & varimp==FALSE & model.ave==FALSE)
  return(list(model.statistics=model.stats,coefficients=coefficients))

if(coef.table==FALSE & varimp==TRUE & model.ave==FALSE)
  return(list(model.statistics=model.stats,
  variable.importance.AICc=varimport.AICc,variable.importance.HP=varimport.HP))

if(coef.table==FALSE & varimp==FALSE & model.ave==FALSE)
  return(list(model.statistics=model.stats))
	
}
beta.plot <-
function(shape1=.5,shape2=.5){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(shape1)*length(shape2)
k<-0
for(i in shape1){
	for(j in shape2){
	  k<-k+1
		curve(pbeta(x,shape1=i,shape2=j),0,1,
			xlab='Z',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
			mtext(paste('(a=',i,', b=',j,')',sep=''),side=3,col='red')
		curve(dbeta(x,shape1=i,shape2=j),0,1,
			xlab='Z',ylab='Probability density',main='Probability Density')
			mtext(paste('(a=',i,', b=',j,')',sep=''),side=3,col='red')
		curve(qbeta(x,shape1=i,shape2=j),0,1,
			xlab='Cumulative probability (quantile)',ylab='Z',main='Quantiles')
			mtext(paste('(a=',i,', b=',j,')',sep=''),side=3,col='red')
		y<-rbeta(1000,shape1=i,shape2=j)
		hist(y,xlab='Z',ylab='Frequency',main='Random Observations')
			mtext(paste('(a=',i,', b=',j,')',sep=''),side=3,col='red')
	  if(k<n.plots) readline('press return for next plot')
		}
	}
}
binom.plot <-
function(size=10,prob=.5){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(size)*length(prob)
k<-0
for(i in size){
	n<-seq(0,i)
	for(j in prob){
	  k<-k+1
		barplot(pbinom(n,size=i,prob=j),names=as.character(n),
			xlab='# Successes (k)',ylab='Cumulative probability (quantile)',
      main='Cumulative Probability')
			mtext(paste('(size=',i,', prob=',j,')',sep=''),side=3,col='red')
		barplot(dbinom(n,size=i,prob=j),names=as.character(n),
			xlab='# Successes (k)',ylab='Probability',main='Probability Mass')
			mtext(paste('(size=',i,', prob=',j,')',sep=''),side=3,col='red')
		barplot(qbinom(seq(0,1,.05),size=i,prob=j),names=seq(0,1,.05),
			xlab='Cumulative probability (quantile)',ylab='# Successes (k)',
      main='Quantiles')
			mtext(paste('(size=',i,', prob=',j,')',sep=''),side=3,col='red')
		y<-rbinom(1000,size=i,prob=j)
    y.table<-as.data.frame(table(y))
    y.levels<-as.data.frame(n)
    names(y.levels)<-'y'
    y.table<-merge(y.levels,y.table,all.x=TRUE)
    y.table$Freq[is.na(y.table$Freq)]<-0
		barplot(y.table$Freq,names=seq(0,i),xlab='# Successes (k)',ylab='Frequency',
      main='Random Observations',col='gray')
			mtext(paste('(size=',i,', prob=',j,')',sep=''),side=3,col='red')
		if(k<n.plots) readline('press return for next plot')
		}
	}
}
box.plots <-
function(x,var='',by='',save.plot=FALSE,
	col='blue',las=1,...){

oldpar<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}
	
if(by==''){ #box-and-whisker w/o groups
	par(mar=c(5,5,4,2)) #graphics settings
	for(i in 1:ncol(y)){ #loop thru variables
		boxplot(y[i],ylab=names(y[i]),col=col,las=las,
		main=paste('Box-and-Whisker Plot of',names(y[i]),sep=' '),...)
		if(save.plot==TRUE){
			dev.print(jpeg,file=paste('box.',names(y[i]),'.jpg',sep=''),width=600,height=800)
			} #end save
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end bw w/o groups

else{ #box-and-whisker w/ groups
	n<-by.names(x,by) #create by variable
	y<-cbind(n,y) #bind with selected variables
	par(mar=c(5,5,4,2)) #graphics settings
	for(i in 3:ncol(y)){ #loop thru variables
		boxplot(y[[i]]~y[[2]],col=col,las=las,
		ylab=names(y[i]),main=paste('Box-and-Whisker Plot of',names(y[i]),sep=' '),...)
		if(save.plot==TRUE){
			dev.print(jpeg,file=paste('box.',names(y[i]),'.jpg',sep=''),width=600,height=800)
			} #end save
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end bw w/groups
par(oldpar)
}
by.names <-
function(infile,by=names(infile)){

x <- unique(infile[,by, drop=FALSE])
z <- as.character(x[,1])
if(1 < length(by))
	for(i in 2:length(by)) {
		z <- paste(z,as.character(x[,i]), sep='.')
	}
z <- data.frame(z)
names(z) <- '..key'
t <- paste(by,collapse='.')
z <- cbind(z,x,..id = 1:dim(z)[1])
infile <- cbind(1:dim(infile)[1],infile); names(infile)[1] <- '..seq'
z <- merge(infile[,c(by,'..seq')],z,by=by)
z <- z[sort(z$..seq,index.return=TRUE)$ix,]
row.names(z) <- z$..seq
z <- z[,c('..id','..key')]
names(z) <- c('id',t)
return(z)
}
chisq.plot <-
function(df=1,ncp=0,xlim=c(0,3)){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(df)
k<-0
for(i in df){
  k<-k+1
	curve(pchisq(x,df=i,ncp=ncp),xlim[1],xlim[2],
		xlab='Z',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
		mtext(paste('(df=',i,')',sep=''),side=3,col='red')
	curve(dchisq(x,df=i,ncp=ncp),xlim[1],xlim[2],
		xlab='Z',ylab='Probability density',main='Probability Density')
		mtext(paste('(df=',i,')',sep=''),side=3,col='red')
	curve(qchisq(x,df=i,ncp=ncp),0,1,
		xlab='Cumulative probability (quantile)',ylab='Z',main='Quantiles')
		mtext(paste('(df=',i,')',sep=''),side=3,col='red')
	y<-rchisq(1000,df=i,ncp=ncp)
	hist(y,xlab='Z',ylab='Frequency',main='Random Observations')
		mtext(paste('(df=',i,')',sep=''),side=3,col='red')
	if(k<n.plots) readline('press return for next plot')
	}
}
ci.lines <-
function(model){
xm<-mean(model[[12]][,2])
n<-length(model[[12]][,2])
ssx<-sum(model[[12]][,2]^2)-sum(model[[12]][,2])^2/n
s.t<-qt(0.975,(n-2))
xv<-seq(min(model[[12]][,2]),max(model[[12]][,2]),(max(model[[12]][,2])-min(model[[12]][,2]))/100)
yv<-coef(model)[1]+coef(model)[2]*xv
se<-sqrt(summary(model)[[6]]^2*(1/n+(xv-xm)^2/ssx))
ci<-s.t*se
uyv<-yv+ci
lyv<-yv-ci
lines(xv,uyv,lty=2)
lines(xv,lyv,lty=2)
}
class.monte <-
function(y,grouping='',prop=.5,type='qda',perm=1000,prior='',...){

y<-as.data.frame(y)
if(!grouping==''){
	grp<-as.factor(grouping)
	}
else{
	grp<-as.factor(y[,1]) #groups assumed to be in first column
	y<-y[,-1]
	}

G<-length(levels(grp))
z<-matrix(0,perm,G+2) #create blank matrix

for(i in 1:perm){
	ran.split(y,grp,prop)
	if(type=='qda'){
		y.da<-qda(training.data,grouping=training.grp,prior=prior,...)
		}
	else if(type=='lda'){
		y.da<-lda(training.data,grouping=training.grp,prior=prior,...)
		}
	y.pred<-predict(y.da,newdata=validate.data)
	y.table<-table(validate.grp,y.pred$class)
	for(j in 1:G){
		z[i,j]<-y.table[j,j]/sum(y.table[j,])
		}
		z[i,G+1]<-sum(diag(y.table))/sum(y.table)
		if(!prior==''){
			z[i,G+2]<-tau(y.table,prior=prior)
			}		
		else{
			z[i,G+2]<-cohen.kappa(y.table)
			}
	}

z2<-as.data.frame(apply(z,2,function(x){ #calculate stats
	z2<-c(quantile(x,probs=c(0,.05,.5,.95,1),na.rm=TRUE),
    mean(x,na.rm=TRUE))
    names(z2)<-c('minimum','5th percentile','median','95th percentile','maximum','mean')
	z2<-round(z2,3) #round elements to 3 decimal places
	}))

if(!prior==''){
	colnames(z2)<-c(levels(grp),'Total','Tau')
	}
else{
	colnames(z2)<-c(levels(grp),'Total','Kappa')
	}

cat('Split-sample Cross-Validation Classification Summary\n')
cat('Correct Classification Rate:\n')

return(z2)	
}
clus.composite <-
function(x,grp){

groups<-as.factor(grp)
groups<-levels(groups)

#compute mean by cluster
x.by<-by(x,grp,mean)
	
#create composite clusters
z<-matrix(0,length(x),length(groups)) #create blank matrix
for(i in 1:length(groups)){
	z[,i]<-x.by[[i]]
	}
z<-as.data.frame(z)
rownames(z)<-colnames(x)
colnames(z)<-groups
z<-t(z)
return(z)
}
clus.stats <-
function(x,grp){

cv<<-function(x,na.rm) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100 
x<-as.data.frame(x)
groups<-as.factor(grp)
groups<-levels(groups)

#compute kruskal-wallis rank sum test p-value
p.value<-rep(0,length(x))
for(i in 1:length(x)){
	p<-kruskal.test(x[,i],grp)
	p.value[i]<-p$p.value
	}
p.value<-round(p.value,3)

#compute number of obs. per cluster
cl.nobs<-rep(0,length(groups))
for(i in 1:length(groups)){
	cl.nobs[i]<-sum(grp==i)
	}

#compute stats by cluster
x.mean<-by(x,grp,mean)
x.cv<-by(x,grp,cv)

#create summary table
cl.mean<-matrix(0,length(x),length(groups)) #create blank matrix
cl.cv<-matrix(0,length(x),length(groups)) #create blank matrix
for(i in 1:length(groups)){
	cl.mean[,i]<-round(x.mean[[i]],3)
	cl.cv[,i]<-round(x.cv[[i]],0)
	}
cl.mean<-as.data.frame(cl.mean)
rownames(cl.mean)<-colnames(x)
colnames(cl.mean)<-groups
cl.mean<-cbind(cl.mean,p.value)

cl.cv<-as.data.frame(cl.cv)
rownames(cl.cv)<-colnames(x)
colnames(cl.cv)<-groups

z<-list(cl.nobs,cl.mean,cl.cv)
names(z)<-c('Cluster.nobs','Cluster.mean','Cluster.cv')
return(z)
}
cohen.kappa <-
function(y){

N<-sum(y)
ccr<-sum(diag(y))/sum(y)

p<-apply(y,1,sum)/N
q<-apply(y,2,sum)/N
num<-ccr-sum(p*q)
den<-1-sum(p*q)
kappa<-num/den
kappa[kappa<0]<-0
return(kappa)		
}
col.summary <-
function(x,var=NULL,by=NULL,outfile=NULL,...){

options(warn=-1)

#statistical functions
nobs<-function(x) length(x)
cv<-function(x,na.rm) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100 
zeros<-function(x,na.rm) sum(x==0,na.rm=TRUE)
pct.zeros<-function(x,na.rm) sum(x==0,na.rm=TRUE)/length(x)*100
nobs.missing<-function(x,na.rm) sum(is.na(x))
pct.missing<-function(x,na.rm) sum(is.na(x))/length(x)*100 
se<-function(x,na.rm) sd(x,na.rm=TRUE)/sqrt(length(x)-sum(is.na(x)))
se.ratio<-function(x,na.rm) se(x)/mean(x,na.rm=TRUE)*100
richness<-function(x,na.rm) nobs(x)-zeros(x)-nobs.missing(x)

#select subset of variables
if(!is.null(var)){
	y<-subset(x,select=eval(parse(text=var)))
	}
else{y<-x}

#summary table w/o groups
if(is.null(by)){ 
	z<-data.frame(apply(y,2,function(x){ #calculate stats
		z<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
		    mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sum(x,na.rm=TRUE),
		    sd(x,na.rm=TRUE),cv(x),zeros(x),pct.zeros(x),nobs.missing(x),
			pct.missing(x),se(x),se.ratio(x),richness(x))
	    names(z)<-c('nobs','min','max','mean',
		    'median','sum','sd','cv','zeros','pct.zeros',
		    'nobs.missing','pct.missing','se','se.ratio',
		    'richness') #create col names
		z<-round(z,3) #round elements to 3 decimal places
		}))
	if(!is.null(outfile)){ #save results to output file
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		}
	}

#summary table by groups
else{ 
	n<-by.names(x,by) #create by variable
	m<-levels(n[,2]) #create object with group levels
	z<-vector('list',length(m))	#create list object for output
	names(z)<-m #assign names to list components
	for(i in 1:length(m)){ #loop thru by groups
		y1<-y[n[,2]==m[i],] #select records within group
		zz<-data.frame(apply(y1,2,function(x){ #calculate stats
			zz<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
			    mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sum(x,na.rm=TRUE),
			    sd(x,na.rm=TRUE),cv(x),zeros(x),pct.zeros(x),nobs.missing(x),
				pct.missing(x),se(x),se.ratio(x),richness(x))
		    names(zz)<-c('nobs','min','max','mean',
			    'median','sum','sd','cv','zeros','pct.zeros',
			    'nobs.missing','pct.missing','se','se.ratio',
			    'richness') #create col names
			zz<-round(zz,3) #round elements to 3 decimal places
			}))
		z[[i]]<-zz
		if(!is.null(outfile)){ #save results to output file
			write(m[i],file=paste(outfile,'.csv',sep=''),append=TRUE)
			write.table(zz,file=paste(outfile,'.csv',sep=''),quote=FALSE,append=TRUE,sep=',')
			}
		} #end loop thru groups
	} #end summary table w/ groups

return(z)
}
commission.mc <-
function(fit,test=FALSE,reps=1000,sensitivity=.95,plot=TRUE,...){

owarn <- options("warn")
on.exit(options(warn=owarn$warn))
options(warn=-1)

##################################################################################
#functions
##################################################################################

#function for computing false positive (commission) error rate on glm object
FP<-function(fit,sensitivity=.95){
	temp<-as.data.frame(cbind(fit$y,fit$fitted.values))
	names(temp)<-c('observed','fitted')
	cutpoint<-quantile(temp$fitted[temp$observed==1],prob=1-sensitivity)
	FP<-(nrow(temp[temp$fitted>=cutpoint & temp$observed==0,])/
		nrow(temp[temp$observed==0,]))
	z<-list('cutpoint'=cutpoint,'FP'=FP)
	return(z)
	}
##################################################################################

#compute observed FP and cutpoint
temp<-FP(fit,sensitivity)
FP.obs<-temp$FP
cutpoint<-temp$cutpoint

#null (random) distribution of FP
if(test==TRUE){
	FP.null<-NULL
	if(is.null(fit$na.action)==TRUE){
		z2<-fit$data
		}
	else{
		z2<-fit$data[-fit$na.action,]
		}
	nobs<-nrow(z2)
	for(i in 1:reps){
		z2[,c(as.character(fit$formula[[2]]))]<-z2[sample(1:nobs,replace=FALSE),
			c(as.character(fit$formula[[2]]))]
		e<-new.env()
		e$fit<-fit
		for(i in names(z2)) assign(i,z2[[i]],envir=e)
		temp<-glm(fit$formula,family=fit$family,offset=fit$offset,weights=fit$prior.weights,data=e)
		temp<-as.data.frame(cbind(temp$y,temp$fitted.values))
		names(temp)<-c('observed','fitted')
		cut.perm<-quantile(temp$fitted[temp$observed==1],prob=1-sensitivity)
		FP.perm<-(nrow(temp[temp$fitted>=cut.perm & temp$observed==0,])/nrow(temp[temp$observed==0,]))
		FP.null<-c(FP.null,FP.perm)
		}

	#compute p-value
	p.value<-sum(FP.null<=FP.obs)
	p.value<-(1+p.value)/(1+reps) 
	}

if(plot==TRUE){
	xmin<-min(FP.obs,min(FP.null))
	xmax<-max(FP.obs,max(FP.null))
	hist(FP.null,col='blue',xaxs='i',yaxs='i',xlim=c(xmin,xmax),
		xlab='False positive error rate',cex.lab=1.5,
		main='Randomization test of false positive error rate',...)
	abline(v=FP.obs,col='black',lty=2,lwd=2,...)
	legend('topleft',inset=c(.1,.1),legend='Observed',lty=2,lwd=2,bty='n')
	legend('topleft',inset=c(.1,.2),legend='Null',fill='blue',bty='n')
	}

z<-as.data.frame(rbind('Cutpoint'=cutpoint,'Observed commission rate'=round(FP.obs,3),
	'P-value'=format.pval(p.value,digits=3,eps=.001)))
names(z)=''
return(z)
}
concordance <-
function(use.data,avail.data,formula,method='glm',
  paired=FALSE,use.strata=NULL,avail.strata=NULL,v=10,bins=NULL,
  obex.plot=TRUE,cc.plot=TRUE,...){

#add method for split-sample validation

#check on library dependencies
require(epiR) 
if(method=='lmer' | method=='glmer') require(lme4)

#create n.vars object
if(method=='glm') n.vars<-length(attr(terms(formula),"term.labels"))
else if(method=='lmer' | method=='glmer') n.vars<-length(attr(terms(formula),"term.labels"))-1

##################################################################
cc.bins<-function(bins){
  
  #create storage objects
  bin.number<-as.data.frame(seq(1,bins))
  names(bin.number)<-'bin.id'
  bin.mids<-seq(1/(bins*2),1-(1/(bins*2)),by=1/(bins))
  bin.ranges<-seq(0,1,by=1/bins)
  
  #tally validation points by bin
  bin.id<-cut(valid.out,breaks=bin.ranges,labels=FALSE)    
  temp<-as.data.frame(table(bin.id))
  temp$bin.id<-as.numeric(as.character(temp$bin.id))
  temp<-merge(temp,bin.number,all=TRUE)
  temp$Freq[is.na(temp$Freq)]<-0
  N.obs<-temp$Freq
  P.obs<-N.obs/sum(N.obs)
  
  #tally available points by bin
  bin.id<-cut(avail.out,breaks=bin.ranges,labels=FALSE)    
  temp<-as.data.frame(table(bin.id))
  temp$bin.id<-as.numeric(as.character(temp$bin.id))
  temp<-merge(temp,bin.number,all=TRUE)
  temp$Freq[is.na(temp$Freq)]<-0
  N.avail<-temp$Freq
  
  #compute expected utilization and counts by bin
  denom<-sum(bin.mids*N.avail)
  P.exp<-(bin.mids*N.avail)/denom
  N.exp<-P.exp*sum(N.obs)
  
  #compute coefficient of concordance
  #cc<-epi.ccc(N.exp,N.obs,ci="z-transform",conf.level=0.95)$rho.c
  cc<-epi.ccc(P.exp,P.obs,ci="z-transform",conf.level=0.95)$rho.c

  #create list object
  z<-list(v,bins,bin.mids,N.obs,P.obs,N.avail,N.exp,P.exp,cc)
  names(z)<-c('v-folds','number of bins','bin midpoint probabilities',
    'observed counts','observed proportions','available counts',
    'expected counts','expected proportions','coefficient of concordance')
  return(z) 
}

##################################################################

#for paired logistic regression
if(paired==TRUE){

  ##create indictor variable for v-fold cross-validation
  #for random v-folds
  if(is.null(use.strata)){

    if(v>1){
      
      #shuffle use data
      use.data<-use.data[sample(1:nrow(use.data),replace=FALSE),]
    
      #create subset indicator
      use.data$v<-cut(1:nrow(use.data),breaks=v,labels=FALSE)
      }
    
    else use.data$v<-1
    }

  #for designated v-strata
  else{

    #for use strata only
    if(is.null(avail.strata)){
      use.data$v<-as.numeric(as.factor(use.strata))
      v<-length(unique(use.strata))
      }

    #for both use and avail strata
    else{
      v.use<-unique(use.strata)
      v.avail<-unique(avail.strata)
      if(!all(v.use %in% v.avail)) stop('use and available must have the same strata')
      if(!all(v.avail %in% v.use)) stop('use and available must have the same strata')
      v<-length(unique(use.strata))
      use.data$v<-as.numeric(as.factor(use.strata))
      avail.data$v<-as.numeric(as.factor(avail.strata))
      }
    }
  
  #create storage objects
  coefs<-as.data.frame(matrix(NA,nrow=v,ncol=n.vars))
  se.coefs<-as.data.frame(matrix(NA,nrow=v,ncol=n.vars))
  valid.out<-NULL
  avail.out<-NULL
  
  #conduct v-fold cross-validation
  for(i in 1:v){ 
    
    #create training data subset
    if(v>1) train.use<-use.data[use.data$v!=i,]    
    else train.use<-use.data
    
    #fit model to training data subset
    if(method=='glm') fit<-glm(formula,data=train.use,family='binomial')  
    else if(method=='lmer' | method=='glmer') fit<-glmer(formula,data=train.use,family='binomial')
    
    #save coefficients and standard errors
    if(method=='glm'){
      coefs[i,]<-coef(fit)
      names(coefs)<-names(coef(fit))
      Vcov<-vcov(fit,useScale=F)
      se.coefs[i,]<-sqrt(diag(Vcov))
      names(se.coefs)<-names(coef(fit))
    }
    else if(method=='lmer' | method=='glmer'){
      coefs[i,]<-fixef(fit)
      names(coefs)<-names(fixef(fit))
      Vcov<-vcov(fit,useScale=F)
      se.coefs[i,]<-sqrt(diag(Vcov))
      names(se.coefs)<-names(fixef(fit))
      }
    
    #create validation use dataset
    if(v>1) valid.use<-use.data[use.data$v==i,]
    else valid.use<-use.data
    
    #predict response for validation use points
    mm<-model.matrix(terms(fit),valid.use)
    if(method=='glm') newdata=mm %*% coef(fit)
    else if(method=='lmer' | method=='glmer') newdata=mm %*% fixef(fit)
    valid.out<-c(valid.out,plogis(newdata))
    
    #create validation available dataset
    if(!is.null(avail.strata)) valid.avail<-avail.data[avail.data$v==i,]  
    else valid.avail<-avail.data
    
    #predict response for available points
    mm<-model.matrix(terms(fit),valid.avail)
    if(method=='glm') newdata=mm %*% coef(fit)
    else if(method=='lmer' | method=='glmer') newdata=mm %*% fixef(fit)
    avail.out<-c(avail.out,plogis(newdata))
    }
  }
  
#for unpaired logistic regression
else{

  ##create indictor variable for v-fold cross-validation
  #for random v-folds
  if(is.null(use.strata)){
    
    if(v>1){
      
      #shuffle data
      use.data<-use.data[sample(1:nrow(use.data),replace=FALSE),]
      avail.data<-avail.data[sample(1:nrow(avail.data),replace=FALSE),]
      
      #create subset indicator
      use.data$v<-cut(1:nrow(use.data),breaks=v,labels=FALSE)
      avail.data$v<-cut(1:nrow(avail.data),breaks=v,labels=FALSE)
      }
    
    else{
      use.data$v<-1
      avail.data$v<-1
      }
    }
  
  #for designated v-strata
  else{
    
    if(is.null(avail.strata)) stop('must provide avail.strata for unpaired design')
    v.use<-unique(use.strata)
    v.avail<-unique(avail.strata)
    if(!all(v.use %in% v.avail)) stop('use and available must have the same strata')
    if(!all(v.avail %in% v.use)) stop('use and available must have the same strata')
    v<-length(unique(use.strata))
    use.data$v<-as.numeric(as.factor(use.strata))
    avail.data$v<-as.numeric(as.factor(avail.strata))
    }
  
  #create storage objects
  coefs<-as.data.frame(matrix(NA,nrow=v,ncol=n.vars+1))
  se.coefs<-as.data.frame(matrix(NA,nrow=v,ncol=n.vars+1))
  valid.out<-NULL
  avail.out<-NULL

  #conduct v-fold cross-validation
  for(i in 1:v){ 
    
    #create training datasets
    if(v>1){
      train.use<-use.data[use.data$v!=i,]
      train.avail<-avail.data[avail.data$v!=i,]
      train<-rbind(train.use,train.avail)
      }
    else{
      train.use<-use.data
      train.avail<-avail.data
      train<-rbind(train.use,train.avail)
      }
    
    #fit model to training data subset
    if(method=='glm') fit<-glm(formula,data=train,family='binomial')  
    else if(method=='lmer' | method=='glmer') fit<-glmer(formula,data=train,family='binomial')
    
    #save coefficients and standard errors
    if(method=='glm'){
      coefs[i,]<-coef(fit)
      names(coefs)<-names(coef(fit))
      Vcov<-vcov(fit,useScale=F)
      se.coefs[i,]<-sqrt(diag(Vcov))
      names(se.coefs)<-names(coef(fit))
      }
    else if(method=='lmer' | method=='glmer'){
      coefs[i,]<-fixef(fit)
      names(coefs)<-names(fixef(fit))
      Vcov<-vcov(fit,useScale=F)
      se.coefs[i,]<-sqrt(diag(Vcov))
      names(se.coefs)<-names(fixef(fit))
      }

    #create validation use dataset
    if(v>1) valid.use<-use.data[use.data$v==i,]
    else valid.use<-use.data
    
    #predict response for validation use points
    mm<-model.matrix(terms(fit),valid.use)
    if(method=='glm') newdata=mm %*% coef(fit)
    else if(method=='lmer' | method=='glmer') newdata=mm %*% fixef(fit)
    valid.out<-c(valid.out,plogis(newdata))

    #create validation available dataset
    if(v>1) valid.avail<-avail.data[avail.data$v==i,]  
    else valid.avail<-avail.data
    
    #predict response for available points
    mm<-model.matrix(terms(fit),valid.avail)
    if(method=='glm') newdata=mm %*% coef(fit)
    else if(method=='lmer' | method=='glmer') newdata=mm %*% fixef(fit)
    avail.out<-c(avail.out,plogis(newdata))
    }
  }  

#brute force optimization of bin number
if(is.null(bins)){
  
  #create objects
  bin.seq<-seq(5,20)
  cc.out<-matrix(NA,nrow=length(bin.seq),ncol=3)
  cc.results<-vector('list',length=length(bin.seq))

  #loop over bins
  for(i in 1:length(bin.seq)){
    cc.results[[i]]<-cc.bins(bins=bin.seq[i])
    cc.out[i,]<-as.numeric(cc.results[[i]][[9]])
    }

  #find bin number for max cc
  cc.out<-as.data.frame(cbind(bin.seq,cc.out))
  names(cc.out)<-c('bins','est','lower','upper')
  
  #create output
  z<-cc.results[cc.out$est==max(cc.out$est)]
  z<-z[[1]]
  bins<-z[[2]]
  
  #optional cc plot
  if(cc.plot){
    matplot(cc.out$bins,cc.out[,-1],typ='l',lty=c(1,2,2),col=c(1,1,1),
      xlab='Number of Bins',ylab='Coefficient of Concordance',...)
    points(cc.out$bins[cc.out$est==max(cc.out$est)],
      cc.out$est[cc.out$est==max(cc.out$est)],pch=19,col='red')
    readline('Press return for next plot')
    }
  }

#user specified bin number
else{
  z<-cc.bins(bins=bins)  
  }

#plot observed vs expected proportions
if(obex.plot){
  P.exp<-as.numeric(z[[8]])
  P.obs<-as.numeric(z[[5]])
  plot(P.exp,P.obs,pch=19,
    xlab='Expected Proportion',ylab='Observed Proportion',...)
  abline(0,1,lty=2,lwd=2)
  text(mean(P.exp),max(P.obs),
    paste("Concord. coef. = ",round(z[[9]][1],2),' (',round(z[[9]][2],2),
      ' - ',round(z[[9]][3],2),'); ',
      v,'-fold; ',bins,' bins',sep=''))
  }

#create final list object
zz<-list(coefficients=coefs,standard.errors=se.coefs,coefficient.concordance=z)

return(zz)
}
contrast.matrix <-
function(grp){

grp<-as.factor(grp)
N<-length(grp)
grp.mat<-matrix(1,nrow=N,ncol=N)
matched<-function(irow,icol,grp){
	    grp[irow]==grp[icol]
		}
irow<-row(matrix(nrow=N,ncol=N))
icol<-col(matrix(nrow=N,ncol=N))
grp.mat[matched(irow,icol,grp)]<-0
z<-as.dist(grp.mat)
return(z)
}
cov.test <-
function(x,groups,var='',method='bartlett',...){

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

#create summary table
z<-matrix(0,ncol(y),2) #create blank matrix
for(i in 1:ncol(y)){ #loop thru variables
	if(method=='bartlett'){
		temp<-bartlett.test(y[,i],groups,...)
		}
	else if(method=='fligner'){
		temp<-fligner.test(y[,i],groups,...)
		}	
	z[i,1]<-temp$statistic
	z[i,2]<-temp$p.value
	}
z<-round(z,3)
rownames(z)<-colnames(y)

if(method=='bartlett'){
	colnames(z)<-c('Bartletts K-squared','p-value')
	cat('Bartlett Test of Homogeneity of Variances:\n')
	}
else if(method=='fligner'){
	colnames(z)<-c('Median chi-squared','p-value')
	cat('Fligner-Killeen Test of Homogeneity of Variances:\n')
	}

return(z)
}
data.dist <-
function(x,method,var='',cor.method='pearson',abs=FALSE,
	outfile='',binary=FALSE,diag=FALSE,upper=FALSE,na.rm=TRUE,...){

library(vegan) #load vegan library
library(MASS) #load MASS library

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

#ecological distance calculations
if(method=='correlation'){ #compute correlation distance
	y<-t(y) #transpose dataset
	if(abs==TRUE) z<-as.dist(1-abs(cor(y,method=cor.method,use='complete.obs',...))) #compute distance
	else z<-as.dist((1-cor(y,method=cor.method,use='complete.obs',...))/2) #compute distance
	if(!outfile==''){ #save outfile
		zz<-as.matrix(z)
		write.table(zz,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		} #end save outfile
	} #end correlation distance
else if(method=='mahalanobis'){ #compute pairwise mahalanobis distance
  z<-matrix(nrow=nrow(y),ncol=nrow(y))
  for(i in 1:nrow(y)){
    for(j in 1:nrow(y)){
      z[i,j]<-mahalanobis(as.numeric(y[i,]),center=as.numeric(y[j,]),cov=cov(y))
    }
  }
  z<-as.dist(sqrt(z))  
  if(!outfile==''){ #save outfile
    zz<-as.matrix(z)
    write.table(zz,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
    } #end save outfile
  } #end mahalanobis distance
else { #compute all other distance methods
	z<-vegdist(y,method=method,binary=binary,diag=diag,
		upper=upper,na.rm=na.rm,...) #vegdist function
	if(!outfile==''){ #save outfile
		zz<-as.matrix(z)
		write.table(zz,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		} #end save outfile
	} #end all other distance methods
return(z)
}
data.stand <-
function(x,method,var='',margin='column',
	outfile='',plot=TRUE,save.plot=FALSE,na.rm=TRUE,
	col.hist='blue',col.line='black',las=1,lab=c(5,5,4),...){

#things to do: add option for standardizing within groups

library(vegan) #load vegan library

if(plot==TRUE){
	old.par<-par(no.readonly=TRUE)
	}

x<-as.data.frame(x)

if(!var==''){
	y1<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y2<-subset(x,select=-eval(parse(text=var))) #select remaining variables
	t1<-y1 #copy to work file for transformations
	}
else{
	y1<-x #original variables
	t1<-x #copy to work file for transformations
	}

#paired histogram function for comparing raw and standardized variable
paired.hist<-function(raw=y1,trans=t1){
	for(i in 1:ncol(raw)){ #loop thru selected variables
		par(mfrow=c(2,1),mar=c(5,5,4,2)) #graphics settings
		hist(raw[[i]],prob=TRUE,
		col=col.hist,las=las,lab=lab,
		xaxs='i',yaxs='i',xlab=names(raw[i]),
		main=paste('Histogram of',names(raw[i]),sep=' '),...)
		par(new=TRUE)
		plot(density(raw[[i]],na.rm=TRUE),
		col=col.line,las=las,lab=lab,
		xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',main='',...)
		if(method=='wisconsin'){ #wisconsin plot
			hist(trans[[i]],prob=TRUE,
			col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlab=names(trans[i]),
			main=paste(method,'standardization',sep=' '),...)
			par(new=TRUE)
			plot(density(trans[[i]],na.rm=TRUE),
			col=col.line,las=las,lab=lab,
			xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',main='',...)
			} #end wisconsin plot
		else{ #all other standardizations
			hist(trans[[i]],prob=TRUE,
			col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlab=names(trans[i]),
			main=paste(margin,method,'standardization',sep=' '),...)
			par(new=TRUE)
			plot(density(trans[[i]],na.rm=TRUE),
			col=col.line,las=las,lab=lab,
			xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',main='',...)
			} #end all other plots
		if(save.plot==TRUE){ #save plot to file
			dev.print(jpeg,file=paste('shist.',names(raw[i]),'.jpg',sep=''),width=800,height=600)
			} #end save	plot
		readline("Press return for next plot ")
		} #end loop thru variables
	} #end paired histogram function

#standardizations
if (method=='wisconsin'){ #wisconsin standardization
	t1<-decostand(t1,method='max',na.rm=TRUE,...) #column max standardization
	t1<-decostand(t1,method='total',na.rm=TRUE,...) #row total standardization
	if(plot==TRUE){ #plot paired histograms
		paired.hist(y1,t1)
		} #end plot paired histograms
	if(!var==''){
		z<-cbind(y2,t1) #bind result to remaining variables
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end wisconsin standardization
else if (method=='chi.square'){ #chi-square standardization
	t1<-decostand(t1,method=method,MARGIN=1,na.rm=na.rm,...) #decostand function
	if(plot==TRUE){ #plot paired histograms
		paired.hist(y1,t1)
		} #end plot paired histograms
	if(!var==''){
		z<-cbind(y2,t1) #bind result to remaining variables
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end chi.square standardization
else { #all other standardizations
	if (margin=='row') {MARGIN=1} #assign MARGIN parameter for decostand(vegan)
	else if (margin=='column') {MARGIN=2} #assign MARGIN parameter for decostand(vegan)
	t1<-decostand(t1,method=method,MARGIN=MARGIN,na.rm=na.rm,...) #decostand function
	if(plot==TRUE){ #plot paired histograms
		paired.hist(y1,t1)
		} #end plot paired histograms
	if(!var==''){
		z<-cbind(y2,t1) #bind result to remaining variables
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end all other standardizations

if(plot==TRUE){
	par(old.par)
	}
return(z)
}
data.trans <-
function(x,method,var='',exp=1,outfile='',
	plot=TRUE,save.plot=FALSE,col.hist='blue',col.line='black',
	las=1,lab=c(5,5,4),...){

if(plot==TRUE){
	old.par<-par(no.readonly=TRUE)
	}

if(!var==''){
	y1<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y2<-subset(x,select=-eval(parse(text=var))) #select remaining variables
	t1<-y1 #copy to work file for transformations
	}
else{
	y1<-x #original variables
	t1<-x #copy to work file for transformations
	}

#paired histogram function for comparing raw and transformed variable
paired.hist<-function(raw,trans){
	for(i in 1:ncol(raw)){ #loop thru selected variables
		par(mfrow=c(2,1),mar=c(5,5,4,2)) #graphics settings
		hist(raw[[i]],prob=TRUE,col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlab=names(raw[i]),
			main=paste('Histogram of',names(raw[i]),sep=' '),...)
		par(new=TRUE)
		plot(density(raw[[i]],na.rm=TRUE),
			xaxs='i',yaxs='i',col=col.line,las=las,lab=lab,
			axes=FALSE,xlab='',ylab='',main='',...)
		if(method=='power'){ #plot title for power transform
			hist(trans[[i]],prob=TRUE,col=col.hist,las=las,lab=lab,
				xaxs='i',yaxs='i',xlab=names(trans[i]),
				main=paste(method,'(x^',exp,')','transformation',sep=' '),...)
			} #end plot title for power transform
		else{ #plot title for other transforms
			hist(trans[[i]],prob=TRUE,col=col.hist,las=las,lab=lab,
				xaxs='i',yaxs='i',xlab=names(trans[i]),
				main=paste(method,'transformation',sep=' '),...)
			} #end plot title for other transforms
		par(new=TRUE)
		plot(density(trans[[i]],na.rm=TRUE),
			xaxs='i',yaxs='i',col=col.line,las=las,lab=lab,
			axes=FALSE,xlab='',ylab='',main='',...)
		if(save.plot==TRUE){
			dev.print(jpeg,file=paste('thist.',method,'.',names(raw[i]),'.jpg',sep=''),width=800,height=600)
			} #end save	
		readline("Press return for next plot ")
		} #end loop thru variables
	} #end paired histogram function

#transformations
if(method=='power'){ #power transformations
	if(exp==0){ #binary transformation
		t1[t1!=0]<-1 #compute binary
		t1<-as.data.frame(round(t1,3)) #round to 3 decimal places
		if(plot==TRUE){
			paired.hist(y1,t1)
			}
		if(!var==''){
			z<-cbind(y2,t1) #bind transformed variables with remaining
			}
		else{z<-t1}
		if(!outfile==''){ #save outfile
			write.csv(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE)
			} #end save outfile
		} #end binary transformation
	else{ #other power transformations
		t1<-t1^exp #compute power transformation
		t1<-as.data.frame(round(t1,3)) #round to 3 decimal places
		if(plot==TRUE){
			paired.hist(y1,t1)
			}
		if(!var==''){
			z<-cbind(y2,t1) #bind transformed variables with remaining
			}
		else{z<-t1}
		if(!outfile==''){ #save outfile
			write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
			} #end save outfile
		} #end other power transformations
	} #end power transformation

else if(method=='log'){ #log10 transformation
	q<-min(t1[t1!=0],na.rm=TRUE) #find minimum non-zero value
	c<-trunc(log10(q)) #order of magnitude constant
	d<-10^(c) #decimal constant (inverse log10 of c
	t1<-log10(t1+d)-c #compute log10 transformation
	t1<-as.data.frame(round(t1,3)) #round to 3 decimal places
	if(plot==TRUE){
		paired.hist(y1,t1)
		}
	if(!var==''){
		z<-cbind(y2,t1) #bind transformed variables with remaining
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end log transformation
	
else{ #arcsine squareroot transformation
	t1<-(2/pi)*asin(sqrt(t1)) #compute arcsin square root
	t1<-as.data.frame(round(t1,3)) #round to 3 decimal places
	if(plot==TRUE){
		paired.hist(y1,t1)
		}
	if(!var==''){
		z<-cbind(y2,t1) #bind transformed variables with remaining
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end arcsine squareroot transformation
	
if(plot==TRUE){
	par(old.par)
	}

return(z)
}
dbeta.plot <-
function(shape1=1,shape2=1,ylim=c(0,5),...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
plot(0,0,type='n',xlim=c(0,1),ylim=ylim,xlab='z',ylab='probability density',
  main='Beta Probability Density',cex.lab=1.5,cex.main=1.5,...)
legend<-NULL
temp<-seq(1,length(shape1))
m<-rep(temp,each=length(shape2))
k<-1
for(i in shape1){
	for(j in shape2){
		curve(dbeta(x,shape1=i,shape2=j),add=TRUE,col=m[k],lty=k,...)
		k<-k+1
		t<-paste('a=',i,', b=',j,sep='')
		legend<-rbind(legend,t)
		}
	}
legend(x='top',inset=c(.01,.01),legend=legend,col=m,lty=seq(1,k),cex=1.5,...)
}
dbinom.plot <-
function(size=10,prob=.5,col=NULL,...){
legend<-NULL
out<-NULL
n<-seq(0,size)
for(j in prob){
	t1<-dbinom(n,size=size,prob=j)
	out<-cbind(out,t1)
	t2<-paste('size=',size,', prob=',j,sep='')
	legend<-rbind(legend,t2)
	}
if(is.null(col)) col.bar=seq(1:length(prob)) 
else col.bar=col
barplot(t(out),beside=TRUE,names=as.character(n),col=col.bar,
	xlab='# successes',ylab='probability',main='Binomal Probability Mass',
  cex.main=1.5,cex.lab=1.5,...)
abline(h=0,col='gray')
legend(x='topright',inset=c(.01,.01),legend=legend,fill=col.bar,cex=1,bty='n')
}
dchisq.plot <-
function(df=1,ncp=0,xlim=c(0,10),ylim=c(0,1),...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab='z',ylab='probability density',
  main='Chi-square Probability Density',cex.lab=1.5,cex.main=1.5,...)
legend<-NULL
k<-1
for(i in df){
	curve(dchisq(x,df=i,ncp=ncp),xlim[1],xlim[2],add=TRUE,lty=k,...)
	k<-k+1
	t<-paste('df=',i,sep='')
	legend<-rbind(legend,t)
	}
legend(x='topright',inset=c(.01,.01),legend=legend,lty=seq(1,k),cex=1.5,...)
}
dexp.plot <-
function(rate=1,xlim=c(0,15),ylim=c(0,1),xlab='z',...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab='probability density',
  main='Exponential Probability Density',cex.lab=1.5,cex.main=1.5,...)
legend<-NULL
k<-1
for(i in rate){
	curve(dexp(x,rate=i),xlim[1],xlim[2],add=TRUE,lty=k,...)
	k<-k+1
	t<-paste('rate=',i,sep='')
	legend<-rbind(legend,t)
	}
legend(x='topright',inset=c(.01,.01),legend=legend,lty=seq(1,k),cex=1.5,...)
}
df.plot <-
function(df1=1,df2=10,xlim=c(0,6),ylim=c(0,1),...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab='z',ylab='probability density',
  main='Fishers F Probability Density',cex.lab=1.5,cex.main=1.5,...)
legend<-NULL
temp<-seq(1,length(df1))
m<-rep(temp,each=length(df2))
k<-1
for(i in df1){
	for(j in df2){
		curve(df(x,df1=i,df2=j),xlim[1],xlim[2],add=TRUE,col=m[k],lty=k,...)
		k<-k+1
		t<-paste('df1=',i,', df2=',j,sep='')
		legend<-rbind(legend,t)
		}
	}
legend(x='topright',inset=c(.01,.01),legend=legend,col=m,lty=seq(1,k),cex=1.5,...)
}
dgamma.plot <-
function(shape=1,scale=1,xlim=c(0,25),ylim=c(0,1),...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab='z',ylab='probability density',
  main='Gamma Probability Density',cex.lab=1.5,cex.main=1.5,...)
legend<-NULL
temp<-seq(1,length(shape))
m<-rep(temp,each=length(scale))
k<-1
for(i in shape){
	for(j in scale){
		curve(dgamma(x,shape=i,scale=j),xlim[1],xlim[2],add=TRUE,col=m[k],lty=k,...)
		k<-k+1
		t<-paste('shape=',i,', scale=',j,sep='')
		legend<-rbind(legend,t)
		}
	}
legend(x='topright',inset=c(.01,.01),legend=legend,col=m,lty=seq(1,k),cex=1.5,...)
}
dgeom.plot <-
function(events=25,prob=.5,...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
legend<-NULL
out<-NULL
n<-seq(0,events)
for(i in prob){
		t1<-dgeom(n,prob=i)
		out<-cbind(out,t1)
		t2<-paste('prob=',i,sep='')
		legend<-rbind(legend,t2)
	}
m<-seq(1:length(prob))
barplot(t(out),beside=TRUE,names=as.character(n),col=m,
	xlab='count',ylab='probability',main='Geometric Probability Mass',
  cex.lab=1.5,cex.main=1.5,...)
legend(x='top',inset=c(.01,.01),legend=legend,fill=m,cex=1.5)
}
dist.plots <-
function(x,groups,distance='euclidean',na.rm=TRUE,
	col='blue',col.line='red',las=1,...){

#create vector of dissimilarities for plotting
if(inherits(x,'dist')) 
    dmat<-x
else if(is.matrix(x)&&nrow(x)==ncol(x)&&all(x[lower.tri(x)]==t(x)[lower.tri(x)])){
    dmat<-x
    attr(dmat,'method')<-'user supplied square matrix'
    }
else dmat<-vegdist(x,method=distance)
distance<-attr(dmat,'method')
dmat<-as.matrix(dmat)
diag(dmat)<-NA
x<-as.dist(dmat)
groups<-as.factor(groups)
matched<-function(irow,icol,groups){
    groups[irow]==groups[icol]
	}
diss.vec<-as.vector(x)
N<-attributes(x)$Size
irow<-as.vector(as.dist(row(matrix(nrow=N,ncol=N))))
icol<-as.vector(as.dist(col(matrix(nrow=N,ncol=N))))
within<-matched(irow,icol,groups)
class.vec<-rep("Between",length(diss.vec))
take<-as.numeric(irow[within])
class.vec[within]<-levels(groups)[groups[take]]
class.vec<-factor(class.vec,levels=c("Between",levels(groups)))

#boxplot of between and within-group dissimilarities
boxplot(diss.vec~class.vec,notch=TRUE,varwidth=TRUE,col=col,
	ylab='Dissimilarity',main='Between- and Within-group Dissimilarities',...)
readline("Press return for next plot ")

#histograms of w/i group dissimilarities
grp<-levels(groups) #create vector of group names
s<-floor(length(grp)^.5) #create multi-figure dimensions
s<-c(s,ceiling(length(grp)/s)) #create multi-figure dimensions
old.par<-par(no.readonly=TRUE)
par(mfrow=s,mar=c(5,5,4,2)) #graphics settings
for(j in 1:length(grp)){ #loop thru groups
	z<-diss.vec[class.vec==grp[j]] #select records for group j
	hist(z,prob=TRUE,col=col,las=las,xlim=range(diss.vec),
		xaxs='i',yaxs='i',xlab=grp[j],main='',...)
	par(new=TRUE)
	plot(density(z,na.rm=na.rm),xaxs='i',yaxs='i',
		col=col.line,las=las,xlim=range(diss.vec),
		axes=FALSE,xlab='',ylab='',main='',...)
	} #end loop thru groups
par(old.par)
readline("Press return for next plot ")

#qqnorm plots of w/i group dissimilarities
s<-floor(length(grp)^.5) #create multi-figure dimensions
s<-c(s,ceiling(length(grp)/s)) #create multi-figure dimensions
old.par<-par(no.readonly=TRUE)
par(mfrow=s,mar=c(5,5,4,2)) #graphics settings
for(j in 1:length(grp)){ #loop thru groups
	z<-diss.vec[class.vec==grp[j]] #select records for group j
	qqnorm(z,datax=TRUE,col=col,las=las,main=grp[j],
		ylab='Sample quantiles',...)
	z.IQR<-IQR(z,na.rm=na.rm)
	if(z.IQR>0){ #add qqline
		qqline(z,datax=TRUE,col=col.line,...)
		} #end qqline
	} #end loop thru groups
par(old.par)

}
dlnorm.plot <-
function(meanlog=0,sdlog=1,xlim=c(0,15),ylim=c(0,2),...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab='z',ylab='probability density',
  main='Lognormal Probability Density',cex.lab=1.5,cex.main=1.5,...)
legend<-NULL
temp<-seq(1,length(meanlog))
m<-rep(temp,each=length(sdlog))
k<-1
for(i in meanlog){
	for(j in sdlog){
		curve(dlnorm(x,meanlog=i,sdlog=j),xlim[1],xlim[2],add=TRUE,col=m[k],lty=k,...)
		k<-k+1
		t<-paste('meanlog=',i,', sdlog=',j,sep='')
		legend<-rbind(legend,t)
		}
	}
legend(x='topright',inset=c(.01,.01),legend=legend,col=m,lty=seq(1,k),cex=1.5,...)
}
dnbinom.plot <-
function(events=25,mu=2,size=1,...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
legend<-NULL
out<-NULL
n<-seq(0,events)
for(i in mu){
	for(j in size){
  		t1<-dnbinom(n,mu=i,size=j)
  		out<-cbind(out,t1)
  		t2<-paste('mu=',i,', k=',j,sep='')
  		legend<-rbind(legend,t2)
  		}
	}
m<-seq(1:(length(mu)*length(size)))
barplot(t(out),beside=TRUE,names=as.character(n),col=m,
	xlab='count',ylab='probability',main='Negative Binomal Probability Mass',
  cex.lab=1.5,cex.main=1.5,...)
legend(x='top',inset=c(.01,.01),legend=legend,fill=m,cex=1.5)
}
dnorm.plot <-
function(mean=0,sd=1,xlim=c(-3,3),ylim=c(0,1),...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab='z',ylab='probability density',
  main='Normal Probability Density',cex.lab=1.5,cex.main=1.5,...)
legend<-NULL
temp<-seq(1,length(mean))
m<-rep(temp,each=length(sd))
k<-1
for(i in mean){
	for(j in sd){
		curve(dnorm(x,mean=i,sd=j),xlim[1],xlim[2],add=TRUE,col=m[k],lty=k,...)
		k<-k+1
		t<-paste('mean=',i,', sd=',j,sep='')
		legend<-rbind(legend,t)
		}
	}
legend(x='topright',inset=c(.01,.01),legend=legend,col=m,lty=seq(1,k),cex=1.5,...)
}
dpois.plot <-
function(events=25,lambda=1,...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
legend<-NULL
out<-NULL
n<-seq(0,events)
for(j in lambda){
	t1<-dpois(n,lambda=j)
	out<-cbind(out,t1)
	t2<-paste('lambda=',j,sep='')
	legend<-rbind(legend,t2)
	}
barplot(t(out),beside=TRUE,names=as.character(n),col=seq(1:length(lambda)),
	xlab='# events',ylab='probability',main='Poisson Probability Mass',
  cex.lab=1.5,cex.main=1.5,...)
legend(x='top',inset=c(.01,.01),legend=legend,fill=seq(1:length(lambda)),cex=1.5)
}
drop.var <-
function(x,var='',outfile='',min.cv=0,min.po=0,min.fo=0,
	max.po=100,max.fo=nrow(x),pct.missing=100){

if(!var==''){
	y1<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y2<-subset(x,select=-eval(parse(text=var))) #select remaining variables
	}
else{
	y1<-x #original variables
	}

#statistical functions
cv<<-function(x,na.rm) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100 
freq.occur<<-function(x,na.rm) sum(!x==0,na.rm=TRUE)
pct.occur<<-function(x,na.rm) sum(!x==0,na.rm=TRUE)/length(x)*100
fpct.missing<<-function(x,na.rm) sum(is.na(x))/length(x)*100 

#delete offending variables
z<-as.matrix(apply(y1,2,cv)<min.cv | 
apply(y1,2,pct.occur)>max.po |
apply(y1,2,pct.occur)<min.po |
apply(y1,2,freq.occur)>max.fo |
apply(y1,2,freq.occur)<min.fo | 
apply(y1,2,fpct.missing)>pct.missing)
z<-y1[,z[,1]==FALSE]
if(!var==''){
	z<-cbind(y2,z) #bind reduced set w/ deselected variables
	}
if(!outfile==''){ #save outfile
	write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
	} #end save outfile
return(z)
}
dt.plot <-
function(df=1,ncp=0,xlim=c(-3,3),ylim=c(0,1),...){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab='z',ylab='probability density',
  main='Students t Probability Density',cex.lab=1.5,cex.main=1.5,...)
legend<-NULL
k<-1
for(i in df){
	curve(dt(x,df=i,ncp=ncp),xlim[1],xlim[2],add=TRUE,lty=k,...)
	k<-k+1
	t<-paste('df=',i,sep='')
	legend<-rbind(legend,t)
	}
legend(x='topright',inset=c(.01,.01),legend=legend,lty=seq(1,k),cex=1.5,...)
}
edf.plots <-
function(x,var='',by='',...){

old.par<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

if(by==''){ #ecdf w/o groups
	for(i in 1:ncol(y)){ #loop thru variables
		plot(sort(y[[i]]),type='o',col='blue',
		yaxs='i',xaxs='i',las=1,lab=c(10,10,3),
		xlab='Plot Rank',ylab=names(y[i]),
		main='Empirical Distribution Function')
		line(c(1,length(y[[i]])),range(y[[i]]))
		if(i != ncol(y)) readline("Press return for next plot")
		} #end loop thru variables
	} #end edf w/o groups
	
else{ #ecdf plots w/ groups
	n<-by.names(x,by) #create by variable
	y<-cbind(n,y) #bind with selected variables
	groups<-levels(y[,2]) #create vector of group names
	s<-floor(length(groups)^.5) #create multi-figure dimension
	s<-c(s,ceiling(length(groups)/s)) #create multi-figure dimensions
	for(i in 3:ncol(y)){ #loop thru variables
		par(mfrow=s,mar=c(5,5,4,2))	#graphics settings	
		for(j in 1:length(groups)){ #loop thru groups
			z<-y[y[,2]==groups[j],] #select records for group j
			plot(sort(z[[i]]),type='o',col='blue',
			yaxs='i',xaxs='i',las=1,lab=c(10,10,3),
			xlab='Plot Rank',ylab=names(z[i]),
			main=groups[j])
			line(c(1,length(z[[i]])),range(z[[i]]))
			} #end loop thru groups
			if(i != ncol(y) & j != length(groups)) readline("Press return for next plot")
		} #end loop thru variables
	} #end edf w/ groups
par(old.par)
}
ecdf.plots <-
function(x,var='',by='',...){

old.par<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

if(by==''){ #ecdf w/o groups
	for(i in 1:ncol(y)){ #loop thru variables
		plot(ecdf(y[[i]]),verticals=TRUE,col.01line = "gray70",
		xlab=names(y[i]),ylab='Cumulative probability',
		main='Empirical Cumulative Distribution Function')
		if(i != ncol(y)) readline("Press return for next plot")
		} #end loop thru variables
	} #end ecdf w/o groups
	
else{ #ecdf plots w/ groups
	n<-by.names(x,by) #create by variable
	y<-cbind(n,y) #bind with selected variables
	groups<-levels(y[,2]) #create vector of group names
	s<-floor(length(groups)^.5) #create multi-figure dimension
	s<-c(s,ceiling(length(groups)/s)) #create multi-figure dimensions
	for(i in 3:ncol(y)){ #loop thru variables
		par(mfrow=s,mar=c(5,5,4,2))	#graphics settings	
		for(j in 1:length(groups)){ #loop thru groups
			z<-y[y[,2]==groups[j],] #select records for group j
			plot(ecdf(z[[i]]),verticals=TRUE,col.01line = "gray70",
			xlab=names(z[i]),ylab='Cumulative probability',
			main=groups[j])
			} #end loop thru groups
			if(i != ncol(y)) {readline("Press return for next plot")}
		} #end loop thru variables
	} #end ecdf w/ groups
par(old.par)
}
exp.plot <-
function(rate=1,xlim=c(0,15)){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(rate)
k<-0
for(i in rate){
  k<-k+1
	curve(pexp(x,rate=i),xlim[1],xlim[2],
		xlab='Z',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
		mtext(paste('(rate=',i,')',sep=''),side=3,col='red')
	curve(dexp(x,rate=i),xlim[1],xlim[2],
		xlab='Z',ylab='Probability density',main='Probability Density')
		mtext(paste('(rate=',i,')',sep=''),side=3,col='red')
	curve(qexp(x,rate=i),0,1,
		xlab='Cumulative probability (quantile)',ylab='Z',main='Quantiles')
		mtext(paste('(rate=',i,')',sep=''),side=3,col='red')
	y<-rexp(1000,rate=i)
	hist(y,xlab='Z',ylab='Frequency',main='Random Observations')
		mtext(paste('(rate=',i,')',sep=''),side=3,col='red')
	if(k<n.plots) readline('press return for next plot')
	}
}
f.plot <-
function(df1=1,df2=10,xlim=c(0,6)){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(df1)*length(df2)
k<-0
for(i in df1){
	for(j in df2){
	  k<-k+1
		curve(pf(x,df1=i,df2=j),xlim[1],xlim[2],
			xlab='Z',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
			mtext(paste('(df1=',i,', df2=',j,')',sep=''),side=3,col='red')
		curve(df(x,df1=i,df2=j),xlim[1],xlim[2],
			xlab='Z',ylab='Probability density',main='Probability Density')
			mtext(paste('(df1=',i,', df2=',j,')',sep=''),side=3,col='red')
		curve(qf(x,df1=i,df2=j),0,1,
			xlab='Cumulative probability (quantile)',ylab='Z',main='Quantiles')
			mtext(paste('(df1=',i,', df2=',j,')',sep=''),side=3,col='red')
		y<-rf(1000,df1=i,df2=j)
		hist(y,xlab='Z',ylab='Frequency',main='Random Observations')
			mtext(paste('(df1=',i,', df2=',j,')',sep=''),side=3,col='red')
	  if(k<n.plots) readline('press return for next plot')
		}
	}
}
foa.plots <-
function(x,var='',margin='column',outfile='',na.rm=TRUE,col.hist='blue',
	col.point='blue',col.line='red',col.den='black',las=1,lab=c(10,10,3),...){

old.par<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

z1<-apply(y>0,2,sum,na.rm=na.rm) #compute species frequency of occurrence
if(min(z1)==0) stop("all species must have non-zero sum")
z2<-100*z1/nrow(y) #compute species percent occurrence
z3<-log(z1) #compute log of species frequency of occurrences
z4<-apply(y,2,sum,na.rm=na.rm)/z1 #compute species mean abundance where present
r1<-apply(y>0,1,sum,na.rm=na.rm) #compute number of species per plot
r2<-apply(y,1,sum,na.rm=na.rm) #compute total abundance per plot

#summary table
if(margin=='column'){
	z<-as.data.frame(cbind(z1,z2,z3,z4))
	z<-round(z,3)
	names(z)<-c('spc.pres','spc.perc','spc.log','spc.mean')
	}
else{z<-as.data.frame(cbind(r1,r2))
	z<-round(z,3)
	names(z)<-c('plt.pres','plt.sum')
	}
if(!outfile==''){ #save outfile
	write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
	} #end save outfile
	
#plot1: empirical distribution of species occurrence - raw scale
plot(sort(z1),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,max(z1)),xlim=c(0,length(z1)),yaxs='i',xaxs='i',
xlab='Species Rank',ylab='Frequency of Occurrence',
main='Empirical Distribution of Species Occurrence',...)
abline(h=quantile(z1,.5),lty=1,col=col.line,...)
abline(h=quantile(z1,.05),lty=2,col=col.line,...)
abline(h=quantile(z1,.95),lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot2: empirical distribution of species occurrence - relativized
plot(sort(z2),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,100),xlim=c(0,length(z2)),yaxs='i',xaxs='i',
xlab='Species Rank',ylab='Percent Occurrence',
main='Empirical Distribution of Species Relative Occurrence',...)
abline(h=50,lty=1,col=col.line,...)
abline(h=5,lty=2,col=col.line,...)
abline(h=95,lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot3: histogram of species occurrence - raw scale
hist(z1,prob=TRUE,col=col.hist,las=las,lab=lab,
xaxs='i',yaxs='i',xlab='Species Occurrence',ylab='Frequency',
main='Histogram of Species Occurrence',...)
par(new=TRUE)
plot(density(z1,na.rm=na.rm),xaxs='i',yaxs='i',
col=col.den,las=las,lab=lab,
axes=FALSE,xlab='',ylab='',main='',...)
readline("Press return for next plot ")

#plot4: histogram of species occurrence - log transformation
hist(z3,prob=TRUE,col=col.hist,las=las,lab=lab,
xaxs='i',yaxs='i',xlab='Log(Species Occurrence)',ylab='Frequency',
main='Histogram of Log(Species Occurrence)',...)
par(new=TRUE)
plot(density(z3,na.rm=na.rm),xaxs='i',yaxs='i',
col=col.den,las=las,lab=lab,
axes=FALSE,xlab='',ylab='',main='',...)
readline("Press return for next plot ")

#plot5: empirical distribution of species mean abundance - raw scale
plot(sort(z4),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,max(z4)),xlim=c(0,length(z4)),yaxs='i',xaxs='i',
xlab='Species Rank',ylab='Mean Abundance',
main='Empirical Distribution of Species Mean Abundance',...)
abline(h=quantile(z4,.5),lty=1,col=col.line,...)
abline(h=quantile(z4,.05),lty=2,col=col.line,...)
abline(h=quantile(z4,.95),lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot6: plot of frequency of occurrence versus mean abundance
plot(z1,z4,type='p',col=col.point,lab=lab,las=las,
yaxs='i',xaxs='i',
xlab='Frequency of Occurrence',ylab='Mean Abundance',
main='Species Occurrence vs Mean Abundance',...)
yorn<-readline("Do you want to identify individual species? Y/N :")
    if(yorn=='Y' || yorn=='y') 
        identify(z1,z4,names(y))
readline("Press return for next plot ")

#plot7: plot of frequency of occurrence versus log of mean abundance
plot(z1,z4,type='p',col=col.point,lab=lab,las=las,
yaxs='i',xaxs='i',log='y',
xlab='Frequency of Occurrence',ylab='Log(Mean Abundance)',
main='Species Occurrence vs Log(Mean Abundance)',...)
yorn<-readline("Do you want to identify individual species? Y/N :")
    if(yorn=='Y' || yorn=='y') 
        identify(z1,z4,names(y))
readline("Press return for next plot ")

#plot8: empirical distribution of species per plot
plot(sort(r1),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,max(r1)),xlim=c(0,length(r1)),yaxs='i',xaxs='i',
xlab='Plot Rank',ylab='Plot Richness',
main='Empirical Distribution of Plot Richness',...)
abline(h=quantile(r1,.5),lty=1,col=col.line,...)
abline(h=quantile(r1,.05),lty=2,col=col.line,...)
abline(h=quantile(r1,.95),lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot9: empirical distribution of total plot abundance
plot(sort(r2),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,max(r2)),xlim=c(0,length(r2)),yaxs='i',xaxs='i',
xlab='Plot Rank',ylab='Total Abundance',
main='Empirical Distribution of Plot Total Abundance',...)
abline(h=quantile(r2,.5),lty=1,col=col.line,...)
abline(h=quantile(r2,.05),lty=2,col=col.line,...)
abline(h=quantile(r2,.95),lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot10: plot of plot richness versus plot total abundance
plot(r1,r2,type='p',col=col.point,lab=lab,las=las,
yaxs='i',xaxs='i',
xlab='Plot Richness',ylab='Plot Total Abundance',
main='Plot Richness vs Total Abundance',...)
yorn<-readline("Do you want to identify individual plots? Y/N :")
    if(yorn=='Y' || yorn=='y') 
        identify(r1,r2,rownames(y))

par(old.par)
return(z)
}
gamma.plot <-
function(shape=1,scale=1,xlim=c(0,25)){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(shape)*length(scale)
k<-0
for(i in shape){
	for(j in scale){
	  k<-k+1
		curve(pgamma(x,shape=i,scale=j),xlim[1],xlim[2],
			xlab='Z',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
			mtext(paste('(shape=',i,', scale=',j,')',sep=''),side=3,col='red')
		curve(dgamma(x,shape=i,scale=j),xlim[1],xlim[2],
			xlab='Z',ylab='Probability density',main='Probability Density')
			mtext(paste('(shape=',i,', scale=',j,')',sep=''),side=3,col='red')
		curve(qgamma(x,shape=i,scale=j),0,1,
			xlab='Cumulative probability (quantile)',ylab='Z',main='Quantiles')
			mtext(paste('(shape=',i,', scale=',j,')',sep=''),side=3,col='red')
		y<-rgamma(1000,shape=i,scale=j)
		hist(y,xlab='Z',ylab='Frequency',main='Random Observations')
			mtext(paste('(shape=',i,', scale=',j,')',sep=''),side=3,col='red')
		if(k<n.plots) readline('press return for next plot')
		}
	}
}
geom.plot <-
function(events=25,prob=c(.2,.5,.7)){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(events)*length(prob)
k<-0
n<-seq(0,events)
for(i in prob){
  k<-k+1
	barplot(pgeom(n,prob=i),names=as.character(n),
		xlab='Count',ylab='Cumulative probability (quantile)',main='Cumulative probability')
		mtext(paste('(prob=',i,')',sep=''),side=3,col='red')
	barplot(dgeom(n,prob=i),names=as.character(n),
		xlab='Count',ylab='Probability',main='Probability Mass')
		mtext(paste('(prob=',i,')',sep=''),side=3,col='red')
	barplot(qgeom(seq(0,.95,.05),prob=i),names=seq(0,.95,.05),
		xlab='Cumulative probability (quantile)',ylab='Count',main='Quantiles')
		mtext(paste('(prob=',i,')',sep=''),side=3,col='red')
	y<-rgeom(1000,prob=i)
  y.table<-as.data.frame(table(y))
  y.levels<-as.data.frame(n)
  names(y.levels)<-'y'
  y.table<-merge(y.levels,y.table,all.x=TRUE)
  y.table$Freq[is.na(y.table$Freq)]<-0
  barplot(y.table$Freq,names=as.character(n),xlab='Count',ylab='Frequency',
    main='Random Observations',col='gray')
		mtext(paste('(prob=',i,')',sep=''),side=3,col='red')
  if(k<n.plots) readline('press return for next plot')
	}
}
hclus.cophenetic <-
function(d,hclus,fit='lm',...){

old.par<-par(no.readonly=TRUE)
d.coph<-cophenetic(hclus)
r<-round(cor(d,d.coph),2)
plot(d,d.coph,xlab='Observed dissimilarity',
	ylab='Cophenetic dissimilarity',
	main=paste('Cophenetic Correlation ',
	'(',hclus$dist.method,', ',hclus$method,')',sep=''),...)
	text(max(d),min(d.coph),paste('Cophenetic correlation = ',r,sep=''),col='red',pos=2)
#	title(sub=paste('Cophenetic correlation = ',r,sep=''),col.sub='red',adj=0)
if(fit=='lm'){
	abline(lm(d.coph~d),col='blue')
	}
else if(fit=='rlm'){
	abline(rlm(d.coph~d),col='blue')
	}
else if(fit=='qls'){
	abline(lqs(d.coph~d),col='blue')
	}
	
par(old.par)
return(r)
}
hclus.scree <-
function(x,...){

old.par<-par(no.readonly=TRUE)
z1<-seq(length(x$height),1)
z<-as.data.frame(cbind(z1,sort(x$height)))
plot(z[,1],z[,2],type='o',lwd=1.5,pch=19,col='blue',
	ylab='Dissimilarity',xlab='Number of Clusters',
	main=paste('Scree Plot of Hierarchical Clustering ',
	'(',x$dist.method,', ',x$method,')',sep=''),...)
par(old.par)
}
hclus.table <-
function(x){

z1<-x$dist.method
z2<-x$method
z3<-seq(length(x$labels)-1,1)
z3<-as.data.frame(cbind(z3,x$merge,x$height))
colnames(z3)<-c('no. clusters','entity','entity','distance')
z<-list(z1,z2,z3)
names(z)<-c('dist.method','method','cluster.table')
return(z)
}
hist.plots <-
function(x,var='',by='',save.plot=FALSE,na.rm=TRUE,
	col.hist='blue',las=1,lab=c(5,5,4),...){

old.par<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

if(by==''){ #histogram w/o groups
	par(mar=c(5,5,4,2)) #graphics settings
	for(i in 1:ncol(y)){ #loop thru variables
		par(new=FALSE)
		hist(y[[i]],col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlab=names(y[i]),
			main=paste('Histogram of',names(y[i]),sep=' '),...)
#		par(new=TRUE)
#		lines(density(y[[i]],na.rm=na.rm),xaxs='i',yaxs='i',
#			col=col.line,las=las,lab=lab,
#			axes=FALSE,xlab='',ylab='',main='',...)
		if(save.plot==TRUE){
			dev.print(jpeg,file=paste('hist.',names(y[i]),'.jpg',sep=''),width=800,height=600)
			} #end save
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end histogram w/o groups
		
else{ #histograms w/ groups
	n<-by.names(x,by) #create by variable
	y<-cbind(n,y) #bind with selected variables
	groups<-levels(y[,2]) #create vector of group names
	s<-floor(length(groups)^.5) #create multi-figure dimensions
	s<-c(s,ceiling(length(groups)/s)) #create multi-figure dimensions
	for(i in 3:ncol(y)){ #loop thru variables
		par(mfrow=s,mar=c(5,5,4,2)) #graphics settings
		for(j in 1:length(groups)){ #loop thru groups
			z<-y[y[,2]==groups[j],] #select records for group j
			y.min<-min(y[i],na.rm=na.rm)
			y.max<-max(y[i],na.rm=na.rm)
			hist(z[[i]],col=col.hist,las=las,lab=lab,
				xaxs='i',yaxs='i',xlab=names(y[i]),main=groups[j],
				xlim=c(y.min,y.max),...)
#			par(new=TRUE)
#			lines(density(z[[i]],na.rm=na.rm),xaxs='i',yaxs='i',
#				col=col.line,las=las,lab=lab,xlim=c(y.min,y.max),
#				axes=FALSE,xlab='',ylab='',main='',...)
			if(j==length(groups) & save.plot==TRUE){
				dev.print(jpeg,file=paste('hist.',names(y[i]),'.jpg',sep=''),width=800,height=600)
				} #end save
			} #end loop thru groups
			if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end histogram w/groups
par(old.par)
}
intrasetcor <-
function(object){ 

    w <- weights(object)
    lc <- sweep(object$CCA$u, 1, sqrt(w), "*")
    cor(qr.X(object$CCA$QR), lc)
}
kappa.mc <-
function(fit,test=FALSE,reps=1000,plot=TRUE,...){

owarn <- options("warn")
on.exit(options(warn=owarn$warn))
options(warn=-1)

########################################################################
#functions
########################################################################

#Kappa function

cohen.kappa<-function(y){
	N<-sum(y)
	ccr<-sum(diag(y))/sum(y)
	p<-apply(y,1,sum)/N
	q<-apply(y,2,sum)/N
	num<-ccr-sum(p*q)
	den<-1-sum(p*q)
	k<-num/den
	k[k<0]<-0
	return(k)		
	}

#Kappa optimization function
kappa.opt<-function(threshold){
	obs<-out[,1]>0
	if(threshold==0){
		pred<-out[,2]>=threshold
		}
	else {
		pred<-out[,2]>threshold
		}
	temp<-c(sum(pred&obs),sum(!pred&obs),sum(pred&!obs),sum(!pred&!obs))
	temp<-as.table(matrix(temp,nrow=2))
	cohen.kappa(temp)
	}
#########################################################################

#compute observed FP and cutpoint
out<-as.data.frame(cbind(fit$y,fit$fitted.values))
names(out)<-c('observed','fitted')
temp<-optimize(kappa.opt,interval=c(min(out$fitted),max(out$fitted)),maximum=TRUE)
cutpoint<-temp$maximum
Kappa.obs<-temp$objective

#null (random) distribution of Kappa
if(test==TRUE){
	Kappa.null<-NULL
	if(is.null(fit$na.action)==TRUE){
		z2<-fit$data
		}
	else{
		z2<-fit$data[-fit$na.action,]
		}
	nobs<-nrow(z2)
	for(i in 1:reps){
		z2[,c(as.character(fit$formula[[2]]))]<-z2[sample(1:nobs,replace=FALSE),
			c(as.character(fit$formula[[2]]))]
		e<-new.env()
		e$fit<-fit
		for(i in names(z2)) assign(i,z2[[i]],envir=e)
		temp<-glm(fit$formula,family=fit$family,offset=fit$offset,weights=fit$prior.weights,data=e)
		out<-as.data.frame(cbind(temp$y,temp$fitted.values))
		names(out)<-c('observed','fitted')
		temp<-optimize(kappa.opt,interval=c(min(out$fitted),max(out$fitted)),maximum=TRUE)
		Kappa.null<-c(Kappa.null,round(temp$objective,3))
		}

	#compute p-value
	p.value<-sum(Kappa.null>=Kappa.obs)
	p.value<-(1+p.value)/(1+reps) 
	}

if(test==TRUE & plot==TRUE){
	xmin<-min(Kappa.obs,min(Kappa.null))
	xmax<-max(Kappa.obs,max(Kappa.null))*1.2
	hist(Kappa.null,col='blue',xaxs='i',yaxs='i',xlim=c(xmin,xmax),
		xlab='Kappa',cex.lab=1.5,
		main='Randomization test of Kappa',...)
	abline(v=Kappa.obs,col='black',lty=2,lwd=2,...)
	legend('topright',inset=c(0,0),legend='Observed',lty=2,lwd=2,bty='n')
	legend('topright',inset=c(0,.05),legend='Null',fill='blue',bty='n')
	}

if(test==TRUE){
	z<-as.data.frame(rbind('Cutpoint'=cutpoint,'Observed Kappa'=round(Kappa.obs,3),
		'P-value'=format.pval(p.value,digits=3,eps=.001)))
	names(z)=''
	}
else{
	z<-as.data.frame(rbind('Cutpoint'=cutpoint,'Observed Kappa'=round(Kappa.obs,3)))
	names(z)=''
	}

return(z)
}
lda.structure <-
function(x.lda,x,dim=ncol(x.lda),
	digits=3,cutoff=0){

#check for dim limit
if(dim>ncol(x.lda)){
    cat("Only",ncol(x.lda),"axes available\n")
    dim<-length(x.lda)
	}

#calculate structure correlations
z<-cor(x,x.lda[,1:dim])

#print results 
cat("\nStructure Correlations:\n")
z<-round(z,digits=digits)
z[abs(z)<cutoff]<-substring('',1,nchar(z[1,1]))
z<-as.data.frame(z)
return(z)
}
lnorm.plot <-
function(meanlog=0,sdlog=1,xlim=c(0,15)){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(meanlog)*length(sdlog)
k<-0
for(i in meanlog){
	for(j in sdlog){
	  k<-k+1
		curve(plnorm(x,meanlog=i,sdlog=j),xlim[1],xlim[2],
			xlab='Z',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
			mtext(paste('(meanlog=',i,', sdlog=',j,')',sep=''),side=3,col='red')
		curve(dlnorm(x,meanlog=i,sdlog=j),xlim[1],xlim[2],
			xlab='Z',ylab='Probability density',main='Probability Density')
			mtext(paste('(meanlog=',i,', sdlog=',j,')',sep=''),side=3,col='red')
		curve(qlnorm(x,meanlog=i,sdlog=j),0,1,
			xlab='Cumulative probability (quantile)',ylab='Z',main='Quantiles')
			mtext(paste('(meanlog=',i,', sdlog=',j,')',sep=''),side=3,col='red')
		y<-rlnorm(1000,meanlog=i,sdlog=j)
		hist(y,xlab='Z',ylab='Frequency',main='Random Observations')
			mtext(paste('(meanlog=',i,', sdlog=',j,')',sep=''),side=3,col='red')
		if(k<n.plots) readline('press return for next plot')
		}
	}
}
mantel2 <-
function(xdis,ydis,method='pearson',permutations=1000,strata=NULL){

  xdis<-as.dist(xdis)
  ydis<-as.vector(as.dist(ydis))
  tmp<-cor.test(as.vector(xdis),ydis,method=method)
  statistic<-as.numeric(tmp$estimate)
  variant<-tmp$method
  if (permutations) {
    N<-attr(xdis,"Size")
    perm<-rep(0,permutations)
    xmat<-as.matrix(xdis)
    asdist<-row(xmat)>col(xmat)

    for(i in 1:permutations) {
      take<-vegan:::permuted.index(N,strata)
      permvec<-(xmat[take,take])[asdist]
      perm[i]<-cor(permvec,ydis,method=method)
      }
    signif<-(sum(perm>=statistic) + 1)/(permutations + 1)
    }

    else{
        signif<-NA
        perm<-NULL
	    }
    res<-list(call=match.call(),method=variant,statistic=statistic, 
        signif=signif,perm=perm,permutations=permutations,xdis=xdis,ydis=ydis)

  if(!missing(strata)){
     res$strata<-deparse(substitute(strata))
     res$stratum.values<-strata
	   }
  class(res)<-'mantel'
  res
}
mantel.part <-
function(y,x1,x2,x3='',p='', digits=3,
	ydist='euclidean',xdist='euclidean',pdist='euclidean',
	yscale=FALSE,xscale=TRUE,pscale=FALSE,...){

library(ecodist)

#optionally standardize data matrices
if(!class(y)=='dist' && yscale==TRUE) y<-scale(y)
if(!class(x1)=='dist' && xscale==TRUE) x1<-scale(x1)
if(!class(x2)=='dist' && xscale==TRUE) x2<-scale(x2)
if(!x3=='' && !class(x3)=='dist' && xscale==TRUE) x3<-scale(x3)
if(!p=='' && !class(p)=='dist' && pscale==TRUE) p<-scale(p)

#optionally create distance matrices
if(!class(y)=='dist') y<-data.dist(y,method=ydist,...)
if(!class(x1)=='dist') x1<-data.dist(x1,method=xdist,...)
if(!class(x2)=='dist') x2<-data.dist(x2,method=xdist,...)
if(!x3=='' && !class(x3)=='dist') x3<-data.dist(x3,method=xdist,...)
if(!p=='' && !class(p)=='dist') p<-data.dist(p,method=pdist,...)

#mantel tests
if(!p==''){
	z1<-mantel(y~x1+p,...)
	z2<-mantel(y~x2+p,...)
	z1.2<-mantel(y~x1+x2+p,...)
	z2.1<-mantel(y~x2+x1+p,...)
	if(!x3==''){
		z3<-mantel(y~x3+p,...)
		z1.3<-mantel(y~x1+x3+p,...)
		z1.23<-mantel(y~x1+x2+x3+p,...)
		z2.3<-mantel(y~x2+x3+p,...)
		z2.13<-mantel(y~x2+x1+x3+p,...)
		z3.1<-mantel(y~x3+x1+p,...)
		z3.2<-mantel(y~x3+x2+p,...)
		z3.12<-mantel(y~x3+x1+x2+p,...)
		}
	}
else{
	z1<-mantel(y~x1,...)
	z2<-mantel(y~x2,...)
	z1.2<-mantel(y~x1+x2,...)
	z2.1<-mantel(y~x2+x1,...)
	if(!x3==''){
		z3<-mantel(y~x3,...)
		z1.3<-mantel(y~x1+x3,...)
		z1.23<-mantel(y~x1+x2+x3,...)
		z2.3<-mantel(y~x2+x3,...)
		z2.13<-mantel(y~x2+x1+x3,...)
		z3.1<-mantel(y~x3+x1,...)
		z3.2<-mantel(y~x3+x2,...)
		z3.12<-mantel(y~x3+x1+x2,...)
		}
	}
		
out<-list()

if(x3==''){
	##partioning call summary
	cat("Call:\n")
	out$call<-match.call()
	print(out$call)

	##conditional table
	if(!p==''){
		out$ptable<-deparse(out$call[[5]])
		cat("\nConditional table (included in all models):\n")
		cat(out$ptable,'\n')
		}
	
	##list of explanatory tables
	mx<-rep('',2)
	for(i in 1:2) mx[i]<-deparse(out$call[[i+2]])
	out$xtables<-mx
	cat("\nExplanatory tables:\n")
	cat(paste(paste(paste("X", 1:length(mx), sep = ""), 
	    ":  ", mx, sep = ""), collapse = "\n"), "\n")

	##marginal Mantel tests
	cat('\nMarginal Mantel tests:','\n')
	    out2<-rbind(z1,z2)
	    out$marginal<-as.data.frame(out2)
	    print(round(out$marginal,digits=digits))
	
	##partial Mantel tests
	cat('\nPartial Mantel tests:','\n')
	    out3<-rbind(z1.2,z2.1)
	    out$partial<-as.data.frame(out3)
	    print(round(out$partial,digits=digits))
	}
	
else{
	##partioning call summary
	cat("Call:\n")
	out$call<-match.call()
	print(out$call)
	
	##conditional table
	if(!p==''){
		out$ptable<-deparse(out$call[[6]])
		cat("\nConditional table (included in all models):\n")
		cat(out$ptable,'\n')
		}

	##explanatory tables	
	mx<-rep('',3)
	for(i in 1:3) mx[i]<-deparse(out$call[[i+2]])
	out$xtables<-mx
	cat("\nExplanatory tables:\n")
	cat(paste(paste(paste("X", 1:length(mx), sep = ""), 
	    ":  ", mx, sep = ""), collapse = "\n"), "\n")

	##marginal Mantel tests
	cat('\nMarginal Mantel tests:','\n')
	    out2<-rbind(z1,z2,z3)
	    out$marginal<-as.data.frame(out2)
	    print(round(out$marginal,digits=digits))
	
	##partial Mantel tests
	cat('\nPartial Mantel tests:','\n')
	    out3<-rbind(z1.2,z1.3,z2.1,z2.3,
	    	z3.1,z3.2,z1.23,z2.13,z3.12)
	    out$partial<-as.data.frame(out3)
	    print(round(out$partial,digits=digits))
	}
	
class(out)<-'mantel.part'
out
}
mrpp2 <-
function(dat,grouping,permutations=1000,distance='euclidean',
  weight.type=1,strata){

  classmean<-function(ind,dmat,indls) {
    sapply(indls,function(x) mean(c(dmat[ind==x,ind==x]),na.rm = TRUE))
    }
  mrpp.perms<-function(ind,dmat,indls,w){
  	weighted.mean(classmean(ind,dmat,indls),w=w,na.rm=TRUE)
	  }
    
  if(inherits(dat,'dist')) 
    dmat<-dat
  else if(is.matrix(dat)&&nrow(dat)==ncol(dat)&&all(dat[lower.tri(dat)]==t(dat)[lower.tri(dat)])){
    dmat<-dat
    attr(dmat,'method')<-'user supplied square matrix'
	  }
  else dmat<-vegdist(dat,method=distance)
  if (any(dmat < -sqrt(.Machine$double.eps))) 
    stop("dissimilarities must be non-negative")
  
  distance<-attr(dmat,'method')
  dmat<-as.matrix(dmat)
  diag(dmat)<-NA
  N<-nrow(dmat)
  grouping<-factor(grouping)
  indls<-levels(grouping)
  ncl<-sapply(indls,function(x) sum(grouping==x))
  w<-switch(weight.type,ncl,ncl-1,ncl*(ncl-1)/2)
  classdel<-classmean(grouping,dmat,indls)
  names(classdel)<-names(ncl)<-indls
  del<-weighted.mean(classdel,w=w,na.rm=TRUE)
  E.del<-mean(dmat,na.rm=TRUE)
  CS<-NA
  if (missing(strata))
	  strata<-NULL
  perms<-sapply(1:permutations,function(x) grouping[vegan:::permuted.index(N,strata=strata)])
  m.ds<-numeric(permutations)
  m.ds<-apply(perms,2,function(x) mrpp.perms(x,dmat,indls,w))
  p<-(1+sum(del>=m.ds))/(permutations+1)
  r2<-1-del/E.del
    
  x<-as.dist(dmat)

  matched<-function(irow,icol,grouping){
	  grouping[irow]==grouping[icol]
		}
  x.vec<-as.vector(x)
	N<-attributes(x)$Size
	irow<-as.vector(as.dist(row(matrix(nrow=N,ncol=N))))
	icol<-as.vector(as.dist(col(matrix(nrow=N,ncol=N))))
	within<-matched(irow,icol,grouping)
	cl.vec<-rep("Between",length(x))
	take<-as.numeric(irow[within])
	cl.vec[within]<-levels(grouping)[grouping[take]]
	cl.vec<-factor(cl.vec,levels=c("Between",levels(grouping)))
    
  out<-list(call=match.call(),delta=del,E.delta=E.del, 
    CS=CS,n=ncl,classdelta=classdel,Pvalue=p,A=r2,distance=distance,
    weight.type=weight.type,boot.deltas=m.ds,permutations=permutations,
    class.vec=cl.vec,diss.vec=x.vec)
  if (!is.null(strata)){
    out$strata<-deparse(substitute(strata))
    out$stratum.values<-strata
		}
  class(out)<-'mrpp'
  out
}
mv.outliers <-
function(x,method,var='',cor.method='pearson',
	outfile='',sd.limit=3,alpha=.001,plot=TRUE,save.plot=FALSE,
	use='complete.obs',na.rm=TRUE,col.hist='blue',col.line='black',
	col.point='blue',las=1,lab=c(5,5,4),...){

library(vegan)
library(MASS)

old.par<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

#paired histogram function for avedist and sd(avedist)
paired.hist<-function(file=d5){
	par(mfrow=c(2,1),mar=c(5,5,4,2)) #graphics settings
	hist(file[[1]],prob=TRUE,col=col.hist,las=las,lab=lab,
		xaxs='i',yaxs='i',xlab='average distance',
		main=paste('Histogram of',method,'distance',sep=' '),...)
	par(new=TRUE)
	plot(density(file[[1]],na.rm=na.rm),
		xaxs='i',yaxs='i',col=col.line,las=las,lab=lab,
		axes=FALSE,xlab='',ylab='',main='',...)
	hist(file[[2]],prob=TRUE,col=col.hist,las=las,lab=lab,
		xaxs='i',yaxs='i',xlab='standard deviation of average distance',
		main='',...)
	par(new=TRUE)
	plot(density(file[[2]],na.rm=na.rm),
		xaxs='i',yaxs='i',col=col.line,las=las,lab=lab,
		axes=FALSE,xlab='',ylab='',main='',...)
	if(save.plot==TRUE){ #save plot to file
		dev.print(jpeg,file=paste('dhist.',method,'.jpg',sep=''),width=800,height=600)
		} #end save	plot
	} #end paired histogram function

#ecological distance calculations
if(method=='mahalanobis'){ #compute mahalanobis distance
	stopifnot(mahalanobis(y,0,diag(ncol(y)))==rowSums(y*y)) #Here,D^2=usual squared Euclidean distances
#	y.pca<-prcomp(y,center=TRUE,scale=TRUE)
#	Sy<-cov.rob(y.pca$x)
	Sy<-cov.rob(y)
	z<-mahalanobis(y,Sy$center,Sy$cov)
	n<-nrow(x)
	p<-ncol(x)
	critical.chi<-qchisq(1-alpha,p) # asymptotic chi-square method	 
	if(plot==TRUE){
		par(mfrow=c(2,1),mar=c(5,5,4,2)) #graphics settings
		plot(density(z),col=col.line,las=las,lab=lab,
			main=paste('Squared Mahalanobis distances', 
			'(n=',nrow(y),'p=',ncol(y),')',sep=' '),...)
		rug(z)
		abline(v=critical.chi,col='red',lty=2)
		qqplot(qchisq(ppoints(nrow(y)),df=ncol(y)),z,
			col=col.point,las=las,lab=lab,
			xlab='Theoretical Quantiles',ylab='Sample Quantiles',
			main=expression("Q-Q plot of Mahalanobis" * ~D^2 *
	        " vs. quantiles of" * ~ chi[3]^2),...)
	    abline(0,1,col='red',lty=2)
	    }
    z<-as.data.frame(z)
    names(z)<-c('mahalanobis distance')
	t1<-z>=critical.chi | z<=-sd.limit #identify extreme values
	t2<-as.matrix(t1)
	t3<-z[apply(t2,1,FUN='any',na.rm=na.rm),,drop=FALSE] #select rows with extremes
	z<-round(t3,3) #round results to 3 decimal places
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		} #end save outfile
 	} #end mahalanobis distance

else{ #all other distances 
	if(method=='correlation'){ #compute correlation distance
		d1<-as.matrix((1-cor(t(y),method=cor.method,use=use))/2) #compute distance
		} #end correlation distance
	else{ #compute other distance methods
		d1<-as.matrix(vegdist(y,method=method,na.rm=na.rm,...)) #vegdist function
		} #end other distance methods
	diag(d1)<-NA
	d2<-apply(d1,2,mean,na.rm=na.rm)
	d3<-(d2-mean(d2))/sd(d2)
	d4<-c('avedist','sd')
	d5<-cbind(d2,d3)
	colnames(d5)<-d4
	d5<-as.data.frame(d5)
	if(plot==TRUE){ #plot paired histograms
		paired.hist(d5)
		} #end plot paired histograms
	t1<-d3>=sd.limit | d3<=-sd.limit #identify extreme values
	t2<-as.matrix(t1)
	t3<-d5[apply(t2,1,FUN='any',na.rm=na.rm),,drop=FALSE] #select rows with extremes
	z<-round(t3,3) #round results to 3 decimal places
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		} #end save outfile
	} #end all other distances
	
par(old.par)
return(z)
}
nbinom.plot <-
function(events=25,mu=2,size=1){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(mu)*length(size)
k<-0
n<-seq(0,events)
for(i in mu){
	for(j in size){
	  k<-k+1
		barplot(pnbinom(n,mu=i,size=j),names=as.character(n),
			xlab='Count',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
			mtext(paste('(mu=',i,', k=',j,')',sep=''),side=3,col='red')
		barplot(dnbinom(n,mu=i,size=j),names=as.character(n),
			xlab='Count',ylab='Probability',main='Probability Mass')
			mtext(paste('(mu=',i,', k=',j,')',sep=''),side=3,col='red')
		barplot(qnbinom(seq(0,.95,.05),mu=i,size=j),names=seq(0,.95,.05),
			xlab='Cumulative probability (quantile)',ylab='Count',main='Quantiles')
			mtext(paste('(mu=',i,', k=',j,')',sep=''),side=3,col='red')
		y<-rnbinom(1000,mu=i,size=j)
	  y.table<-as.data.frame(table(y))
	  y.levels<-as.data.frame(n)
	  names(y.levels)<-'y'
	  y.table<-merge(y.levels,y.table,all.x=TRUE)
	  y.table$Freq[is.na(y.table$Freq)]<-0
	  barplot(y.table$Freq,names=as.character(n),xlab='Count',ylab='Frequency',
      main='Random Observations',col='gray')
			mtext(paste('(mu=',i,', k=',j,')',sep=''),side=3,col='red')
	  if(k<n.plots) readline('press return for next plot')
		}
	}
}
nhclus.scree <-
function(x,max.k,...){

library(cluster)
old.par<-par(no.readonly=TRUE)

#nonhierarchical clustering
s<-rep(0,max.k)
d<-rep(0,max.k)
for(i in 2:max.k){
	temp<-pam(x,k=i,...)
	s[i]<-temp$silinfo$avg.width #get result (ave silhouette width)
	d[i]<-temp$objective[2] #get result (min sum of diss to medoid)
	}

#merge results into data.frame
s<-s[-1]
d<-d[-1]
n<-c(2:max.k)
y<-round(as.data.frame(cbind(n,d,s)),3)
names(y)<-c('no. clusters','sum within-clus diss','ave silhouette width')

#create scree plot
par(mar=c(5.1,4.5,4.1,4.5))
plot(y[,1],y[,2],type='o',lwd=1.5,pch=19,col='blue',
	ylab='Within-cluster Dissimilarity',xlab='Number of Clusters')
par(new=TRUE)
plot(y[,1],y[,3],type='o',lty=2,lwd=1.5,pch=19,col='red',
	xlab='',ylab='',yaxt='n')
	sil.tick<-pretty(range(y[,3]))
	axis(side=4,at=sil.tick,srt=90)
	mtext('Average Silhouette Width',side=4,line=3)
legend(x='topright',inset=c(.05,.05),legend=c('Within-cluster dissimilarity','Average silhoutte width'),col=c('blue','red'),lty=c(1,3),lwd=c(1.5,1.5))
title(main='Scree Plot of K-means Clustering')

par(old.par)
return(y)
}
nmds.monte <-
function(x,k,distance='bray',trymax=50,autotransform=FALSE,
	trace=0,zerodist='add',perm=100,col.hist='blue',col.line='red',
	lty=2,las=1,lab=c(5,5,4),...){

library(vegan)
library(MASS)

z<-metaMDS(comm=x,k=k,distance=distance,trymax=trymax,
	autotransform=autotransform,trace=trace,...) #nmds analysis
z.stress<-z$stress #get stress
y.stress<-rep(0,perm)

for(i in 1:perm){
	y<-apply(x,2,sample) #permute data matrix
	y<-metaMDS(comm=y,k=k,distance=distance,trymax=trymax,
		autotransform=autotransform,trace=trace,...) #nmds analysis
	y.stress[i]<-y$stress #get stress
	}
n<-sum(y.stress<=z.stress) #compute number of random runs with stress < observed
p.value<-(1+n)/(1+perm) #compute p-value

xmin<-min(z.stress,min(y.stress))
xmax<-max(z.stress,max(y.stress))
hist(y.stress,col=col.hist,las=las,lab=lab,
	xaxs='i',yaxs='i',xlim=c(xmin,xmax),xlab='Stress',
	main=paste('Random Permutation Distribution of Stress for',k,'Dimensions',sep=' '),...)
abline(v=z.stress,col=col.line,lty=lty,lwd=2,...)

cat('Randomization Test of Stress:\n')
cat('Permutation stress values:\n')
print(y.stress)
z<-rbind('Observed stress'=z.stress,'P-value'=p.value)
return(z)
}
nmds.scree <-
function(x,distance='bray',k=6,trymax=50,
	autotransform=FALSE,trace=0,...){

library(vegan)
library(MASS)

old.par<-par(no.readonly=TRUE)

nmds.stress<-rep(0,k)
nmds.dim<-c(1:k)
	for(i in 1:k){
	y<-metaMDS(x,distance=distance,k=i,trymax=trymax,
		autotransform=autotransform,trace=trace,...)
	nmds.stress[i]<-y$stress
	}
plot(nmds.dim,nmds.stress,type='o',pch=19,col='blue',
	ylab='Stress',xlab='Ordination Axis',
	main='Scree Plot of Stress vs. Dimension',...)

par(old.par)
}
norm.plot <-
function(mean=0,sd=1,xlim=c(-3,3)){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(mean)*length(sd)
k<-0
for(i in mean){
	for(j in sd){
	  k<-k+1
		curve(pnorm(x,mean=i,sd=j),xlim[1],xlim[2],
			xlab='Z',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
			mtext(paste('(mean=',i,', sd=',j,')',sep=''),side=3,col='red')
		curve(dnorm(x,mean=i,sd=j),xlim[1],xlim[2],
			xlab='Z',ylab='Probability density',main='Probability Density')
			mtext(paste('(mean=',i,', sd=',j,')',sep=''),side=3,col='red')
		curve(qnorm(x,mean=i,sd=j),0,1,
			xlab='Cumulative probability (quantile)',ylab='Z',main='Quantiles')
			mtext(paste('(mean=',i,', sd=',j,')',sep=''),side=3,col='red')
		y<-rnorm(1000,mean=i,sd=j)
		hist(y,xlab='Z',ylab='Frequency',main='Random Observations')
			mtext(paste('(mean=',i,', sd=',j,')',sep=''),side=3,col='red')
	  if(k<n.plots) readline('press return for next plot')
		}
	}
}
norm.test <-
function(x,groups='',var='',method='ad',...){

library(nortest)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

if(!groups==''){
	y.resid<-matrix(0,nrow(y),ncol(y))
	grp<-as.factor(groups)
	for(i in 1:ncol(y)){ #loop thru variables
		y.aov<-aov(y[,i]~grp)
		y.resid[,i]<-y.aov$residuals
		}
	y.resid<-as.data.frame(y.resid)
	rownames(y.resid)<-rownames(y)
	colnames(y.resid)<-colnames(y)
	y<-y.resid
	}
	
#create summary table
z<-matrix(0,ncol(y),2) #create blank matrix
for(i in 1:ncol(y)){ #loop thru variables
	if(method=='ad'){
		temp<-ad.test(y[,i],...)
		}
	else if(method=='sf'){
		temp<-sf.test(y[,i],...)
		}		
	if(method=='cvm'){
		temp<-cvm.test(y[,i],...)
		}
	else if(method=='lillie'){
		temp<-lillie.test(y[,i],...)
		}
	else if(method=='pearson'){
		temp<-pearson.test(y[,i],...)
		}
	z[i,1]<-round(temp$statistic,3)
	z[i,2]<-temp$p.value
	}
rownames(z)<-colnames(y)
if(method=='ad'){
	colnames(z)<-c('Anderson-Darling A','p-value')
	cat('Anderson-Darling Test of Normality:\n')
	}
if(method=='sf'){
	colnames(z)<-c('Shapiro-Francia','p-value')
	cat('Shapiro-Francia Test of Normality:\n')
	}
if(method=='cvm'){
	colnames(z)<-c('Cramer-von Mises W','p-value')
	cat('Cramer-van Mises Test of Normality:\n')
	}
if(method=='lillie'){
	colnames(z)<-c('Lilliefors D','p-value')
	cat('Lilliefors (Kolmogorov-Smirnov) Test of Normality:\n')
	}
if(method=='pearson'){
	colnames(z)<-c('Pearson chi-square','p-value')
	cat('Pearson chi-square Test of Normality:\n')
	}
z<-as.data.frame(z)
z[,2]<-format.pval(z[,2],digits=3,eps=.001)		
return(z)
}
ordi.monte <-
function(x,ord,dim=length(x),perm=1000,center=TRUE,
	scale=TRUE,digits=3,plot=TRUE,col.hist='blue',col.line='red',
	lty=2,las=1,lab=c(5,5,4),...){

p<-length(x)
if(dim>p){
    cat("Only",p,"axes available\n")
    dim<-p
	}

if(ord=='pca'){
	z<-prcomp(x,center=center,scale=scale) #prcomp analysis
	z.eig<-z$sdev[1:dim]^2 #compute eigenvalues
	z.teig<-t(z.eig) #transpose eigenvalues
	z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
	write('',file='y.csv') #empty outfile if it exists
	for(i in 1:perm){
		y<-apply(x,2,sample) #permute data matrix
		y<-prcomp(y,center=center,scale=scale) #prcomp on permuted matrix
		y<-y$sdev[1:dim]^2 #compute eigenvalues
		y<-as.data.frame(t(y)) #coerce to data frame and transpose
		write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
		}
	y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
	p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
	p.value<-p.value/perm #compute p-value
	names<-colnames(z$rotation[,1:dim]) #add 'PC#' names
	}

else if(ord=='ca'){
	library(vegan)
	z<-cca(x) #correspondence analysis
	z.eig<-z$CA$eig[1:dim] #get eigenvalues
	z.teig<-t(z.eig) #transpose eigenvalues
	z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
	write('',file='y.csv') #empty outfile if it exists
	for(i in 1:perm){
		y<-apply(x,2,sample) #permute data matrix
		y<-cca(y) #CA on permuted matrix
		y<-y$CA$eig[1:dim] #get eigenvalues
		y<-as.data.frame(t(y)) #coerce to data frame and transpose
		write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
		}
	y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
	p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
	p.value<-p.value/perm #compute p-value
	names<-names(z$CA$eig[1:dim]) #add 'CA#' names
	}

else if(ord=='dca'){
	library(vegan)
	if(dim>4){
    cat("Only",4,"axes available\n")
    dim<-4
	}
	z<-decorana(x,...) #detrended correspondence analysis
	z.eig<-z$evals[1:dim] #get eigenvalues
	z.teig<-t(z.eig) #transpose eigenvalues
	z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
	write('',file='y.csv') #empty outfile if it exists
	for(i in 1:perm){
		y<-apply(x,2,sample) #permute data matrix
		y<-decorana(y,...) #DCA on permuted matrix
		y<-y$evals[1:dim] #get eigenvalues
		y<-as.data.frame(t(y)) #coerce to data frame and transpose
		write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
		}
	y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
	p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
	p.value<-p.value/perm #compute p-value
	names<-names(z$eval[1:dim]) #add 'CA#' names
	}

if(plot==TRUE){
	for(i in 1:dim){
		xmin<-min(min(y[[i]],z.eig[i]))
		xmax<-max(max(y[[i]],z.eig[i]))
		hist(y[[i]],col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlim=c(xmin,xmax),
			xlab='Eigenvalue',
			main=paste('Random Permutation Distribution of Eigenvalues for',names[i],sep=' '),...)
		abline(v=z.eig[i],col=col.line,lty=lty,lwd=2,...)
	    readline("Press return for next plot ")
		}
	}	

cat('Randomization Test of Eigenvalues:\n')
z<-rbind('Eigenvalue'=z.eig,'P-value'=p.value)
colnames(z)<-names
z<-round(z,digits=digits)
return(z)
}
ordi.overlay <-
function(x.ord,x,var='',fit=TRUE,choices=c(1,2),expand=5,
	tau=.95,pch=19,...){

library(cobs)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
scores<-scores(x.ord,display='site',choices=choices)
names<-colnames(scores)

for(i in 1:length(y)){
	if(fit==TRUE){
		set.panel(2,2,relax=TRUE)
	
		#top left plot
		out<-cobs(scores[,2],y[,i],tau=tau,...)
		plot(scores[,2],y[,i],xlab=names[2],ylab=names(y[i]),pch=pch)
		orderx<-order(scores[,2])
		lines(scores[,2][orderx],out$fitted[orderx],lty=1,col='blue')	
		
		#top right plot
		plot(scores[,1],scores[,2],col='blue',pch=pch,
		     xlab=names[1],ylab=names[2],
		     cex=expand*y[,i]/max(y[,i]),main=names(y[i]))
		
		#bottom left plot
		frame()
		
		#bottom right plot
		out<-cobs(scores[,1],y[,i],tau=tau,...)
		plot(scores[,1],y[,i],xlab=names[1],ylab=names(y[i]),pch=pch)
		orderx<-order(scores[,1])
		lines(scores[,1][orderx],out$fitted[orderx],lty=1,col='blue')	
		  
		readline("Press return for next plot ")    
		}
	else{
		plot(scores[,1],scores[,2],col='blue',pch=pch,
			xlab=names[1],ylab=names[2],
			cex=expand*y[,i]/max(y[,i]),main=names(y[i]))
		readline("Press return for next plot ")
		}
	}
}
ordi.part <-
function(y,x1,x2,x3='',p='',method='rda',model='reduced',digits=3,...){

#create datasets
x12<-cbind(x1,x2)
x1p<-cbind(x1,p)
x2p<-cbind(x2,p)
x12p<-cbind(x1,x2,p)
if(!x3==''){
	x13<-cbind(x1,x3)
	x23<-cbind(x2,x3)
	x123<-cbind(x1,x2,x3)
	x3p<-cbind(x3,p)
	x13p<-cbind(x1,x3,p)
	x23p<-cbind(x2,x3,p)
	}
##ordination models
if(method=='rda'){
	if(!p==''){
		z1<-rda(y,x1,p,...)
		z2<-rda(y,x2,p,...)
		z12<-rda(y,x12,p,...)
		z1.2<-rda(y,x1,x2p,...)
		z2.1<-rda(y,x2,x1p,...)
		if(!x3==''){
			z3<-rda(y,x3,p,...)
			z13<-rda(y,x13,p,...)
			z23<-rda(y,x23,p,...)
			z1.3<-rda(y,x1,x3p,...)
			z1.23<-rda(y,x1,x23p,...)
			z2.3<-rda(y,x2,x3p,...)
			z2.13<-rda(y,x2,x13p,...)
			z3.1<-rda(y,x3,x1p,...)
			z3.2<-rda(y,x3,x2p,...)
			z3.12<-rda(y,x3,x12p,...)
			z123<-rda(y,x123,p,...)
			}
		}
	else{
		z1<-rda(y,x1,...)
		z2<-rda(y,x2,...)
		z12<-rda(y,x12,...)
		z1.2<-rda(y,x1,x2,...)
		z2.1<-rda(y,x2,x1,...)
		if(!x3==''){
			z3<-rda(y,x3,...)
			z13<-rda(y,x13,...)
			z23<-rda(y,x23,...)
			z1.3<-rda(y,x1,x3,...)
			z1.23<-rda(y,x1,x23,...)
			z2.3<-rda(y,x2,x3,...)
			z2.13<-rda(y,x2,x13,...)
			z3.1<-rda(y,x3,x1,...)
			z3.2<-rda(y,x3,x2,...)
			z3.12<-rda(y,x3,x12,...)
			z123<-rda(y,x123,...)
			}
		}
	}
else if(method=='cca'){
	if(!p==''){
		z1<-cca(y,x1,p,...)
		z2<-cca(y,x2,p,...)
		z12<-cca(y,x12,p,...)
		z1.2<-cca(y,x1,x2p,...)
		z2.1<-cca(y,x2,x1p,...)
		if(!x3==''){
			z3<-cca(y,x3,p,...)
			z13<-cca(y,x13,p,...)
			z23<-cca(y,x23,p,...)
			z1.3<-cca(y,x1,x3p,...)
			z1.23<-cca(y,x1,x23p,...)
			z2.3<-cca(y,x2,x3p,...)
			z2.13<-cca(y,x2,x13p,...)
			z3.1<-cca(y,x3,x1p,...)
			z3.2<-cca(y,x3,x2p,...)
			z3.12<-cca(y,x3,x12p,...)
			z123<-cca(y,x123,p,...)
			}
		}
	else{
		z1<-cca(y,x1,...)
		z2<-cca(y,x2,...)
		z12<-cca(y,x12,...)
		z1.2<-cca(y,x1,x2,...)
		z2.1<-cca(y,x2,x1,...)
		if(!x3==''){
			z3<-cca(y,x3,...)
			z13<-cca(y,x13,...)
			z23<-cca(y,x23,...)
			z1.3<-cca(y,x1,x3,...)
			z1.23<-cca(y,x1,x23,...)
			z2.3<-cca(y,x2,x3,...)
			z2.13<-cca(y,x2,x13,...)
			z3.1<-cca(y,x3,x1,...)
			z3.2<-cca(y,x3,x2,...)
			z3.12<-cca(y,x3,x12,...)
			z123<-cca(y,x123,...)
			}
		}
	}

##permutation tests
if(x3==''){
	x1.pv<-anova(z1,model=model,...)[1,5]
	x2.pv<-anova(z2,model=model,...)[1,5]
	v1.pv<-anova(z1.2,model=model,...)[1,5]
	v2.pv<-anova(z2.1,model=model,...)[1,5]
	}
else{
	x1.pv<-anova(z1,model=model,...)[1,5]
	x2.pv<-anova(z2,model=model,...)[1,5]
	x3.pv<-anova(z3,model=model,...)[1,5]
	x12.pv<-anova(z12,model=model,...)[1,5]
	x13.pv<-anova(z13,model=model,...)[1,5]
	x23.pv<-anova(z23,model=model,...)[1,5]
	v1.pv<-anova(z1.23,model=model,...)[1,5]
	v2.pv<-anova(z2.13,model=model,...)[1,5]
	v3.pv<-anova(z3.12,model=model,...)[1,5]
	}
		
out<-list()

if(x3==''){
	##partioning call
	cat("Call:\n")
	out$call<-match.call()
	print(out$call)

	##conditional table
	if(!p==''){
		out$ptable<-deparse(out$call[[5]])
		cat("\nConditional table (included in all models):\n")
		cat(out$ptable,'\n')
		}
	
	##list of explanatory tables
	mx<-rep('',2)
	for(i in 1:2) mx[i]<-deparse(out$call[[i+2]])
	out$xtables<-mx
	cat("\nExplanatory tables:\n")
	cat(paste(paste(paste("X", 1:length(mx), sep = ""), 
	    ":  ", mx, sep = ""), collapse = "\n"), "\n")

	##total effects
	cat('\nPartitioning of ',z12$inertia,':\n',sep = '')
	    out1<-c(Total=z12$tot.chi,
	     	Conditioned=z12$pCCA$tot.chi,
	        Constrained=z12$CCA$tot.chi,
	        Unconstrained=z12$CA$tot.chi)
	    out$total<-as.data.frame(cbind(Inertia=out1,Proportion=out1/out1[1]))
	    print(round(out$total,digits=digits))
	
	##marginal effects
	cat('\nMarginal effects:','\n')
	    out2<-c(x1=z1$CCA$tot.chi,
	    x2=z2$CCA$tot.chi)
	    out$marginal<-as.data.frame(cbind(Inertia=out2,
	    	Prop.Total=out2/z12$tot.chi,
	    	Prop.Constrained=out2/z12$CCA$tot.chi,
	    	Pvalue=c(x1.pv,x2.pv)))
	    print(round(out$marginal,digits=digits))
	
	##component effects
	#partition components
	cat('\nComponents:','\n')
	    out3<-c(v1=z1.2$CCA$tot.chi,
	    v2=z2.1$CCA$tot.chi,
	    v12=z1$CCA$tot.chi-z1.2$CCA$tot.chi)
	    out$components<-as.data.frame(cbind(Inertia=out3,
	    	Prop.Total=out3/z12$tot.chi,
	    	Prop.Constrained=out3/z12$CCA$tot.chi,
	    	Pvalue=c(v1.pv,v2.pv,NA)))
	    print(round(out$components,digits=digits))
	#other components
	if(!p=='') out$conditional<-c(z12$pCCA$tot.chi,z12$pCCA$tot.chi/z12$tot.chi)
	out$residual<-c(z12$CA$tot.chi,z12$CA$tot.chi/z12$tot.chi)
	}
	
else{
	##partioning call summary
	cat("Call:\n")
	out$call<-match.call()
	print(out$call)

	##conditional table
	if(!p==''){
		out$ptable<-deparse(out$call[[6]])
		cat("\nConditional table (included in all modesl):\n")
		cat(out$ptable,'\n')
		}

	##explanatory tables	
	mx<-rep('',3)
	for(i in 1:3) mx[i]<-deparse(out$call[[i+2]])
	out$xtables<-mx
	cat("\nExplanatory tables:\n")
	cat(paste(paste(paste("X", 1:length(mx), sep = ""), 
	    ":  ", mx, sep = ""), collapse = "\n"), "\n")

	##total effects
	cat('\nPartitioning of ',z123$inertia,':\n',sep = '')
	    out1<-c(Total=z123$tot.chi,
	     	Conditioned=z123$pCCA$tot.chi,
	        Constrained=z123$CCA$tot.chi,
	        Unconstrained=z123$CA$tot.chi)
	    out$total<-as.data.frame(cbind(Inertia=out1,Proportion=out1/out1[1]))
	    print(round(out$total,digits=digits))
	
	##marginal effects
	cat('\nMarginal effects:','\n')
	    out2<-c(x1=z1$CCA$tot.chi,
	    x2=z2$CCA$tot.chi,
	    x3=z3$CCA$tot.chi,
	    x12=z12$CCA$tot.chi,
	    x13=z13$CCA$tot.chi,
	    x23=z23$CCA$tot.chi)
	    out$marginal<-as.data.frame(cbind(Inertia=out2,
	    	Prop.Total=out2/z123$tot.chi,
	    	Prop.Constrained=out2/z123$CCA$tot.chi,
	    	Pvalue=c(x1.pv,x2.pv,x3.pv,x12.pv,x13.pv,x23.pv)))
	    print(round(out$marginal,digits=digits))
	
	##component effects
	#calculations
	v123<-z1.23$CCA$tot.chi+(z1$CCA$tot.chi-z1.2$CCA$tot.chi)+(z1$CCA$tot.chi-z1.3$CCA$tot.chi)-z1$CCA$tot.chi
	#partition components
	cat('\nComponents:','\n')
	    out3<-c(v1=z1.23$CCA$tot.chi,
	    v2=z2.13$CCA$tot.chi,
	    v3=z3.12$CCA$tot.chi,
	    v12=z1$CCA$tot.chi-z1.2$CCA$tot.chi-v123,
	    v13=z1$CCA$tot.chi-z1.3$CCA$tot.chi-v123,
	    v23=z2$CCA$tot.chi-z2.3$CCA$tot.chi-v123,
	    v123=v123)
	    out$components<-as.data.frame(cbind(Inertia=out3,
	    	Prop.Total=out3/z123$tot.chi,
	    	Prop.Constrained=out3/z123$CCA$tot.chi,
	    	Pvalue=c(v1.pv,v2.pv,v3.pv,NA,NA,NA,NA)))
	    print(round(out$components,digits=digits))
	#other components
	if(!p=='') out$conditional<-c(z123$pCCA$tot.chi,z123$pCCA$tot.chi/z123$tot.chi)
	out$residual<-c(z123$CA$tot.chi,z123$CA$tot.chi/z123$tot.chi)
	}
	
class(out)<-'ordi.part'
out
}
ordi.scree <-
function(x,ord,...){

library(vegan)
library(MASS)

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

if(ord=='pca'){
	if(class(x)=='prcomp'||class(x)=='princomp'){
		eig<-x$sdev^2
		if(!x$scale==''){
			inertia<-'correlations'
			}
		else{
			inertia<-'variance'
			}
		}
	else if(any(class(x)=='rda')){ #from rda()
		eig<-x$CA$eig
		inertia==x$inertia
		}
	#calculate broken stick values for PCA
	if(inertia=='correlations'){
		p<-length(eig)
		y<-rep(0,p)
		for(i in 1:p) y[i]<-1/i
		for(i in 1:p) y[i]<-sum(y[i:p]) #broken stick values
		}
	}
else if(ord=='ca'){ #from cca()
	eig<-x$CA$eig
	}
else if(ord=='mds'){ #from cmdscale()
	eig<-x$eig
	}			
else if(ord=='rda'||ord=='cca'||ord=='cmds'){ #from rda(),cca(),or capscale()
	eig<-x$CCA$eig
	}
	
#create scree plot
if(ord=='pca'&&inertia=="correlations"){ #PCA with broken stick
	plot(eig,type='o',lwd=1.5,pch=19,col='blue',ylim=c(0,max(max(eig),max(y))),
		ylab='Eigenvalue',xlab='Ordination Axis',
		main='Scree Plot of Unconstrained Eigenvalues',...)
	par(new=TRUE)
	plot(y,type='o',lty=3,lwd=1.5,pch=19,col='green',ylim=c(0,max(max(eig),max(y))),
		xlab='',ylab='',main='')
	abline(h=sum(eig)/p,col='red')
	legend(x='topright',inset=c(.1,.1),legend=c('Observed','Broken-stick','Latent root'),col=c('blue','green','red'),lty=c(1,3,1))
	}
else{ #all others
	plot(eig,type='o',lwd=1.5,pch=19,col='blue',ylim=c(0,max(eig)),
		ylab='Eigenvalue',xlab='Ordination Axis',...)
		if(ord=='pca'||ord=='ca'||ord=='mds'){
			title(main='Scree Plot of Unconstrained Eigenvalues')
			}
		else if(ord=='rda'||ord=='cca'||ord=='cmds'){
			title(main='Scree Plot of Constrained Eigenvalues')
			}
	}
#create cumulative scree plot	
readline("Press return for next plot ")	
plot(cumsum(eig/sum(eig)),type='o',pch=19,col='blue',
	ylab='Cumulative Proportion',xlab='Ordination Axis',...)
	if(ord=='pca'||ord=='ca'||ord=='mds'){
		title(main='Cumulative Scree Plot of Unconstrained Eigenvalues')
		}
	else if(ord=='rda'||ord=='cca'||ord=='cmds'){
		title(main='Cumulative Scree Plot of Constrained Eigenvalues')
		}	
}
panel.cor <-
function(x, y, digits=2, prefix="", cex.cor,...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x,y,...)
    r.abs <- abs(cor(x,y,...))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r.abs)
}
pbinom.plot <-
function(size=10,prob=.5,col=NULL,...){
legend<-NULL
out<-NULL
n<-seq(0,size)
for(j in prob){
	t1<-pbinom(n,size=size,prob=j)
	out<-cbind(out,t1)
	t2<-paste('size=',size,', prob=',j,sep='')
	legend<-rbind(legend,t2)
	}
if(is.null(col)) col.bar=seq(1:length(prob)) 
else col.bar=col
barplot(t(out),beside=TRUE,names=as.character(n),col=col.bar,
	xlab='# successes',ylab='probability',
  main='Cumulative Probability',
  cex.main=1.5,cex.lab=1.5,...)
abline(h=0,col='gray')
legend(x='topleft',inset=c(.01,.01),legend=legend,fill=col.bar,cex=1,bty='n')
}
pca.communality <-
function(x.pca,x,dim=length(x.pca$sdev),
	digits=3){
#don't know why communality for all dimensions doesn't always
#equal one, perhaps singularity issue.

#check for dim limit
if(dim>length(x.pca$sdev)){
    cat("Only",length(x.pca$sdev),"axes available\n")
    dim<-length(x.pca$sdev)
	}

#calculate final communality estimates
z<-cor(x,x.pca$x[,1:dim])
z<-z^2
c<-apply(z,1,sum)

#print results 
cat("\nFinal Communalities -",dim,'dimensions:\n',sep=' ')
c<-round(c,digits=digits)
c<-as.data.frame(c)
return(c)
}
pca.eigenval <-
function(x.pca,dim=length(x.pca$sdev),digits=7){

#check for dim limit
if(dim>length(x.pca$sdev)){
    cat("Only",length(x.pca$sdev),"axes available\n")
    dim<-length(x.pca$sdev)
	}

#calculate some variables
names<-colnames(x.pca$rotation[,1:dim])
var<-x.pca$sdev^2
trace<-sum(var)
prop.var<-var/trace

#broken-stick distribution
p<-length(x.pca$sdev)
y<-rep(0,p)
for(i in 1:p) y[i]<-1/i
for(i in 1:p) y[i]<-sum(y[i:p])
y<-y[1:dim]

#print results    
cat('Importance of components:\n')
z<-rbind('Variance(eigenvalue)'=var[1:dim],
	'Proportion of Variance'=prop.var[1:dim],
	'Cumulative Proportion'=cumsum(prop.var[1:dim]),
	'Broken-stick value'=y)
colnames(z)<-names
z<-round(z,digits=digits)
return(z)
}
pca.eigenvec <-
function(x.pca,dim=length(x.pca$sdev),
	digits=7,cutoff=0){

#check for dim limit	
if(dim>ncol(x.pca$rotation)){
    cat("Only",ncol(x.pca$rotation),"axes available\n")
    dim<-ncol(x.pca$rotation)
    }

#print results    
cat("\nEigenvectors:\n")
z<-format(round(x.pca$rotation[,1:dim],digits=digits))
z[abs(x.pca$rotation[,1:dim])<cutoff]<-substring('',1,nchar(z[1,1]))
z<-as.data.frame(z)
return(z)
}
pca.structure <-
function(x.pca,x,dim=length(x.pca$sdev),
	digits=3,cutoff=0){

#check for dim limit
if(dim>length(x.pca$sdev)){
    cat("Only",length(x.pca$sdev),"axes available\n")
    dim<-length(x.pca$sdev)
	}

#calculate structure correlations
z<-cor(x,x.pca$x[,1:dim])

#print results 
cat("\nStructure Correlations:\n")
z<-round(z,digits=digits)
z[abs(z)<cutoff]<-substring('',1,nchar(z[1,1]))
z<-as.data.frame(z)
return(z)
}
plot.anosim <-
function(x,title1='ANOSIM (within- vs between-group rank dissimilarities)',
	title2='ANOSIM (observed vs expected R)',col='blue',...){

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

#boxplot between and within-group dissimilarities
boxplot(x$dis.rank~x$class.vec,notch=TRUE,varwidth=TRUE,col=col,
	ylab='Rank Dissimilarity',...)
title(title1)
if(x$permutations){
	pval<-format.pval(x$signif,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste("R = ",round(x$statistic,3),', ','P = ',pval),3)

readline("Press return for next plot ")

#histogram of permutation distribution
r<-x$statistic
E.r<-mean(x$perm)
perm<-x$perm
hist(perm,col=col,xlim=range(c(r,perm)),xlab='R',main='')
abline(v=r,col='red',lwd=2,lty=3)
title(title2)
if(x$permutations){
	pval<-format.pval(x$signif,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste('Observed R = ',round(x$statistic,3),',',
	'Expected R = ',round(E.r,3),',','P = ',pval),3)

}
plot.mantel <-
function(x,title1='MANTEL Scatterplot)',
	title2='MANTEL (observed vs expected R)',col='blue',...){

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

#setup graphics device
par(mfrow=c(2,1))

#scatterplot of dissimilarity matrices
plot(x$xdis,x$ydis,col=col,xlab='X-matrix dissimilarities',
	ylab='Y-matrix dissimilarities',...)
title(title1)
if(x$permutations){
	pval<-format.pval(x$signif,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste("R = ",round(x$statistic,3),', ','P = ',pval),3)

#histogram of permutation distribution
r<-x$statistic
E.r<-mean(x$perm)
perm<-x$perm
hist(perm,col=col,xlim=range(c(r,perm)),xlab='r',main='')
abline(v=r,col='red',lwd=2,lty=3)
title(title2)
if(x$permutations){
	pval<-format.pval(x$signif,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste('Observed r = ',round(x$statistic,3),',',
	'Expected r = ',round(E.r,3),',','P = ',pval),3)

}
plot.mrpp <-
function(x,title1='MRPP (within- vs between-group dissimilarities)',
	title2='MRPP (observed vs expected delta)',col='blue',...){

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

#boxplot of between and within-group dissimilarities
boxplot(x$diss.vec~x$class.vec,notch=TRUE,varwidth=TRUE,col=col,
	ylab='Dissimilarity',...)
title(title1)
if(x$permutations){
	pval<-format.pval(x$Pvalue,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste('A = ',round(x$A,3),',','P = ',pval),3)

readline("Press return for next plot ")

#histogram of permutation distribution
d<-x$delta
bd<-x$boot.deltas
hist(bd,col=col,xlim=range(c(d,bd)),xlab='Delta',main='')
abline(v=d,col='red',lwd=2,lty=3)
title(title2)
if(x$permutations){
	pval<-format.pval(x$Pvalue,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste('Observed delta = ',round(x$delta,3),',',
	'Expected delta = ',round(x$E.delta,3),',','P = ',pval),3)

}
plot.ordi.part <-
function(x,which='total',digits=1,...){

#plot template
parts<-length(x$xtables)
rad<-0.725
cp<-switch(parts, c(0, 0), c(0, 0, 1, 0), c(0, 0, 1, 0, 
    0.5, -sqrt(3/4)), c(-0.5, 0.3, 0.5, 0.3, 0, -sqrt(3/4) + 
    0.3))
cp<-matrix(cp, ncol = 2, byrow = TRUE)
plot(cp, axes = FALSE, xlab = "", ylab = "", asp = 1, type = "n", 
    xlim = (range(cp[, 1]) + c(-rad, rad)), ylim = (range(cp[, 
        2]) + c(-rad, rad)))
box()
symbols(cp, circles = rep(rad, min(parts, 3)), inches = FALSE, 
    add = TRUE, ...)

#components
if(which=='total') labels<-round(x$components[,2]*100,digits)
else if (which=='constrained') labels<-round(x$components[,3]*100,digits)
labels<-as.character(labels)
switch(parts, text(0, 0, labels, ...), text(rbind(cp[1,], 
	cp[2, ], colMeans(cp)), labels, ...), text(rbind(cp, 
    colMeans(cp[1:2, ]), colMeans(cp[c(1,3), ]), colMeans(cp[2:3, ]),
    colMeans(cp)), labels, ...))

#partition labels
labels<-x$xtables
labels<-as.character(labels)
switch(parts, text(0, .5, labels, ...), text(rbind(c(0,.5), 
	c(1,.5)), labels, ...), text(rbind(c(0,.5),c(1,.5),c(.5,-1.3)),
	labels, ...))

#other components and labels
xy<-par('usr')
if(nrow(x$total)==4){
	text(xy[2]-0.7*diff(xy[1:2]), xy[3]+0.05*diff(xy[3:4]), 
		paste('Conditional =', round(x$conditional[2]*100,digits)), pos=2, ...)
	}
text(xy[2]-0.05*diff(xy[1:2]), xy[3]+0.05*diff(xy[3:4]), 
	paste('Residual =', round(x$residual[2]*100,digits)), pos=2, ...)

#plot title
if(which=='total') title('Partition = Percent of Total Variance')
else if (which=='constrained') title('Partition = Percent of Explained Variance')
invisible()
}
pois.plot <-
function(events=25,lambda=1,xlim=c(0,25)){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(lambda)
k<-0
n<-seq(0,events)
for(i in lambda){
  k<-k+1
	barplot(ppois(n,lambda=i),names=as.character(n),
		xlab='# Events',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
		mtext(paste('(lambda=',i,')',sep=''),side=3,col='red')
	barplot(dpois(n,lambda=i),names=as.character(n),
		xlab='# Events',ylab='Probability',main='Probability Mass')
		mtext(paste('(lambda=',i,')',sep=''),side=3,col='red')
	barplot(qpois(seq(0,.95,.05),lambda=i),names=seq(0,.95,.05),
		xlab='Cumulative probability (quantile)',ylab='# Events',
    main='Quantiles')
		mtext(paste('(lambda=',i,')',sep=''),side=3,col='red')
	y<-rpois(1000,lambda=i)
  y.table<-as.data.frame(table(y))
  y.levels<-as.data.frame(n)
  names(y.levels)<-'y'
  y.table<-merge(y.levels,y.table,all.x=TRUE)
  y.table$Freq[is.na(y.table$Freq)]<-0
  barplot(y.table$Freq,names=as.character(n),xlab='# Events',ylab='Frequency',
    main='Random Observations',col='gray')
    mtext(paste('(lambda=',i,')',sep=''),side=3,col='red')
  if(k<n.plots) readline('press return for next plot')
	}
}
qqnorm.plots <-
function(x,var='',by='',save.plot=FALSE,
	na.rm=TRUE,col.point='blue',col.line='red',las=1,...){

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

if(by==''){ #QQnorm plots w/o groups
	for(i in 1:ncol(y)){ #loop thru variables
		qqnorm(y[,i],datax=TRUE,col=col.point,las=las,
			main=paste('Normal Q-Q Plot of',names(y[i]),sep=' '),...)
		y.IQR<-IQR(y[,i],na.rm=na.rm)
		if(y.IQR>0){ #add qqline
			qqline(y[,i],datax=TRUE,col=col.line,...)
			} #end qqline
		if(save.plot==TRUE){ #save plot
			dev.print(jpeg,file=paste('qqnorm.',names(y[i]),'.jpg',sep=''),width=800,height=600)
			} #end save plot
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end qqnorm w/o groups
	
else{ #QQnorm plots w/ groups
	n<-by.names(x,by) #create by variable
	y<-cbind(n,y) #bind with selected variables
	groups<-levels(y[,2]) #create vector of group names
	s<-floor(length(groups)^.5) #create multi-figure dimension
	s<-c(s,ceiling(length(groups)/s)) #create multi-figure dimensions
	for(i in 3:ncol(y)){ #loop thru variables
		par(mfrow=s,mar=c(5,5,4,2))	#graphics settings	
		for(j in 1:length(groups)){ #loop thru groups
			z<-y[y[,2]==groups[j],] #select records for group j
			qqnorm(z[,i],datax=TRUE,col=col.point,las=las,
				main=groups[j],
				ylab=paste('Sample quantiles of',names(z[i]),sep=' '),...)
			z.IQR<-IQR(z[,i],na.rm=na.rm)
			if(z.IQR>0){ #add qqline
				qqline(z[,i],datax=TRUE,col=col.line,...)
				} #end qqline
			if(j==length(groups) & save.plot==TRUE){
				dev.print(jpeg,file=paste('qqnorm.',names(z[i]),'.jpg',sep=''),width=800,height=600)
				} #end save
			} #end loop thru groups
			if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end qqnorm w/ groups
}
ran.split <-
function(x,grouping='',prop=.5){

if(!grouping==''){
	y<-cbind(grouping,x)
	}
else{
	y<-x #groups assumed to be in first column
	}
	
N<-nrow(y)
g<-runif(N)
g[g<prop]<-0
g[g>=prop]<-1
y<-cbind(g,y)
y<-split(y,g)
y1<-y[[1]]
y2<-y[[2]]
training.data<<-y1[,-c(1,2)]
validate.data<<-y2[,-c(1,2)]
training.grp<<-y1[,2]
validate.grp<<-y2[,2]

z1<-c(nrow(y1),nrow(y2))
z2<-round(c(nrow(y1)/N,nrow(y2)/N),2)
z1<-cbind(z1,z2)
colnames(z1)<-c('Number of samples','Proportion')
rownames(z1)<-c('Training set','Validation set')
z2<-table(training.grp)
z3<-table(validate.grp)
z<-list(z1,z2,z3,training.data,validate.data)
names(z)<-c('Random Subset Summary:','Training Table','Validation Table',
  'Training dataset','Validation dataset')
return(z)
}
redun.plot <-
function(x,var='',perm=1000,quantiles=c(.025,.975),...){

old.par<-par(no.readonly=TRUE)

if(!var==''){

	y<-subset(x,select=eval(parse(text=var))) #select variables to plot
	z.obs<-as.matrix(cor(y,use='complete.obs'))
	diag(z.obs)<-NA

	for(i in 1:ncol(y)){ #loop thru variables

		y.obs<-sort(z.obs[,i])
		y.ran<-matrix(0,perm,length(y.obs))
		for(j in 1:perm){
			t1<-apply(y,2,sample) #permute data matrix
			t2<-as.matrix(cor(t1,use='complete.obs'))
			diag(t2)<-NA
			t3<-sort(t2[,i])
			if(length(t3)==length(y.obs)) y.ran[j,]<-t3
			}
		y.ran.lower<-apply(y.ran,2,quantile,probs=quantiles[1])
		y.ran.upper<-apply(y.ran,2,quantile,probs=quantiles[2])

		par(new=FALSE)
		plot(y.obs,type='l',lwd=2,col='blue',
			xaxs='i',yaxs='i',ylim=c(-1,1),
			xlab='Rank order of pairwise correlations',
			ylab='Correlation',
			main=paste('Redundancy of ',names(y[i]),' vs. random data',sep=''),...)
		par(new=TRUE)
		plot(y.ran.lower,type='l',lwd=2,col='green',
			xaxs='i',yaxs='i',ylim=c(-1,1),
			xlab='',ylab='',main='',...)
		par(new=TRUE)
		plot(y.ran.upper,type='l',lwd=2,col='green',
			xaxs='i',yaxs='i',ylim=c(-1,1),
			xlab='',ylab='',main='',...)
		abline(h=0,col='red',lty=3,lwd=1,...)
		legend(x='bottomright',inset=c(.1,.1),legend=c('Actual','Random (95CI)'),
			col=c('blue','green'),lty=c(1,1),lwd=c(2,1))
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	}

else{
	z.obs<-sort(as.vector(as.dist(cor(x,use='complete.obs'))))
	z.ran<-matrix(0,perm,length(z.obs))
	for(i in 1:perm){
		t1<-apply(x,2,sample) #permute data matrix
		t2<-sort(as.vector(as.dist(cor(t1,use='complete.obs'))))
		if(length(t2)==length(z.obs)) z.ran[i,]<-t2
		}
	z.ran.lower<-apply(z.ran,2,quantile,probs=quantiles[1])
	z.ran.upper<-apply(z.ran,2,quantile,probs=quantiles[2])
	plot(z.obs,type='l',lwd=2,col='blue',
		xaxs='i',yaxs='i',ylim=c(-1,1),
		xlab='Rank order of pairwise correlations',
		ylab='Correlation',
		main='Redundancy of actual vs. random data',...)
	par(new=TRUE)
	plot(z.ran.lower,type='l',lwd=1,col='green',
		xaxs='i',yaxs='i',ylim=c(-1,1),
		xlab='',ylab='',main='',...)
	par(new=TRUE)
	plot(z.ran.upper,type='l',lwd=1,col='green',
		xaxs='i',yaxs='i',ylim=c(-1,1),
		xlab='',ylab='',main='',...)
	abline(h=0,col='red',lty=3,lwd=1,...)
	legend(x='bottomright',inset=c(.1,.1),legend=c('Actual','Random (95CI)'),
		col=c('blue','green'),lty=c(1,1),lwd=c(2,1))
	}

par(old.par)
}
replace.missing <-
function(x,var='',method='median',outfile=''){

if(!var==''){
	y1<-subset(x,select=-eval(parse(text=var))) #select all other variables
	y2<-subset(x,select=eval(parse(text=var))) #select all specified variables
	}
else{
	y2<-x #original variables
	}

if(method=='mean'){ #for method=mean
	for(i in names(y2)){ #loop through selected variables
		t<-round(mean(y2[[i]],na.rm=TRUE),3) #compute mean for each variable
		y2[is.na(y2[[i]]),i]<-t #assign mean value to missing value
		} 
	}
if(method=='median'){ #for method=median
	for(i in names(y2)){ #loop through selected variables
		t<-median(y2[[i]],na.rm=TRUE) #compute median for each variable
		y2[is.na(y2[[i]]),i]<-t #assign median value to missing value
		} 
	}

if(!var==''){
	z<-cbind(y1,y2) #combine deselected and (modified) selected variables
	}
else{z<-y2}
if(!outfile==''){ #save outfile
	write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
	} #end save outfile
return(z)
}
scatter.plots <-
function(data,y,x,fit='lowess',col.line='red',cex.main=2,...){

oldpar<-par(no.readonly=TRUE)

y<-as.data.frame(subset(data,select=eval(parse(text=y)))) #select y variable
x<-as.data.frame(subset(data,select=eval(parse(text=x)))) #select y variable
 
#loop thru y variables
for(i in 1:ncol(y)){ 

	#loop thru x variables
	for (j in 1:ncol(x)){
	
		#scatter plot
		plot(x[,j],y[,i],
			xlab=names(x[j]),ylab=names(y[i]),
			main=paste(names(x[j]),'vs',names(y[i]),sep='  '),
			cex.main=cex.main)
		if(fit=='lowess'){
			ok<-is.finite(x[,j]) & is.finite(y[,i]) 
 			if(any(ok))	lines(lowess(x[,j][ok],y[,i][ok]),col=col.line)
			}
		else if(fit=='linear'){
			ok<-is.finite(x[,j]) & is.finite(y[,i]) 
 			if(any(ok))	abline(reg=lm(y[,i][ok]~x[,j][ok]),col=col.line)
			}
			
		if(!j==ncol(x) | !i==ncol(y)) {readline("Press return for next plot ")}		

		} #end loop thru x variables
	} #end loop thru y variables
	
par(oldpar)
}
sum.stats <-
function(x,var='',by='',margin='column',...){

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

variable<-colnames(y)
sample<-rownames(y)

#statistical functions
nobs<-function(x) length(x)
cv<-function(x,na.rm) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100 
zeros<-function(x,na.rm) sum(x==0,na.rm=TRUE)
pct.zeros<-function(x,na.rm) sum(x==0,na.rm=TRUE)/length(x)*100
nobs.missing<-function(x,na.rm) sum(is.na(x))
pct.missing<-function(x,na.rm) sum(is.na(x))/length(x)*100 
se<-function(x,na.rm) sd(x,na.rm=TRUE)/sqrt(length(x)-sum(is.na(x)))
se.ratio<-function(x,na.rm) se(x)/mean(x,na.rm=TRUE)*100
richness<-function(x,na.rm) nobs(x)-zeros(x)-nobs.missing(x)
sh.diversity<-function(x,na.rm) -sum(((x)/sum(x,na.rm=TRUE))*log(((x)/sum(x,na.rm=TRUE))),na.rm=TRUE)
sh.evenness<-function(x,na.rm) sh.diversity(x)/log(richness(x))
si.diversity<-function(x,na.rm){
	if(richness(x)==0) 0
	else 1-sum(((x)/sum(x,na.rm=TRUE))*((x)/sum(x,na.rm=TRUE)),na.rm=TRUE)
	}
si.evenness<<-function(x,na.rm) si.diversity(x)/(1-(1/richness(x)))

if(by==''){ #summary table w/o groups
	if(margin=='column'){ #summary table by columns
		z1<-data.frame(apply(y,2,function(x){ #calculate stats
			z1<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
		    mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sum(x,na.rm=TRUE),
		    sd(x,na.rm=TRUE),cv(x),zeros(x),pct.zeros(x),nobs.missing(x),
			pct.missing(x),se(x),se.ratio(x),richness(x),sh.diversity(x),
			sh.evenness(x),si.diversity(x),si.evenness(x))
		    names(z1)<-c('nobs','min','max','mean',
		    'median','sum','sd','cv','zeros','pct.zeros',
		    'nobs.missing','pct.missing','se','se.ratio',
		    'richness','sh.diversity','sh.evenness',
		    'si.diversity','si.evenness') #create col names
			z1<-round(z1,3) #round elements to 3 decimal places
			}))
		z2<-data.frame(t(apply(z1,1,function(x){ #calculate stats on stats
			z2<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
			mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sd(x,na.rm=TRUE),cv(x))
		    names(z2)<-c('nobs','min','max','mean',
		    'median','sd','cv') #create row names
			z2<-round(z2,3) #round elements to 3 decimal places
			})))
		z<-list(z1,z2) #create list with col stats and sum stats
		names(z)<-c('Column.Summary','Table.Summary')
		} #end summary table by columns

	else{ #summary table by rows
		z1<-data.frame(t(apply(y,1,function(x){ #calculate stats
			z1<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
		    mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sum(x,na.rm=TRUE),
		    sd(x,na.rm=TRUE),cv(x),zeros(x),pct.zeros(x),nobs.missing(x),
			pct.missing(x),se(x),se.ratio(x),richness(x),sh.diversity(x),
			sh.evenness(x),si.diversity(x),si.evenness(x))
		    names(z1)<-c('nobs','min','max','mean',
		    'median','sum','sd','cv','zeros','pct.zeros',
		    'nobs.missing','pct.missing','se','se.ratio',
		    'richness','sh.diversity','sh.evenness',
		    'si.diversity','si.evenness') #create col names
			z1<-round(z1,3) #round elements to 3 decimal places
			})))
		z2<-data.frame(apply(z1,2,function(x){ #calculate stats on stats
			z2<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
			mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sd(x,na.rm=TRUE),cv(x))
		    names(z2)<-c('nobs','min','max','mean',
		    'median','sd','cv') #create row names
			z2<-round(z2,3) #round elements to 3 decimal places
			}))
		z<-list(z1,z2) #create list with row stats and sum stats
		names(z)<-c('Row.Summary','Table.Summary')
		} #end summary table by rows
	} #end summary table w/o groups

else{ #summary table w/ groups
#	write('',file=paste(outfile,'.csv',sep='')) #empty outfile if it exists
	fns<-c('nobs','min','max','mean',
	    'median','sum','sd','cv','zeros','pct.zeros',
	    'nobs.missing','pct.missing','se','se.ratio',
	    'richness','sh.diversity','sh.evenness',
	    'si.diversity','si.evenness') #create names vector
	n<-by.names(x,by) #create by variable
	for(i in 1:length(fns)){ #loop thru by groups
		cat(t<-paste(strtrim(paste('--',fns[i],paste(rep('-',80),collapse='')),80),'\n')) #create line break
		q<-list(n[,2]) #create a list of group names
		names(q)<-names(n)[2] #assign by name to q
		z1<-aggregate(y,q,fns[i],na.rm=TRUE) #calculate stats
		zz1<-round(z1[,2:ncol(z1)],3) #round stats to 3 decimal places
		g<-z1[,1,,drop=FALSE] #select group variable
		z1<-cbind(g,zz1) #bind group variable with selected variables
		z2<-data.frame(t(apply(z1[,-1],1,function(x){ #calculate stats on stats
			z2<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
			mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sd(x,na.rm=TRUE),cv(x))
		    names(z2)<-c('nobs','min','max','mean',
		    'median','sd','cv') #create row names
			z2<-round(z2,3) #round elements to 3 decimal places
			})))
		z<-cbind(z1,z2) #bind col stats with summary stats
		print(z) #print to console
		} #end loop thru groups
	} #end summary table w/ groups
return(z)
}
tau <-
function(y,prior){

z<-matrix(0,nrow(y),ncol(y)) #create blank matrix
N<-sum(y)
ccr<-sum(diag(y))
n<-apply(y,1,sum)

num<-ccr-sum(prior*n)
den<-N-sum(prior*n)
tau<-num/den
tau[tau<0]<-0
return(tau)		
}
t.plot <-
function(df=1,ncp=0,xlim=c(-3,3)){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))
n.plots<-length(df)
k<-0
for(i in df){
  k<-k+1
	curve(pt(x,df=i,ncp=ncp),xlim[1],xlim[2],
		xlab='Z',ylab='Cumulative probability (quantile)',main='Cumulative Probability')
		mtext(paste('(df=',i,')',sep=''),side=3,col='red')
	curve(dt(x,df=i,ncp=ncp),xlim[1],xlim[2],
		xlab='Z',ylab='Probability density',main='Probability Density')
		mtext(paste('(df=',i,')',sep=''),side=3,col='red')
	curve(qt(x,df=i,ncp=ncp),0,1,
		xlab='Cumulative probability (quantile)',ylab='Z',main='Quantiles')
		mtext(paste('(df=',i,')',sep=''),side=3,col='red')
	y<-rt(1000,df=i,ncp=ncp)
	hist(y,xlab='Z',ylab='Frequency',main='Random Observations')
		mtext(paste('(df=',i,')',sep=''),side=3,col='red')
	if(k<n.plots) readline('press return for next plot')
	}
}
uv.outliers <-
function(x,id,var,by=NULL,outfile=NULL,sd.limit=3,digits=3){

#sd outliers w/o groups
if(is.null(by)){ 
	y1<-as.data.frame(subset(x,select=eval(parse(text=id)))) #select plot.id variables
	y2<-as.data.frame(subset(x,select=eval(parse(text=var)))) #select variables to standardize
	t1<-scale(y2) #calculate sd's
	t2<-abs(t1)>=sd.limit
	row.vector<-apply(t2,1,FUN='any',na.rm=TRUE)#select rows with extremes
	col.vector<-apply(t2,2,FUN='any',na.rm=TRUE)#select cols with extremes
    if(sum(row.vector)>0){	
		t3<-t1[row.vector,col.vector,drop=FALSE]
		t3[abs(t3)<sd.limit]<-NA
		t3<-round(t3,digits)
		t4<-as.data.frame(y1[row.vector,,drop=FALSE])
		z<-cbind(t4,t3)
		if(!is.null(outfile)){ #write table to outfile
			write.table(z,file=paste(outfile,'.csv',sep=''),row.names=FALSE,quote=FALSE,sep=',')
			} #end save outfile
		}
	else stop('No outliers exist')
	} #end sd outliers w/o groups

#sd outliers w/ groups
else{
	if(!is.null(outfile)) write('',file=outfile) #empty outfile if it exists
	n<-by.names(x,by) #create by variable
	y<-cbind(n,x)
	m<-levels(n[,2]) #create object with group levels
	z<-vector('list',length(m))	#create list object for output
	names(z)<-m #assign names to list components
	for(i in 1:length(m)){ #loop thru by groups
		t0<-y[n[,2]==m[i],,drop=FALSE] #select records within group
		y1<-as.data.frame(subset(t0,select=eval(parse(text=id)))) #select plot.id variables
		y2<-as.data.frame(subset(t0,select=eval(parse(text=var)))) #select variables to standardize
		t1<-scale(y2) #calculate sd's
		t2<-abs(t1)>=sd.limit
		row.vector<-apply(t2,1,FUN='any',na.rm=TRUE)#select rows with extremes
		col.vector<-apply(t2,2,FUN='any',na.rm=TRUE)#select cols with extremes
	    if(sum(row.vector)>0){	
			t3<-t1[row.vector,col.vector,drop=FALSE]
			t3[abs(t3)<sd.limit]<-NA
			t3<-round(t3,digits)
			t4<-as.data.frame(y1[row.vector,,drop=FALSE])
			z[[i]]<-cbind(t4,t3)
			if(!is.null(outfile)){ #write table to outfile
				write(m[i],file=paste(outfile,'.csv',sep=''),append=TRUE)
				write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,append=TRUE,sep=',')
				} #end save outfile
			}
		else z[[i]]<-NULL
		} #end loop thru groups
	} #end sd outliers w/groups
	
return(z)
}
uv.plots <-
function(x,var=NULL,col.fill='blue',col.point='black',
	col.line='red',...){

oldpar<-par(no.readonly=TRUE)

if(!is.null(var)){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)} #graphics settings

#layout plot
layout(matrix(c(1,1,2,3,4,5),
	nrow=3,ncol=2,byrow=TRUE),
	heights=c(.1,.45,.45),widths=c(.5,.5))
par(mar=c(1,1,1,1))

#loop thru variables
for(i in 1:ncol(y)){ 

	#plot title
	plot.new()
	text(.5,.5,labels=names(y[i]),
		font=2,cex=3,xpd=TRUE)
	
	#histogram (upper left panel)
	par(mar=c(5,5,4,2))
	hist(y[[i]],prob=TRUE,col=col.fill,
		xaxs='i',yaxs='i',xlab=names(y[i]),
		main='Histogram',...)
	lines(density(y[[i]]))
			
	#box-and-whisker plot (upper right panel)
	par(mar=c(5,5,4,2))
	boxplot(y[i],col=col.fill,ylab=names(y[i]),
		main='Box-and-Whisker Plot',...)

	#empirical distribution function plot (lower left panel)
	par(mar=c(5,5,4,2))
	plot(sort(y[[i]]),type='o',col=col.point,yaxs='i',xaxs='i',
		xlab='Plot Number',ylab=names(y[i]),
		main='Empirical Distribution',...)

	#normal quantile-quantile plot (lower right panel)
	par(mar=c(5,5,4,2))
	qqnorm(y[,i],datax=TRUE,col=col.point,
		main='Normal Q-Q Plot',...)
	y.IQR<-IQR(y[,i],na.rm=TRUE)
	if(y.IQR>0)	qqline(y[,i],datax=TRUE,col=col.line,...)

	par(mar=c(1,1,1,1))		
	if(!i==ncol(y)) {readline("Press return for next plot ")}		
		
	} #end loop thru variables
	
par(oldpar)
}
weibull.plot <-
function(shape=1,mri=10){
old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
par(mfrow=c(2,2))

#define time since event function
tweibull<-function(x,shape=shape,scale=scale){
  1-pweibull(x,shape=shape,scale=scale)
}

hweibull<-function(x,shape,mri){
  scale<-round(mri/gamma(1+(1/shape)),0)
  z<-(shape*(x^(shape-1)))/(scale^shape)
  return(z)  
}

#define scale as a function of mri and shape
scale<-round(mri/gamma(1+(1/shape)),0)

#create plots
curve(pweibull(x,shape,scale),from=0,to=2.5*scale,ylim=c(0,1),
	xlab='Time (between events)',
  ylab='Cumulative probability (quantile)',
  main='Cumulative Probability')
	mtext(paste('(shape=',shape,', mri=',mri,')',sep=''),side=3,col='red')
curve(dweibull(x,shape,scale),from=0,to=2.5*scale,
	xlab='Time (between events)',
  ylab='Probability density',
  main='Probability Density')
	mtext(paste('(shape=',shape,', mri=',mri,')',sep=''),side=3,col='red')
curve(tweibull(x,shape,scale),from=0,to=2.5*scale,ylim=c(0,1),
	xlab='Time (since last event)',
  ylab='Prob (not having an event up to time t)',
  main='Time Since Event')
	mtext(paste('(shape=',shape,', mri=',mri,')',sep=''),side=3,col='red')
curve(hweibull(x,shape,mri),from=0,to=5*scale,
  xlab='Time(between events)',
  ylab='Hazard',
  main='Hazard Function')
	mtext(paste('(shape=',shape,', mri=',mri,')',sep=''),side=3,col='red')
}
