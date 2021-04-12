sampleSelCoefs <- function(L=5000,dist='point',params=1){
    ## recover()
    if(dist == 'point'){
        if(is.matrix(params)){
            weights <- params[1,]/sum(params[1,])
            gammas <- sample(
                params[2,],
                size=L,
                replace=TRUE,
                prob=weights
            )
        } else {
            gammas <- rep(params,L)
        }
    } else if(dist=='exp'){
        gammas <- rexp(L,1/params)
    }
    return(gammas)
}
sampleEffects <- function(wn,gammas,polarize=TRUE){
    L <- length(gammas)
    effects <- rnorm(L,0,sqrt(wn*gammas))
    if(polarize){
        ## under stabilizing selection the sign of the allele's effect does
        ## not influence its dynamics, so we can just orient them all to be
        ## positive if we want without changing anything
        return(abs(effects))
    } else
        return(effects)
}
freqSpec <- function(x,gamma){
    if(gamma==0){
        (x*(1-x))^-1
    } else {
        exp(-gamma*x*(1-x))/(x*(1-x))
    }
}
sampleFrequencies <- function(gammas,pop.size=10000,maf.cutoff=0.01){
    ## currently this function assumes for simplicity that sample freq = pop freq
    ## note that this function is going to be really slow when there are a large
    ## number of sites and every site has a different selection coefficient, because
    ## there's a bunch of numerical integrals to compute.
    ## probably can be made faster, but not clear if necessary...

    mac.cutoff <- 2*pop.size*maf.cutoff
    my.macs <- mac.cutoff:(2*pop.size-mac.cutoff)
    gamma.counts <- table(gammas)
    unique.gammas <- sort(unique(gammas))
    num <- numeric()
    pop.freqs <- numeric(length(gammas))
    pb <- txtProgressBar(min=0,max=length(unique.gammas)*length(my.macs),style=3)
    j <- 1
    for( k in seq_along(gamma.counts)){


        denom <- integrate(
            f=function(X) freqSpec(X,unique.gammas[k]),
            lower=maf.cutoff-1/(4*pop.size),
            upper=1-maf.cutoff+1/(4*pop.size)
        )$value
        for(i in seq_along(my.macs)){
            num[i] <- integrate(
                f=function(X) freqSpec(X,unique.gammas[k]),
                lower=(my.macs[i]-1/2)/(2*pop.size),
                upper=(my.macs[i]+1/2)/(2*pop.size)
            )$value
            j <- j+1
            setTxtProgressBar(pb,j)
        }
        pmf <- num/denom
        pop.freqs[gammas==unique.gammas[k]] <- sample(
            my.macs/(2*pop.size),
            size=gamma.counts[k],
            replace=TRUE,
            prob=pmf
        )
    }
    return(pop.freqs)
}



wn <- 10 # w^2/n...for the purposes of simulated data, this parameter is basically
## just an arbitrary number that sets the scale of the simulated effect sizes.
## In an estimation setting, it will be scaled according to however the effect
## sizes are scaled, and derives it's meaning from that scale


## simulation of an architecture with 99% of loci having a scaled selection coefficient
## of gamma=3, and 1% having gamma=50
point.gammas <- sampleSelCoefs(params=matrix(c(99,1,5,50),nrow=2))
point.effects <- sampleEffects(10,point.gammas)
point.freqs <- sampleFrequencies(point.gammas,maf.cutoff=0.01)
point.vars <- point.effects^2*point.freqs*(1-point.freqs)
point.thr <- quantile(point.vars,probs=0.7)
above.point.thr <- point.vars>point.thr

pdf('../figures/pointSelCoefs.pdf',width=8,height=6)
plot(
    x=point.freqs[above.point.thr],
    y=point.effects[above.point.thr],
    bty='n',
    pch=20,
    xlim=c(0,1),
    xlab='Frequency of Trait Increasing Allele',
    ylab='Effect Size',
    ylim=c(0,max(my.effects)*1.05)
)
dev.off()



## simulation of architecture with an exponential distribution of scaled selection coefficients sizes
my.gammas <- sampleSelCoefs(L=5000,dist='exp',params=1)
my.effects <- sampleEffects(10,my.gammas)
my.freqs <- sampleFrequencies(my.gammas,maf.cutoff=0.01)
vars <- my.effects^2*my.freqs*(1-my.freqs)
my.thr <- quantile(vars,probs=0.7)
above.thr <- vars>my.thr

pdf('../figures/expSelCoefs.pdf',width=8,height=6)
plot(
    x=my.freqs[above.thr],
    y=my.effects[above.thr],
    bty='n',
    pch=20,
    xlim=c(0,1),
    xlab='Frequency of Trait Increasing Allele',
    ylab='Effect Size',
    ylim=c(0,max(my.effects)*1.05)
)
dev.off()

## assume everything is neutral and effects come from a single Gaussian
L <- 5000
neut.gammas <- rep(0,L)
neut.effects <- rnorm(L,0,wn)
neut.freqs <- sampleFrequencies(neut.gammas)
neut.vars <- neut.effects^2*neut.freqs*(1-neut.freqs)
neut.thr <- quantile(neut.vars,probs=0.7)
above.neut.thr <- neut.vars > neut.thr

pdf('../figures/neutSelCoefs.pdf',width=8,height=6)
plot(
    x=neut.freqs[above.neut.thr],
    y=neut.effects[above.neut.thr],
    bty='n',
    pch=20,
    xlim=c(0,1),
    xlab='Frequency of Trait Increasing Allele',
    ylab='Effect Size',
    ylim=c(0,max(neut.effects)*1.05)
)
dev.off()
