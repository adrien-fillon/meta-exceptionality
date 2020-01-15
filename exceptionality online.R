
###############################
## Install and load packages ##
###############################

install.packages(c("plyr", "dplyr", "dplyr", "rstudioapi", "foreign", "metafor","pequod", "csv", "ggplot2", "xlsx",  "psych", 
                   "compute.es", "devtools", "broom", "MBESS", "formatR",
                   "MAd", "data.table", "weightr", "powerAnalysis", "metaforest", "forcats", "rmarkdown","esc"))
install.packages("stringi", repos="http://cran.rstudio.com/", dependencies=TRUE)

Sys.setenv(LANG = "en")
library(plyr)
library(dplyr)
library(foreign)
library(rstudioapi)
library(metafor)
library(pequod)
library(csv)
library(ggplot2)
library(xlsx)
library(psych)
library(compute.es)
library(devtools)
library(MBESS)
library(formatR)
library(data.table)
library(weightr)
library(powerAnalysis)
library(metaforest)
library(forcats)
library(Hmisc)
library(esc)
library(broom)
library (MAd)
library(rmarkdown)

# setting formatting options
options(scipen=999, digits =3)


# Setting the working directory to the same folder as the .R file is in
# this.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(this.dir)


# Code obtained on 02/10/2016 from:
# http://journals.sagepub.com/doi/suppl/10.1177/1745691616662243 

init.value <- c(0.5, 0.2, 0.10)
alpha <- 0.05		# Size of two-sided test

estimate.onestep.selection.heterogeneous <- function(z.obs,n1,n2,alpha,theta.init){
  p.value <- 1 - pnorm(z.obs)
  sel <- z.obs > 0 & p.value < alpha
  
  if(sum(sel)==0 | sum(sel)==length(sel)){
    if(sum(sel)==length(sel)){ w <- 1/(length(sel)+2) }		# Not identified; use Wilson-like
    if(sum(sel)==0){ w <- 1 - 1/(length(sel)+2) }			# estimator.
    e <- estimate.known.onestep.selection.heterogeneous(z.obs,n1,n2,w,alpha, theta.init[1:2])
    tmp <- e[[2]]
    tmp <- cbind(tmp, c(NA,NA))
    tmp <- rbind(tmp, c(NA,NA,NA))
    return(list(est=c(e[[1]],w), est.var=tmp, ll=e[[3]]))
  }
  
  if(sum(sel) > 0 & sum(sel)<length(sel)){
    theta.init <- c(theta.init[1], log(theta.init[2:3]))
    tmpf <- function(theta,z.obs,n1,n2,alpha){ onestep.heterogeneous.nll(c(theta[1],exp(theta[2:3])),z.obs,n1,n2,alpha) }
    tmpg <- function(theta,z.obs,n1,n2,alpha){ onestep.heterogeneous.nll(c(theta[1],theta[2:3]),z.obs,n1,n2,alpha) }
    tmpo <- optim(theta.init, tmpf, z.obs=z.obs, n1=n1, n2=n2, alpha=alpha)
    theta.hat <- c(tmpo$par[1], exp(tmpo$par[2:3]))
    tmpv <- matrix(NA, 3, 3)
    suppressWarnings(try( tmpv <- solve(optimHess(theta.hat, tmpg, z.obs=z.obs, n1=n1, n2=n2, alpha=alpha)), silent=TRUE))
    return(list(est=theta.hat, est.var=tmpv, ll=-tmpo$value))
  }	
}

# Loading the master dataset and labeling it "master"
dat = read.csv("past-behavior-meta-coding.csv", skip = 3)
# head(dat)
# str(dat, list.len=ncol(dat))

# let's just keep the data that we need for effect-size calculations
dataset <- dat[, c("Article", 
                   "Article.number",
                   "Study..",
                   "Sample..",
                   "Published",
                   "DV..",
                   "N.sample.size..gross..",
                   "N.sample.size..post.attrition.",
                   "Exceptional.Cell.N","Routine.Cell.N",
                   "IV1.manipulated.",
                   "Design.type",
                   "DV.group",
                   "DV.percentage.",
                   "DV.count.",
                   "DV.scale.", 
                   "Exceptional.Cell.M", 
                   "Routine.Cell.M", 
                   "Exceptional.Cell.SD", 
                   "Routine.Cell.SD", 
                   "F.ANOVA.F", 
                   "F.ANOVA.DF..x..y.", 
                   "t.from.t.test", 
                   "t.test.df..x.", 
                   "p.values", 
                   "Cnt.Cell1.option1", 
                   "Cnt.Cell1.option2", 
                   "p.values.1", 
                   "Calculated.Chisq", 
                   "Reported.Chisq",  
                   "Calculated.Chisq.1", 
                   "X..Cell1", 
                   "X..Cell2", 
                   "Does.effect.need.to.be.reversed.",
                   "Best.effect.estimate", 
                   "Action.inaction", 
                   "Status.quo", 
                   "Controllable.vs..uncontrollable.Events", 
                   "Exceptional.vs..normal.outcome",
                   "Outcome.severity",
                   "Regret.vs..counterfactuals",
                   "Routine.strength",
                   "Regret.vs..Counterfactuals..using.DV.groups.")]

# let's combine the chisquare columns from proportions and counts for clarity
dataset$chisq <- rowSums(dataset[,c("Calculated.Chisq", "Calculated.Chisq.1")], na.rm = TRUE)
dataset$chisq[dataset$chisq == 0] <- NA

# set up an empty frame to hold the effect size
dat2 <- data.frame(matrix(NA, nrow=nrow(dataset), ncol=2)) 
names(dat2) <- c("yi","vi")

# for (i in 1:nrow(dataset)){
for (i in 1:nrow(dataset)) {  
  # if we only have t statistic and N per each condition
  if (!is.na(dataset$Exceptional.Cell.N[i]) & 
      !is.na(dataset$Routine.Cell.N[i]) &
      !is.na(dataset$t.from.t.test[i])) {
    tmp2 <- NA
    tmp2 <- esc_t(t = dataset$t.from.t.test[i], 
                  grp1n = dataset$Exceptional.Cell.N[i], 
                  grp2n = dataset$Routine.Cell.N[i])
    if (dataset$Does.effect.need.to.be.reversed.[i]== "Yes") {
      tmp2$es <- tmp2$es *(-1)
    }
    dat2$yi[i] <- tmp2$es
    dat2$vi[i] <- tmp2$var
  } else 
    # if we only have t statistic and overall N
    if (!is.na(dataset$N.sample.size..post.attrition.[i]) & 
        !is.na(dataset$t.from.t.test[i])) {
      tmp3 <- NA
      tmp3 <- esc_t(t = dat$t.from.t.test[i], 
                    totaln = dat$N.sample.size..post.attrition.[i])
      if (dataset$Does.effect.need.to.be.reversed.[i]== "Yes") {
        tmp2$es <- tmp2$es *(-1)
      }
      dat2$yi[i] <- tmp3$es
      dat2$vi[i] <- tmp3$var
    } else 
      # if we only have F statistic and DF
      if (!is.na(dataset$N.sample.size..post.attrition.[i]) & 
          !is.na(dataset$t.from.t.test[i])) {
        tmp4 <- NA
        tmp4 <- esc_t(t = dat$t.from.t.test[i], 
                      totaln = dat$N.sample.size..post.attrition.[i])
        if (dataset$Does.effect.need.to.be.reversed.[i]== "Yes") {
          tmp2$es <- tmp2$es *(-1)
        }
        dat2$yi[i] <- tmp4$es
        dat2$vi[i] <- tmp4$var
      } else 
        # if we only have Cohen's d and Ns
        if (!is.na(dataset$Best.effect.estimate[i]) & 
            !is.na(dataset$Exceptional.Cell.N[i]) & 
            !is.na(dataset$Routine.Cell.N[i])) {
          dat2$yi[i] <- dataset$Best.effect.estimate[i]
          dat2$vi[i] <- 1/dataset$Exceptional.Cell.N[i]+ 
            1/dataset$Routine.Cell.N[i] + 
            dataset$Best.effect.estimate[i]^2/(2*(dataset$Exceptional.Cell.N[i]+dataset$Routine.Cell.N[i]))
        } else 
          # if we only have Cohen's d and overall
          if (!is.na(dataset$Best.effect.estimate[i]) & 
              !is.na(dataset$N.sample.size..post.attrition.[i])) {
            dat2$yi[i] <- dataset$Best.effect.estimate[i] 
            dat2$vi[i] <- 1/(dataset$N.sample.size..post.attrition.[i]/2) + 
              1/(dataset$N.sample.size..post.attrition.[i]/2) + 
              dataset$Best.effect.estimate[i]^2 / 
              (2*((dataset$N.sample.size..post.attrition.[i]/2)+(dataset$N.sample.size..post.attrition.[i]/2)))
          } else 
            # if we only have chisquare and N
            if (!is.na(dataset$chisq[i]) & 
                !is.na(dataset$N.sample.size..post.attrition.[i])) {
              tmp5 <- NA
              tmp5 <- chies(dataset$chisq[i], 
                            dataset$N.sample.size..post.attrition.[i], 
                            level = 95, dig = 2, verbose = TRUE, id=NULL, data=NULL)
              if (dataset$Does.effect.need.to.be.reversed.[i]== "Yes") {
                tmp5$d <- tmp5$d *(-1)
              }
              dat2$yi[i] <- tmp5$d
              dat2$vi[i] <- tmp5$var.d
            } else 
              # if we have M, SD, N for two conditions
              if (!is.na(dataset$Exceptional.Cell.N[i]) & 
                  !is.na(dataset$Exceptional.Cell.M[i]) & 
                  !is.na(dataset$Exceptional.Cell.SD[i]) & 
                  !is.na(dataset$Routine.Cell.N[i]) &
                  !is.na(dataset$Routine.Cell.M[i]) &
                  !is.na(dataset$Routine.Cell.SD[i])) {
                tmp1 <- NA
                tmp1 <- escalc(measure="SMD", 
                               m1i=dataset$Exceptional.Cell.M[i], 
                               sd1i=dataset$Exceptional.Cell.SD[i], 
                               n1i=dataset$Exceptional.Cell.N[i],
                               m2i=dataset$Routine.Cell.M[i], 
                               sd2i=dataset$Routine.Cell.SD[i], 
                               n2i=dataset$Routine.Cell.N[i], append = FALSE)
                if (dataset$Does.effect.need.to.be.reversed.[i]== "Yes") {
                  tmp1$yi <- tmp1$yi *(-1)
                }
                dat2$yi[i] <- tmp1$yi 
                dat2$vi[i] <- tmp1$vi
              } 
  
}

dataset <- cbind(dataset,dat2)
dataset$d <- dataset$yi
dataset$dvar <- dataset$vi

# convert Cohen's d to Hedge's g
dataset$yi <- d_to_g(dataset$d, dataset$dvar,
                     ifelse(!is.na(dataset$Exceptional.Cell.N),dataset$Exceptional.Cell.N, dataset$N.sample.size..post.attrition./2),
                     ifelse(!is.na(dataset$Routine.Cell.N),dataset$Routine.Cell.N, dataset$N.sample.size..post.attrition./2))[,1]
dataset$vi <- d_to_g(dataset$d, dataset$dvar,
                     ifelse(!is.na(dataset$Exceptional.Cell.N),dataset$Exceptional.Cell.N, dataset$N.sample.size..post.attrition./2),
                     ifelse(!is.na(dataset$Routine.Cell.N),dataset$Routine.Cell.N, dataset$N.sample.size..post.attrition./2))[,2]

## 
# getting ready for moderator analyses
# first, we set up the variables

dataset$Action.inaction[dataset$Action.inaction == 99] <- NA

dataset$Controllable.vs..uncontrollable.Events[dataset$Controllable.vs..uncontrollable.Events == 99] <- NA

### more moderators added by Lucas 8/24/2017

dataset$Exceptional.vs..normal.outcome[dataset$Exceptional.vs..normal.outcome == 99] <- NA

dataset$Status.quo[dataset$Status.quo == 99] <- NA

dataset$Outcome.severity[dataset$Outcome.severity == 99] <- NA



# we'll skip this moderator since we run the analyses separated for DV type
# dataset$Regret.vs..counterfactuals[dataset$Regret.vs..counterfactuals == 3] <- NA

dataset$Routine.strength.recoded <- 
  ifelse(dataset$Routine.strength == "certainty (always/never/each day)", 3, 
         ifelse(dataset$Routine.strength == "high frequency (often/regularly/most times/habit", 2, 
                ifelse(dataset$Routine.strength == "medium frequency (several times/occasionally", 1, NA)))


# str(dataset)
# describe(dataset)

##Results for experimental studies ##

### Exception/routine experimental studies meta-analysis summary ###

##Conducting the meta-analyis##

# collapsing first for manipulated studies
dataset1 <- subset(dataset, IV1.manipulated.=="Manipulated")
dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
dataset1$articlestudy <- factor(dataset1$articlestudy)

collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                 method="BHHR", mod=NULL, data=dataset1)
# str(collapsed)

###
# if we collapse the effects within each study
byarticlestudy <- 
  cbind(
    aggregate(
      yi ~ articlestudy, 
      data=dataset1, 
      FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
    aggregate(
      N.sample.size..post.attrition. ~ articlestudy, 
      data=dataset1, FUN = max)[2]
  )
names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
totalexp <- merge(byarticlestudy, collapsed,by="articlestudy")

allexpcollapsed <- rma(es, var, method = "REML", data = totalexp, slab = articlestudy)
allexpcollapsed


confint(allexpcollapsed, digits=3)

# code from http://www.metafor-project.org/doku.php/analyses:konstantopoulos2011
# same two level analysis with the rma.mv
allexpcollapsedml2 <- rma.mv(yi, vi, random = ~ 1 | articlestudy, data=dataset1)
allexpcollapsedml2
print(allexpcollapsedml2, digits=3)

# code from http://www.metafor-project.org/doku.php/analyses:konstantopoulos2011
# three level, taking into account same article/authors
allexpcollapsedml3 <- rma.mv(yi, vi, random = ~ 1 | Article/articlestudy, data=dataset1)
allexpcollapsedml3
# print(allexpcollapsedml3, digits=3)

par(mfrow=c(2,1))
profile(allexpcollapsedml3, sigma2=1)
profile(allexpcollapsedml3, sigma2=2)

# code from http://www.metafor-project.org/doku.php/analyses:konstantopoulos2011
# Estimated ICC Intraclass Correlation of the True Effects
# The three-level model used in the meta-analysis above allows for the underlying true effects within articles to be correlated.
ICC = round(allexpcollapsedml3$sigma2[1] / sum(allexpcollapsedml3$sigma2), 3)
ICC

allexpcollapsedml4 <- rma.mv(yi, vi, random = ~ factor(articlestudy) | Article, data=dataset1)
allexpcollapsedml4

###Summary###

### Forest plot with Bells and whistles ###

# Create a table of information from the res model
allexpcollapsed.restable<-cbind.data.frame(allexpcollapsed$b, allexpcollapsed$se, 
                                           allexpcollapsed$pval, allexpcollapsed$ci.lb, allexpcollapsed$ci.ub,  allexpcollapsed$k)
colnames(allexpcollapsed.restable)<-c("Hedge's g","SE", "p", "CI Lower", "CI Upper", "k")
row.names(allexpcollapsed.restable)<-paste(X,"-",Y)

# Create a table of heterogeneity information
allexpcollapsed.hettable<-cbind.data.frame(allexpcollapsed$QE, 
                                           round(allexpcollapsed$QEp, digits = 3), allexpcollapsed$I2)
colnames(allexpcollapsed.hettable)<-c("Q", "p", "I2")
row.names(allexpcollapsed.hettable)<-paste(X,"-",Y)

# Create a table of heterogeneity information
allexpcollapsed.hettable<-cbind.data.frame(allexpcollapsed$QE, 
                                           ifelse(allexpcollapsed$QEp < 0.001, "< .001", round(allexpcollapsed$QEp,3)), allexpcollapsed$I2)
colnames(allexpcollapsed.hettable)<-c("Q", "p", "I2")

###Forest plot###

r fig.width=11.5, fig.height=12          
forest(allexpcollapsed, alim=c(-2,4), xlim=c(-4,5), ilab=totalexp$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Past behavior normality and regret (comparisons)", cex=.8)
#add headlines to the forest
op <- par(font=4)
text(-4, allexpcollapsed$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
text( 5, allexpcollapsed$k+2, "Observed [95% CI]",  pos=2, cex=.8)
text(-2, allexpcollapsed$k+2, "Sample size",cex=.8)
par(op)

###Funnel plot###

r fig.width=14, fig.height=6
funnel(allexpcollapsed)

# http://www.metafor-project.org/doku.php/plots:contour_enhanced_funnel_plot
# Note that the funnel is centered not at the model estimate (as is usually done when drawing funnel plots), but at 0 (i.e., at the value under the null hypothesis of no effect). Various levels of statistical significance of the points/studies are indicated by the shaded regions. In particular, the unshaded (i.e., white) region in the middle corresponds to p-values greater than .10, the gray-shaded region corresponds to p-values between .10 and .05, the dark gray-shaded region corresponds to p-values between .05 and .01, and the region outside of the funnel corresponds to p-values below .01. Funnel plots drawn in this way are more useful for detecting publication bias due to the suppression of non-significant findings. See Peters et al. (2008) for more details. Note that, based on Sterne and Egger (2001), the vertical axis represents the standard error (as compared to Peters et al., 2008, who use the inverse of the standard error on the vertical axis).
funnel(allexpcollapsed, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)

### Publication Bias ###


#### Basic tests ####


Please read Schwarzer's book chapter 5 pages 107-138 for for explanation of trimfill and failsafe. 
(Also, can check Borenstein et al. (2009)'s  book, pages 284-287 )

** Trim and Fill**
  
  # run trim and fill
  taf <- trimfill(allexpcollapsed)
  taf
  
  # let's see this in the funnels
  # funnel after trim and fill
  funnel(taf)  
  
  ** Rank test**
    
    
    # to make the assymetry test more statistical rather than visual, let's do a ranktest
    ranktest(allexpcollapsed)
  
  
  
  ** Reg test**
    
    regtest(allexpcollapsed)
  
  ** Fail safe **
    
    # fail safe
    # https://rdrr.io/cran/metafor/man/fsn.html
    fsn(es, var, data = totalexp, type="Rosenthal", alpha=.05)
  
  #### Advanced tests ####
  
  ** PET PEESE **
    
    # From http://daniellakens.blogspot.com/2015/04/why-meta-analysis-of-90-precognition.html 
    # first let's convert from our meta to Daniel's variables
    
    #calculate variance and standard error
    totalexp$SE<-sqrt(totalexp$var)
  
  #PET
  PET<-lm(totalexp$es~totalexp$SE, weights = 1/totalexp$var)
  summary(PET)
  confint(PET)
  print(c(summary(PET)$coefficients[1], confint(PET)[1,1], confint(PET)[1,2]))
  
  #PEESE
  PEESE<-lm(totalexp$es~totalexp$var, weights = 1/totalexp$var)
  summary(PEESE)
  confint(PEESE)
  print(c(summary(PEESE)$coefficients[1], confint(PEESE)[1,1], confint(PEESE)[1,2]))
  
  ** puniform **
    
    if(!require(puniform)){install.packages('puniform')}
  library(puniform)
  
  puniform (yi=totalexp$es, vi=totalexp$var, side="right", method="P", plot = "FALSE")
  
  ** Three-parameter selection model **
    
    # read "Correcting for bias in psychology: A comparison of meta-analytic methods"
    # link: https://osf.io/preprints/psyarxiv/9h3nu
    # Conclusion: "Our results indicate that one method - the three-parameter selection model (Iyengar & Greenhouse, 1988; McShane, B?ckenholt, & Hansen, 2016) - generally performs better than trim-and-fill, p-curve, p-uniform, PET, PEESE, or PET-PEESE,"
    # McShane, B. B., B?ckenholt, U., & Hansen, K. T. (2016). Adjusting for publication bias
    # in meta-analysis: an evaluation of selection methods and some cautionary notes. Per-
    # spectives on Psychological Science, 11(5), 730-749.
    
    # Three-parameter selection model
    # weightr: https://www.rdocumentation.org/packages/weightr/versions/1.0.0/topics/weightfunct
    weightfunct(totalexp$es, totalexp$var, steps=c(.05/2,1))
  
  ** Henmi & Copas (2010) **
    
    # Henmi & Copas (2010)
    # http://www.metafor-project.org/doku.php/analyses:henmi2010
    hc(allexpcollapsed)
  
  ##Results for compared studies ##
  
  ### Exception/routine compared studies meta-analysis summary ###
  
  # collapsing second for compared studies
  dataset1 <- subset(dataset, IV1.manipulated.=="Measured/compared")
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  ###
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  totalcomp <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  allcompcollapsed <- rma(es, var, method = "REML", data = totalcomp, slab = articlestudy)
  allcompcollapsed
  
  ###Summary###
  
  This analysis is based on `r allcompcollapsed$k` compared studies that evaluated the impact of `r X` over `r Y`. 
  
  ### Forest plot with Bells and whistles ###
  
  # Create a table of information from the res model
  allcompcollapsed.restable<-cbind.data.frame(allcompcollapsed$b, allcompcollapsed$se, 
                                              allcompcollapsed$pval, allcompcollapsed$ci.lb, allcompcollapsed$ci.ub,  allcompcollapsed$k)
  colnames(allcompcollapsed.restable)<-c("Hedge's g","SE", "p", "CI Lower", "CI Upper", "k")
  row.names(allcompcollapsed.restable)<-paste(X,"-",Y)
  
  # Create a table of heterogeneity information
  allcompcollapsed.hettable<-cbind.data.frame(allcompcollapsed$QE, 
                                              round(allcompcollapsed$QEp, digits = 3), allcompcollapsed$I2)
  colnames(allcompcollapsed.hettable)<-c("Q", "p", "I2")
  row.names(allcompcollapsed.hettable)<-paste(X,"-",Y)
  
  # Create a table of heterogeneity information
  allcompcollapsed.hettable<-cbind.data.frame(allcompcollapsed$QE, 
                                              ifelse(allcompcollapsed$QEp < 0.001, "< .001", round(allcompcollapsed$QEp,3)), allcompcollapsed$I2)
  colnames(allcompcollapsed.hettable)<-c("Q", "p", "I2")
  
  
  ###Forest plot###
  
  r fig.width=14, fig.height=12         
  forest(allcompcollapsed, alim=c(-2,4), xlim=c(-4,5), ilab=totalcomp$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Past behavior normality and regret (comparisons)", cex=.8)
  #add headlines to the forest
  op <- par(font=4)
  text(-4, allcompcollapsed$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 5, allcompcollapsed$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, allcompcollapsed$k+2, "Sample size",cex=.8)
  par(op)
  
  ###Funnel plot###
  
  r fig.width=14, fig.height=6
  funnel(allcompcollapsed)
  
  
  # http://www.metafor-project.org/doku.php/plots:contour_enhanced_funnel_plot
  # Note that the funnel is centered not at the model estimate (as is usually done when drawing funnel plots), but at 0 (i.e., at the value under the null hypothesis of no effect). Various levels of statistical significance of the points/studies are indicated by the shaded regions. In particular, the unshaded (i.e., white) region in the middle corresponds to p-values greater than .10, the gray-shaded region corresponds to p-values between .10 and .05, the dark gray-shaded region corresponds to p-values between .05 and .01, and the region outside of the funnel corresponds to p-values below .01. Funnel plots drawn in this way are more useful for detecting publication bias due to the suppression of non-significant findings. See Peters et al. (2008) for more details. Note that, based on Sterne and Egger (2001), the vertical axis represents the standard error (as compared to Peters et al., 2008, who use the inverse of the standard error on the vertical axis).
  funnel(allcompcollapsed, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)
  
  ### Publication Bias ###
  
  #### Basic tests ####
  
  Please read Schwarzer's book chapter 5 pages 107-138 for for explanation of trimfill and failsafe. 
(Also, can check Borenstein et al. (2009)'s  book, pages 284-287 )

** Trim and Fill**
  # run trim and fill
  taf <- trimfill(allcompcollapsed)
  taf
  
  # let's see this in the funnels
  # funnel after trim and fill
  funnel(taf)  
  
  ** Rank test**
    
    
    # to make the assymetry test more statistical rather than visual, let's do a ranktest
    ranktest(allcompcollapsed)
  
  ** Reg test**
    
    regtest(allcompcollapsed)
  
  ** Fail safe **
    
    # fail safe
    # https://rdrr.io/cran/metafor/man/fsn.html
    fsn(es, var, data = totalcomp, type="Rosenthal", alpha=.05)
  
  
  #### Advanced tests ####
  
  ** PET PEESE **
    
    # From http://daniellakens.blogspot.com/2015/04/why-meta-analysis-of-90-precognition.html 
    # first let's convert from our meta to Daniel's variables
    
    #calculate variance and standard error
    totalcomp$SE<-sqrt(totalcomp$var)
  
  #PET
  PET<-lm(totalcomp$es~totalcomp$SE, weights = 1/totalcomp$var)
  summary(PET)
  confint(PET)
  print(c(summary(PET)$coefficients[1], confint(PET)[1,1], confint(PET)[1,2]))
  
  #PEESE
  PEESE<-lm(totalcomp$es~totalcomp$var, weights = 1/totalcomp$var)
  summary(PEESE)
  confint(PEESE)
  print(c(summary(PEESE)$coefficients[1], confint(PEESE)[1,1], confint(PEESE)[1,2]))
  
  
  ** puniform **
    
    if(!require(puniform)){install.packages('puniform')}
  library(puniform)
  
  puniform (yi=totalcomp$es, vi=totalcomp$var, side="right", method="P", plot = "FALSE")
  
  ** Three-parameter selection model **
    
    # read "Correcting for bias in psychology: A comparison of meta-analytic methods"
    # link: https://osf.io/preprints/psyarxiv/9h3nu
    # Conclusion: "Our results indicate that one method - the three-parameter selection model (Iyengar & Greenhouse, 1988; McShane, B?ckenholt, & Hansen, 2016) - generally performs better than trim-and-fill, p-curve, p-uniform, PET, PEESE, or PET-PEESE,"
    # McShane, B. B., B?ckenholt, U., & Hansen, K. T. (2016). Adjusting for publication bias
    # in meta-analysis: an evaluation of selection methods and some cautionary notes. Per-
    # spectives on Psychological Science, 11(5), 730-749.
    
    # Three-parameter selection model
    # weightr: https://www.rdocumentation.org/packages/weightr/versions/1.0.0/topics/weightfunct
    weightfunct(totalcomp$es, totalcomp$var, steps=c(.05/2,1))
  
  ** Henmi & Copas (2010) **
    
    # Henmi & Copas (2010)
    # http://www.metafor-project.org/doku.php/analyses:henmi2010
    hc(allcompcollapsed)
  
  
  ## Comparing the two effects: Experimental versus compared ##
  
  compareexptocomp <- data.frame(
    estimate = as.numeric(c(allexpcollapsed$b, allcompcollapsed$b)), 
    stderror = as.numeric(c(allexpcollapsed$se, allcompcollapsed$se)), 
    numstudies = as.numeric(c(allexpcollapsed$k, allcompcollapsed$k)),
    meta = as.character(c("Experimental","Compared")), 
    tau2 = as.numeric(c(c(allexpcollapsed$tau2, allcompcollapsed$tau2))))
  compareexptocomp
  
  comparetwosides <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=compareexptocomp, digits=3)
  comparetwosides
  
  <!-- MODERATORS -->
    
    # Moderator analyses #
    
    ## Results for experimental studies ##
    
    # we'll setup a copy dataaset to clean output a bit for the Metaforest limited functions
    dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Design.type=="Between")),]
  
  datasettmp <- dataset1
  datasettmp$Effect.size <- datasettmp$yi
  
  datasettmp$DV.group <- mapvalues(datasettmp$DV.group, from = c("1 - regret", "2 - counterfactuals", "3 - compensation", "4 - self blame", "5 - avoidability", "6 - other", "7 - offender punishment"), to = c("Regret", "Counterfactuals", "Compensation", "Self blame", "Avoidability", "Other", "Offender punishment"))
  datasettmp$Exceptional.vs..normal.outcome <- mapvalues(datasettmp$Exceptional.vs..normal.outcome, from = c("0", "1"), to = c("Routine", "Exceptional"))
  datasettmp$Outcome.severity <- mapvalues(datasettmp$Outcome.severity, from = c("0", "1"), to = c("Low severity", "High severity"))
  datasettmp$Action.inaction <- mapvalues(datasettmp$Action.inaction, from = c("0", "1"), to = c("Inaction", "Action"))
  datasettmp$Controllable.vs..uncontrollable.Events <- mapvalues(datasettmp$Controllable.vs..uncontrollable.Events, from = c("0", "1"), to = c("Uncontrollable", "Controllable"))
  datasettmp$Published <- mapvalues(datasettmp$Published, from = c("0", "1"), to = c("Unpublished", "Published"))
  datasettmp$Routine.strength.recoded <- mapvalues(datasettmp$Routine.strength.recoded, from = c("1", "2", "3"), to = c("Weak routine", "Medium routine", "Strong routine"))
  datasettmp$Status.quo <- mapvalues(datasettmp$Status.quo, from = c("0", "1"), to = c("Not statusquo", "Status quo"))
  
  # metaforest doesn't accept NAs, just blanks
  datasettmp$Exceptional.vs..normal.outcome<-fct_explicit_na(datasettmp$Exceptional.vs..normal.outcome, na_level = "Unclear")
  datasettmp$Outcome.severity<-fct_explicit_na(datasettmp$Outcome.severity, na_level = "Unclear")
  datasettmp$Controllable.vs..uncontrollable.Events<-fct_explicit_na(datasettmp$Controllable.vs..uncontrollable.Events, na_level = "Unclear")
  datasettmp$Action.inaction<-fct_explicit_na(datasettmp$Action.inaction, na_level = "Unclear")
  datasettmp$Routine.strength.recoded<-fct_explicit_na(datasettmp$Routine.strength.recoded, na_level = "Unclear")
  datasettmp$Status.quo<-fct_explicit_na(datasettmp$Status.quo, na_level = "Unclear")
  
  
  #Conduct random-effects weighted MetaForest analysis
  set.seed(42)
  mf.random <- MetaForest(formula = yi ~ DV.group +
                            Exceptional.vs..normal.outcome + 
                            Outcome.severity +
                            Action.inaction +
                            Controllable.vs..uncontrollable.Events+
                            Published +
                            Routine.strength.recoded +
                            Status.quo, 
                          data = datasettmp,
                          whichweights = "random",
                          num.trees = 500, 
                          method = "REML")
  
  summary(mf.random)
  plot(mf.random)
  VarImpPlot(mf.random)
  
  #Univariate partial dependence plot
  PartialDependence(mf.random, vars = "DV.group")
  PartialDependence(mf.random, vars = "Exceptional.vs..normal.outcome")
  PartialDependence(mf.random, vars = "Outcome.severity")
  PartialDependence(mf.random, vars = "Action.inaction")
  PartialDependence(mf.random, vars = "Controllable.vs..uncontrollable.Events")
  PartialDependence(mf.random, vars = "Published")
  PartialDependence(mf.random, vars = "Routine.strength.recoded")
  PartialDependence(mf.random, vars = "Status.quo")
  
  WeightedScatter(datasettmp, yi="yi",vars = c("DV.group",
                                               "Exceptional.vs..normal.outcome",
                                               "Outcome.severity",
                                               "Action.inaction",
                                               "Controllable.vs..uncontrollable.Events",
                                               "Published",
                                               "Routine.strength.recoded",
                                               "Status.quo"))
  WeightedScatter(datasettmp, yi="Effect.size" ,vars = c("DV.group"))
  WeightedScatter(datasettmp, yi="Effect.size" ,vars = c("Exceptional.vs..normal.outcome"))               
  WeightedScatter(datasettmp, yi="Effect.size" ,vars = c("Outcome.severity"))
  WeightedScatter(datasettmp, yi="Effect.size" ,vars = c("Action.inaction"))
  WeightedScatter(datasettmp, yi="Effect.size" ,vars = c("Controllable.vs..uncontrollable.Events"))
  WeightedScatter(datasettmp, yi="Effect.size" ,vars = c("Published"))
  WeightedScatter(datasettmp, yi="Effect.size" ,vars = c("Routine.strength.recoded"))
  WeightedScatter(datasettmp, yi="Effect.size" ,vars = c("Status.quo"))
  ```
  
  ### DV categories as a moderator ###
  
  # no big difference in the effect size between types of DV
  expmoderatorsdv <- rma(yi, vi, method = "REML", data = dataset, 
                         slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                         subset = (IV1.manipulated.=="Manipulated"), 
                         ni = N.sample.size..post.attrition.,
                         mods = ~ DV.group)
  expmoderatorsdv
  
  #### experimental studies / regret ####
  
  # experimental studies / regret
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$DV.group == "1 - regret")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  ###
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expregret <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expregret
  
  forest(expregret, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expregret$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expregret$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expregret$k+2, "Sample size",                  cex=.8)
  par(op)
  
  #### experimental studies / counterfactuals + avoidability ####
  
  # experimental studies / counterfactuals + avoidability
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & ((dataset$DV.group == "2 - counterfactuals") | (dataset$DV.group == "5 - avoidability"))),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expcounterfactuals <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expcounterfactuals
  
  forest(expcounterfactuals, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expcounterfactuals$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expcounterfactuals$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expcounterfactuals$k+2, "Sample size",                  cex=.8)
  par(op)
  
  
  #### experimental studies / compensation ####
  
  # experimental studies / 3 - compensation
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$DV.group == "3 - compensation")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expcompensation <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expcompensation
  
  forest(expcompensation, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expcompensation$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expcompensation$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expcompensation$k+2, "Sample size",                  cex=.8)
  par(op)
  
  #### experimental studies / 4 - self blame ####
  
  # experimental studies / 4 - self blame
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$DV.group == "4 - self blame")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expblame <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expblame
  
  forest(expblame, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expblame$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expblame$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expblame$k+2, "Sample size",                  cex=.8)
  par(op)
  
  #### experimental studies / 6 - offender punishment ####
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$DV.group == "7 - offender punishment")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expoffenderpunishment <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expoffenderpunishment
  
  forest(expoffenderpunishment, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expoffenderpunishment$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expoffenderpunishment$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expoffenderpunishment$k+2, "Sample size",                  cex=.8)
  par(op)
  
  ### Comparing the different DV type moderator effects for experimental studies (more accurate, since collapsed effects) ###
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(expregret$b, expcounterfactuals$b, expcompensation$b, expblame$b, expoffenderpunishment$b)), 
    stderror = as.numeric(c(expregret$se, expcounterfactuals$se, expcompensation$se, expblame$se, expoffenderpunishment$se)), 
    numstudies = as.numeric(c(expregret$k, expcounterfactuals$k,expcompensation$k, expblame$k, expoffenderpunishment$k)),
    meta = as.character(c("Regret DV","Counterfactuals DV", "Compensation DV", "Blame DV","Offender Punishment DV")), 
    tau2 = as.numeric(c(c(expregret$tau2, expcounterfactuals$tau2, expcompensation$tau2, expblame$tau2, expoffenderpunishment$tau2))))
  comparedvtypes
  
  comparedvs <- rma(estimate, sei=stderror, mods = ~ factor(meta), method="FE", data=comparedvtypes, digits=3)
  comparedvs
  
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(expregret$b, expcounterfactuals$b, expcompensation$b, expblame$b, expoffenderpunishment$b)), 
    stderror = as.numeric(c(expregret$se, expcounterfactuals$se, expcompensation$se, expblame$se, expoffenderpunishment$se)), 
    numstudies = as.numeric(c(expregret$k, expcounterfactuals$k,expcompensation$k, expblame$k, expoffenderpunishment$k)),
    meta = as.character(c("1 - Regret DV","2 - Counterfactuals DV", "4 - Compensation DV", "3 - Blame DV","5 - Offender Punishment DV")), 
    tau2 = as.numeric(c(c(expregret$tau2, expcounterfactuals$tau2, expcompensation$tau2, expblame$tau2, expoffenderpunishment$tau2))))
  comparedvtypes
  
  comparedvs2 <- rma(estimate, sei=stderror, mods = ~ factor(meta), method="FE", data=comparedvtypes, digits=3)
  comparedvs2
  
  ## Results for compared studies ##
  
  #moderators
  # all three theoretical moderators found significant moderators
  expmoderators <- rma(yi, vi, method = "REML", data = dataset, 
                       slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                       subset = (IV1.manipulated.=="Measured/compared"), 
                       ni = N.sample.size..post.attrition.,
                       mods = ~ Action.inaction + 
                         factor(Routine.strength.recoded) +
                         Controllable.vs..uncontrollable.Events)
  
  expmoderators
  
  ### DV categories as a moderator ###
  
  # no big differnece in the effect size between types of DV
  expmoderatorsdv <- rma(yi, vi, method = "REML", data = dataset, 
                         slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                         subset = (IV1.manipulated.=="Measured/compared"), 
                         ni = N.sample.size..post.attrition.,
                         mods = ~ DV.group)
  expmoderatorsdv
  
  #### compared studies / regret ####
  
  # compared studies / regret
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Measured/compared") & (dataset$DV.group == "1 - regret")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expregret <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expregret
  
  forest(expregret, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expregret$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expregret$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expregret$k+2, "Sample size",                  cex=.8)
  par(op)
  
  #### compared studies / counterfactuals + avoidability ####
  
  # compared studies / counterfactuals + avoidability
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Measured/compared") & ((dataset$DV.group == "2 - counterfactuals") | (dataset$DV.group == "5 - avoidability"))),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expcounterfactuals <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expcounterfactuals
  
  forest(expcounterfactuals, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expcounterfactuals$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expcounterfactuals$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expcounterfactuals$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # Comparing the different DV type moderator effects for compared studies (more accurate, since collapsed effects)
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(expregret$b, expcounterfactuals$b)), 
    stderror = as.numeric(c(expregret$se, expcounterfactuals$se)),
    numstudies = as.numeric(c(expregret$k, expcounterfactuals$k)),
    meta = as.character(c("Regret DV","Counterfactuals DV")), 
    tau2 = as.numeric(c(c(expregret$tau2, expcounterfactuals$tau2))))
  comparedvtypes
  
  comparedvs <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=comparedvtypes, digits=3)
  comparedvs
  
  # Our theoretical moderators
  # Exceptional.vs..normal.outcome as a moderator (experimental)
  
  expmoderatorsexceptional <- rma(yi, vi, method = "REML", data = dataset, 
                                  slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                                  subset = (IV1.manipulated.=="Manipulated"), 
                                  ni = N.sample.size..post.attrition.,
                                  mods = ~ Exceptional.vs..normal.outcome)
  expmoderatorsexceptional
  # Exceptional outcome / experimental
  # exp studies / action
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Exceptional.vs..normal.outcome == "1")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expexceptionaloutcome <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expexceptionaloutcome
  forest(expexceptionaloutcome, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expexceptionaloutcome$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expexceptionaloutcome$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expexceptionaloutcome$k+2, "Sample size",                  cex=.8)
  par(op)
  # Routine studies / experimental
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Exceptional.vs..normal.outcome == "0")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  exproutineoutcome <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  exproutineoutcome
  
  forest(exproutineoutcome, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, exproutineoutcome$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, exproutineoutcome$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, exproutineoutcome$k+2, "Sample size",                  cex=.8)
  par(op)
  
  ### Comparing the exceptional/routine outcome moderator effects for exp studies ###
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(expexceptionaloutcome$b, exproutineoutcome$b)), 
    stderror = as.numeric(c(expexceptionaloutcome$se, exproutineoutcome$se)),
    numstudies = as.numeric(c(expexceptionaloutcome$k, exproutineoutcome$k)),
    meta = as.character(c("Exceptional outcome ","Routine outcome")), 
    tau2 = as.numeric(c(c(expexceptionaloutcome$tau2, exproutineoutcome$tau2))))
  comparedvtypes
  
  compareseverity <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=comparedvtypes, digits=3)
  compareseverity
  
  # Outcome.severity as a moderator (experimental) ###
  
  expmoderatorsseverity <- rma(yi, vi, method = "REML", data = dataset, 
                               slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                               subset = (IV1.manipulated.=="Manipulated"), 
                               ni = N.sample.size..post.attrition.,
                               mods = ~ Outcome.severity)
  expmoderatorsseverity
  # high sevetity / experimental
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Outcome.severity == "1")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expseverity1 <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  sum(dataset2$N.sample.size..post.attrition.)
  expseverity1
  
  forest(expseverity1, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expseverity1$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expseverity1$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expseverity1$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # low severity studies / experimental
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Outcome.severity == "0")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expseverity0 <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  sum(dataset2$N.sample.size..post.attrition.)
  expseverity0
  
  forest(expseverity0, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expseverity0$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expseverity0$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expseverity0$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # Comparing the severity moderator effects for exp studies
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(expseverity1$b, expseverity0$b)), 
    stderror = as.numeric(c(expseverity1$se, expseverity0$se)),
    numstudies = as.numeric(c(expseverity1$k, expseverity0$k)),
    meta = as.character(c("High severity studies","Low severity studies")), 
    tau2 = as.numeric(c(c(expseverity1$tau2, expseverity0$tau2))))
  comparedvtypes
  
  compareseverity <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=comparedvtypes, digits=3)
  compareseverity
  
  # Action/inaction as a moderator (experimental)
  
  expmoderatorsaction <- rma(yi, vi, method = "REML", data = dataset, 
                             slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                             subset = (IV1.manipulated.=="Manipulated"), 
                             ni = N.sample.size..post.attrition.,
                             mods = ~ Action.inaction)
  expmoderatorsaction
  
  # action studies / experimental
  
  # exp studies / action
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Action.inaction == "1")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expaction <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expaction
  
  forest(expaction, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expaction$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expaction$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expaction$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # inaction studies / experimental
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Action.inaction == "0")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expinaction <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expinaction
  
  forest(expinaction, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expinaction$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expinaction$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expinaction$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # Comparing the moderator effects for exp studies
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(expaction$b, expinaction$b)), 
    stderror = as.numeric(c(expaction$se, expinaction$se)),
    numstudies = as.numeric(c(expaction$k, expinaction$k)),
    meta = as.character(c("Action studies","Inaction studies")), 
    tau2 = as.numeric(c(c(expaction$tau2, expinaction$tau2))))
  comparedvtypes
  
  compareactioninaction <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=comparedvtypes, digits=3)
  compareactioninaction
  
  # Controllable.vs..uncontrollable.Events as a moderator (experimental)
  
  expmoderatorscontrol <- rma(yi, vi, method = "REML", data = dataset, 
                              slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                              subset = (IV1.manipulated.=="Manipulated"), 
                              ni = N.sample.size..post.attrition.,
                              mods = ~ Controllable.vs..uncontrollable.Events)
  expmoderatorscontrol
  
  # controllable studies / experimental
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Controllable.vs..uncontrollable.Events == "1")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expcontrol <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expcontrol
  
  forest(expcontrol, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expcontrol$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expcontrol$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expcontrol$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # uncontrollable studies / experimental
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Controllable.vs..uncontrollable.Events == "0")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  ###
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expuncontrol <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expuncontrol
  
  forest(expuncontrol, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expuncontrol$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expuncontrol$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expuncontrol$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # Comparing the moderator for exp studies
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(expcontrol$b, expuncontrol$b)), 
    stderror = as.numeric(c(expcontrol$se, expuncontrol$se)),
    numstudies = as.numeric(c(expcontrol$k, expuncontrol$k)),
    meta = as.character(c("Controllable studies","Uncontrollable studies")), 
    tau2 = as.numeric(c(c(expcontrol$tau2, expuncontrol$tau2))))
  comparedvtypes
  
  comparecontrol <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=comparedvtypes, digits=3)
  comparecontrol
  
  # Published as a moderator (experimental)
  
  expmoderatorspublished <- rma(yi, vi, method = "REML", data = dataset, 
                                slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                                subset = (IV1.manipulated.=="Manipulated"), 
                                ni = N.sample.size..post.attrition.,
                                mods = ~ Published)
  expmoderatorspublished
  
  # published studies / experimental studies
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Published=="Yes")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  exppublished <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  exppublished
  
  forest(exppublished, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, exppublished$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, exppublished$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, exppublished$k+2, "Sample size",                  cex=.8)
  par(op)
  
  #### unpublished studies / all studies  ####
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Published=="No")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expunpublished <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expunpublished
  
  forest(expunpublished, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expunpublished$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expunpublished$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expunpublished$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # Comparing the moderator for experimental studies
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(exppublished$b, expunpublished$b)), 
    stderror = as.numeric(c(exppublished$se, expunpublished$se)),
    numstudies = as.numeric(c(exppublished$k, expunpublished$k)),
    meta = as.character(c("Published studies","Unpublished")), 
    tau2 = as.numeric(c(c(exppublished$tau2, expunpublished$tau2))))
  comparedvtypes
  
  comparepublished <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=comparedvtypes, digits=3)
  comparepublished
  
  # Published as a moderator (Measured/compared)
  
  expmoderatorspublished <- rma(yi, vi, method = "REML", data = dataset, 
                                slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                                subset = (IV1.manipulated.=="Measured/compared"), 
                                ni = N.sample.size..post.attrition.,
                                mods = ~ Published)
  expmoderatorspublished
  
  # published studies / Measured/compared studies
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Measured/compared") & (dataset$Published=="Yes")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  exppublished <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  exppublished
  
  forest(exppublished, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, exppublished$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, exppublished$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, exppublished$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # unpublished studies / Measured/compared studies
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Measured/compared") & (dataset$Published=="No")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expunpublished <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expunpublished
  
  forest(expunpublished, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expunpublished$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expunpublished$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expunpublished$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # Comparing the moderator for Measured/compared studies
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(exppublished$b, expunpublished$b)), 
    stderror = as.numeric(c(exppublished$se, expunpublished$se)),
    numstudies = as.numeric(c(exppublished$k, expunpublished$k)),
    meta = as.character(c("Published studies","Unpublished")), 
    tau2 = as.numeric(c(c(exppublished$tau2, expunpublished$tau2))))
  comparedvtypes
  
  comparepublished <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=comparedvtypes, digits=3)
  comparepublished
  
  # Moderator Routine Strength###
  # moderators
  # all three theoretical moderators found significant moderators
  expmoderators <- rma(yi, vi, method = "REML", data = dataset, 
                       slab = paste0(dataset$Article, "/S", dataset$Study.., "-", dataset$Sample.., "-", dataset$DV..),
                       subset = (IV1.manipulated.=="Manipulated"), 
                       ni = N.sample.size..post.attrition.,
                       mods = ~ Routine.strength.recoded)
  
  expmoderators
  # Strong routine studies / experimental
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Routine.strength.recoded == "3")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expstrongroutine <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expstrongroutine
  
  forest(expstrongroutine, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expstrongroutine$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expstrongroutine$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expstrongroutine$k+2, "Sample size",                  cex=.8)
  par(op)
  
  #  Weak routine studies / experimental
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Routine.strength.recoded == "2" | dataset$Routine.strength.recoded == "1")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expweakroutine <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expweakroutine
  
  forest(expweakroutine, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expweakroutine$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expweakroutine$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expweakroutine$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # Comparing the moderator for exp studies
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(expstrongroutine$b, expweakroutine$b)), 
    stderror = as.numeric(c(expstrongroutine$se, expweakroutine$se)),
    numstudies = as.numeric(c(expstrongroutine$k, expweakroutine$k)),
    meta = as.character(c("Strong routine studies","weak routine studies")), 
    tau2 = as.numeric(c(c(expstrongroutine$tau2, expweakroutine$tau2))))
  comparedvtypes
  
  compareroutinestrenght <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=comparedvtypes, digits=3)
  compareroutinestrenght
  
  # Moderator Status Quo
  
  # Status quo studies / experimental
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Status.quo == "1")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  expstatusquo <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expstatusquo
  
  forest(expstatusquo, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expstatusquo$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expstatusquo$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expstatusquo$k+2, "Sample size",                  cex=.8)
  par(op)
  
  #  Non status quo studies / experimental
  
  dataset1 <- dataset[which((dataset$IV1.manipulated.=="Manipulated") & (dataset$Status.quo == "0")),]
  
  # collapsing first for manipulated studies
  dataset1$articlestudy <- paste(dataset1$Article, "/", dataset1$Study.., "/", dataset1$Sample..)
  dataset1$articlestudy <- factor(dataset1$articlestudy)
  
  collapsed <- agg(id=articlestudy, es=yi, var=vi, cor=.5,
                   method="BHHR", mod=NULL, data=dataset1)
  # str(collapsed)
  
  # if we collapse the effects within each study
  byarticlestudy <- 
    cbind(
      aggregate(
        yi ~ articlestudy, 
        data=dataset1, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x), count=length(x))),
      aggregate(
        N.sample.size..post.attrition. ~ articlestudy, 
        data=dataset1, FUN = max)[2]
    )
  names(collapsed)[names(collapsed) == 'id'] <- 'articlestudy'
  dataset2 <- merge(byarticlestudy, collapsed,by="articlestudy")
  
  
  expnonstatusquo <- rma(es, var, method = "REML", data = dataset2, slab = articlestudy)
  expnonstatusquo
  
  forest(expnonstatusquo, alim=c(-2,3), xlim=c(-4,4), ilab=dataset2$N.sample.size..post.attrition., ilab.xpos=-2, xlab = "Intrapersonal exception/routine", cex=.8)
  op <- par(font=4)
  # the columns are from -2 to 2 (xlim)
  text(-4, expnonstatusquo$k+2, "Author(s), Year, and Study #",    pos=4, cex=.8)
  text( 4, expnonstatusquo$k+2, "Observed [95% CI]",  pos=2, cex=.8)
  text(-2, expnonstatusquo$k+2, "Sample size",                  cex=.8)
  par(op)
  
  # Comparing the moderator for exp studies
  
  comparedvtypes <- data.frame(
    estimate = as.numeric(c(expstatusquo$b, expnonstatusquo$b)), 
    stderror = as.numeric(c(expstatusquo$se, expnonstatusquo$se)),
    numstudies = as.numeric(c(expstatusquo$k, expnonstatusquo$k)),
    meta = as.character(c("Status quo studies","Non status quo studies")), 
    tau2 = as.numeric(c(c(expstatusquo$tau2, expnonstatusquo$tau2))))
  comparedvtypes
  
  compareroutinestrenght <- rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=comparedvtypes, digits=3)
  compareroutinestrenght
  
  # write.csv(dataset, "past-behavior-meta-coding-outfile.csv")