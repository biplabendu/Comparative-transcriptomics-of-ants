library(ggplot2)
library(boot)

MES = read.csv("MES.csv", header =T)
attach(MES)
names(MES)
```

#GLM analysis to correlate each module with the four phenotypic traits (caste, queen number, worker sterility and invasiveness) without the genus added effect.

rate.lm <- glm(#module ~ factor(Genus), data = MES)
  summary(rate.lm)
  
  lmfit <- function(data, indices) {
    fit <- lm(#module ~ factor(Genus), data = data[indices, ])
      return(coef(fit))
  }
  
  # Bootstrap results 100 times
  results <- boot(data= MES, statistic = lmfit, R = 1000, strata=MES$Genus)
  
  
  for (i in 2:length(rate.lm$coefficients)) {
    bci <- boot.ci(results, type = "basic", index = i)
    print(sprintf("%s,%.4f,%.4f,%.4f", names(rate.lm$coefficients)[i], results$t0[i],
                  bci$basic[4], bci$basic[5]))
  }
  
  # pvalue
  twosidep<-function(data){
    p1<-sum(data>0)/length(data)
    p2<-sum(data<0)/length(data)
    p<-min(p1,p2)*2
    return(p)
  }
  
  twosidep(results$t[,2])
  twosidep(results$t[,3]) 
  twosidep(results$t[,4])
  twosidep(results$t[,5])
  twosidep(results$t[,6])
  
  # Adjusted corrected pvalue
  
  data = read.csv("pvalue.csv", header =T)
  attach(data)
  names(data)
  
  p.adjust(p, method="fdr", n=length(p))
  