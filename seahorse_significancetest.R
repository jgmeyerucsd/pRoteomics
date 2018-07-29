
#### mesurement number 12
g1<-c(471.63,
      507.02,
      495.66,
      590.61,
      473.22) ### albumax b4 and during, low gluc
g2<-c(431.48,
      491.51,
      319.23,
      500.62,
      362.44) ### albumax b4 and during, high gluc
g3<-c(364.82,
      399.94,
      246.94,
      245.37,
      324.69) ### albumax b4 and during, high gluc +DCA

f5b<-c(6.4,6.48,7.45,3.54,2.64) ### albumax b4 and during, low gluc
f6b<-c(1.8,5.2,0.39,1.5,0.58) ### albumax b4 and during, low gluc


f5m<-c(9.71,10.75,12.63,6.85,5.49) ### albumax b4 and during, low gluc
f6m<-c(6.87,11.27,3.69,1.60,4.46) ### albumax b4 and during, low gluc
t.test(f5b,f6b,alternative=c("greater"),paired=T)
t.test(f5m,f6m,alternative=c("greater"),paired=T)



g2<-c(431.48,
      491.51,
      319.23,
      500.62,
      362.44) ### albumax b4 and during, high gluc
g3<-c(364.82,
      399.94,
      246.94,
      245.37,
      324.69) ### albumax b4 and during, high gluc +DCA





aov(g1,g2,g3)
t.test(g1,g2)

t.test(g1,g2,alternative=c("greater"),paired=T)
t.test(g1,g3,alternative=c("greater"),paired=T)
t.test(g2,g3,alternative=c("greater"),paired=T)



s1<-c(296.91,
      555.43,
      285.64,
      349.54,
      246.72)
s2<-c(234.43,
      142.02,
      297.78,
      297.97,
      268.26)
s3<-c(270,
      447.27,
      408.18,
      274.80,
      246.01)

t.test(s1,s2,alternative=c("greater"),paired=T)
t.test(s2,s3,alternative=c("less"),paired=T)
t.test(g2,g3,alternative=c("greater"),paired=T)


