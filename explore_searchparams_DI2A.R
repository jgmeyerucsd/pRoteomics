
setwd("D:/FAIMS/20180822_DI2A/")
############### PEPTIDE LEVEL FDRs
## 3 dalton isolation, 60k OT resolution
mso60k3i_3p10f<-peplvlfdr(msplitresults="mso3p10f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_3p20f<-peplvlfdr(msplitresults="mso3p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_3p30f<-peplvlfdr(msplitresults="mso3p30f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_3p40f<-peplvlfdr(msplitresults="mso3p40f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_3p50f<-peplvlfdr(msplitresults="mso3p50f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_3p60f<-peplvlfdr(msplitresults="mso3p60f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
## 3 dalton isolation, 15k OT resolution
mso15k3i_3p10f<-peplvlfdr(msplitresults="mso3p10f_DIA2_OT_15k_3i1o_1.txt") ## < 1% FDR
mso15k3i_3p20f<-peplvlfdr(msplitresults="mso3p20f_DIA2_OT_15k_3i1o_1.txt") ## < 1% FDR
mso15k3i_3p30f<-peplvlfdr(msplitresults="mso3p30f_DIA2_OT_15k_3i1o_1.txt") ## < 1% FDR
mso15k3i_3p40f<-peplvlfdr(msplitresults="mso3p40f_DIA2_OT_15k_3i1o_1.txt") ## < 1% FDR
mso15k3i_3p50f<-peplvlfdr(msplitresults="mso3p50f_DIA2_OT_15k_3i1o_1.txt") ## < 1% FDR
mso15k3i_3p60f<-peplvlfdr(msplitresults="mso3p60f_DIA2_OT_15k_3i1o_1.txt") ## < 1% FDR
## p8 dalton isolation, 60k OT resolution
mso60kp8i_p8p10f<-peplvlfdr(msplitresults="msop8p10f_DIA2_OT_60k.txt") ## < 1% FDR
mso60kp8i_p8p20f<-peplvlfdr(msplitresults="msop8p20f_DIA2_OT_60k.txt") ## < 1% FDR
mso60kp8i_p8p30f<-peplvlfdr(msplitresults="msop8p30f_DIA2_OT_60k.txt") ## < 1% FDR
mso60kp8i_p8p40f<-peplvlfdr(msplitresults="msop8p40f_DIA2_OT_60k.txt") ## < 1% FDR
mso60kp8i_p8p50f<-peplvlfdr(msplitresults="msop8p50f_DIA2_OT_60k.txt") ## < 1% FDR
mso60kp8i_p8p60f<-peplvlfdr(msplitresults="msop8p60f_DIA2_OT_60k.txt") ## < 1% FDR
## p8 dalton isolation, 15k OT resolution
mso15kp8i_p8p10f<-peplvlfdr(msplitresults="msop8p10f_DIA2_OT_15k.txt") ## < 1% FDR
mso15kp8i_p8p20f<-peplvlfdr(msplitresults="msop8p20f_DIA2_OT_15k.txt") ## < 1% FDR
mso15kp8i_p8p30f<-peplvlfdr(msplitresults="msop8p30f_DIA2_OT_15k.txt") ## < 1% FDR
mso15kp8i_p8p40f<-peplvlfdr(msplitresults="msop8p40f_DIA2_OT_15k.txt") ## < 1% FDR
mso15kp8i_p8p50f<-peplvlfdr(msplitresults="msop8p50f_DIA2_OT_15k.txt") ## < 1% FDR
mso15kp8i_p8p60f<-peplvlfdr(msplitresults="msop8p60f_DIA2_OT_15k.txt") ## < 1% FDR

#### 30k resolution
## 3 mz isolation
mso30k3i_3p10f<-peplvlfdr(msplitresults="mso3p10f_DIA2_OT_30k_3i1o_1.txt") ## < 1% FDR
mso30k3i_3p20f<-peplvlfdr(msplitresults="mso3p20f_DIA2_OT_30k_3i1o_1.txt") ## < 1% FDR
mso30k3i_3p30f<-peplvlfdr(msplitresults="mso3p30f_DIA2_OT_30k_3i1o_1.txt") ## < 1% FDR
mso30k3i_3p40f<-peplvlfdr(msplitresults="mso3p40f_DIA2_OT_30k_3i1o_1.txt") ## < 1% FDR
mso30k3i_3p50f<-peplvlfdr(msplitresults="mso3p50f_DIA2_OT_30k_3i1o_1.txt") ## < 1% FDR
mso30k3i_3p60f<-peplvlfdr(msplitresults="mso3p60f_DIA2_OT_30k_3i1o_1.txt") ## < 1% FDR
## 0.8 mz isolation
mso30kp8i_p8p10f<-peplvlfdr(msplitresults="msop8p10f_DIA2_OT_30k.txt") ## < 1% FDR
mso30kp8i_p8p20f<-peplvlfdr(msplitresults="msop8p20f_DIA2_OT_30k.txt") ## < 1% FDR
mso30kp8i_p8p30f<-peplvlfdr(msplitresults="msop8p30f_DIA2_OT_30k.txt") ## < 1% FDR
mso30kp8i_p8p40f<-peplvlfdr(msplitresults="msop8p40f_DIA2_OT_30k.txt") ## < 1% FDR
mso30kp8i_p8p50f<-peplvlfdr(msplitresults="msop8p50f_DIA2_OT_30k.txt") ## < 1% FDR
mso30kp8i_p8p60f<-peplvlfdr(msplitresults="msop8p60f_DIA2_OT_30k.txt") ## < 1% FDR


#### 50k resolution
## 3 mz isolation
mso50k3i_3p10f<-peplvlfdr(msplitresults="mso3p10f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i_3p20f<-peplvlfdr(msplitresults="mso3p20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i_3p30f<-peplvlfdr(msplitresults="mso3p30f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i_3p40f<-peplvlfdr(msplitresults="mso3p40f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i_3p50f<-peplvlfdr(msplitresults="mso3p50f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i_3p60f<-peplvlfdr(msplitresults="mso3p60f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
## 0.8 mz isolation
mso50kp8i_p8p10f<-peplvlfdr(msplitresults="msop8p10f_DIA2_OT_50k.txt") ## < 1% FDR
mso50kp8i_p8p20f<-peplvlfdr(msplitresults="msop8p20f_DIA2_OT_50k.txt") ## < 1% FDR
mso50kp8i_p8p30f<-peplvlfdr(msplitresults="msop8p30f_DIA2_OT_50k.txt") ## < 1% FDR
mso50kp8i_p8p40f<-peplvlfdr(msplitresults="msop8p40f_DIA2_OT_50k.txt") ## < 1% FDR
mso50kp8i_p8p50f<-peplvlfdr(msplitresults="msop8p50f_DIA2_OT_50k.txt") ## < 1% FDR
mso50kp8i_p8p60f<-peplvlfdr(msplitresults="msop8p60f_DIA2_OT_50k.txt") ## < 1% FDR
