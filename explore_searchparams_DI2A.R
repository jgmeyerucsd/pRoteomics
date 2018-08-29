
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

### PRECURSOR TOLERANCES
# 60k, 3mz iso
mso60k3i_1p20f<-peplvlfdr(msplitresults="mso1p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_1p5p20f<-peplvlfdr(msplitresults="mso1p5p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_2p20f<-peplvlfdr(msplitresults="mso2p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_2p5p20f<-peplvlfdr(msplitresults="mso2p5p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_3.50f<-peplvlfdr(msplitresults="mso3p5p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_420f<-peplvlfdr(msplitresults="mso4p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_4.5p20f<-peplvlfdr(msplitresults="mso4p5p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i_5p20f<-peplvlfdr(msplitresults="mso5p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i<-peplvlfdr(msplitresults="mso7p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR
mso60k3i<-peplvlfdr(msplitresults="mso5p5p20f_DIA2_OT_60k_3i1o_1.txt") ## < 1% FDR


# 60k, 0.8 mz iso
mso60k<-peplvlfdr(msplitresults="mso1p2p20f_DIA2_OT_60k.txt") ## < 1% FDR
mso60k<-peplvlfdr(msplitresults="mso1p6p20f_DIA2_OT_60k.txt") ## < 1% FDR
mso60k<-peplvlfdr(msplitresults="mso2p4p20f_DIA2_OT_60k.txt") ## < 1% FDR


# 50k, 3mz iso
mso50k3i<-peplvlfdr(msplitresults="mso3p5_20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i<-peplvlfdr(msplitresults="mso0p0_20f_DIA2_OT_50k_3i1o_1.txt") ## 5.0 typo
mso50k3i<-peplvlfdr(msplitresults="mso6p5_20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
# repeat 3 more
mso50k3i<-peplvlfdr(msplitresults="mso5p5p_20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i<-peplvlfdr(msplitresults="mso6p0p_20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i<-peplvlfdr(msplitresults="mso7p0p_20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR


mso50k3i<-peplvlfdr(msplitresults="mso4p20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i<-peplvlfdr(msplitresults="mso4p5p20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i<-peplvlfdr(msplitresults="mso5p20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i<-peplvlfdr(msplitresults="mso5p5p20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR
mso50k3i<-peplvlfdr(msplitresults="mso6p20f_DIA2_OT_50k_3i1o_1.txt") ## < 1% FDR


######################################################33
#############3 NO FAIMS
#############################################################

mso15k3mz<-peplvlfdr(msplitresults="D:/FAIMS/20180823_noFAIMS/mso_3p_20f_noFAIMS_OT_15k_3i1o_1.txt") ## < 1% FDR
mso30k3mz<-peplvlfdr(msplitresults="D:/FAIMS/20180823_noFAIMS/mso_3p_20f_noFAIMS_OT_30k_3i1o_1.txt")

mso15kp8mz<-peplvlfdr(msplitresults="D:/FAIMS/20180823_noFAIMS/mso_0p8p_20f_noFAIMS_OT_15k_p8ip4o_1.txt") ## < 1% FDR
mso30kp8mz<-peplvlfdr(msplitresults="D:/FAIMS/20180823_noFAIMS/mso_0p8p_20f_noFAIMS_OT_30k_p8ip4o_1.txt") ## < 1% FDR
mso50kp8mz<-peplvlfdr(msplitresults="D:/FAIMS/20180823_noFAIMS/mso_0p8p_20f_noFAIMS_OT_50k_p8ip4o_1.txt") ## < 1% FDR
mso60kp8mz<-peplvlfdr(msplitresults="D:/FAIMS/20180823_noFAIMS/mso_0p8p_20f_noFAIMS_OT_60k_p8ip4o_1.txt") ## < 1% FDR

mso60kp8mz<-peplvlfdr(msplitresults="D:/FAIMS/20180822_DI2A/mso0p8p_500f_DIA2_ITp8_30m_p1f.txt")
msoitp8mz<-peplvlfdr(msplitresults="D:/FAIMS/20180822_DI2A/mso0p8p_500f_DIA2_ITp8.txt")

