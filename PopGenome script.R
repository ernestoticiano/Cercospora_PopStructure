

library("PopGenome")
ernesto.popgenome <- readData("C:/Users//jseari2/Desktop/Ernesto_playground", format="VCF", include.unknown=TRUE)


setwd("C:/Users/jseari2/Desktop/Ernesto_playground")
ernesto.popgenome <- readData("VCF_Ernesto", format="VCF", include.unknown=TRUE)
ernesto.popgenome <- set.populations(ernesto.popgenome, list(
                                     c("DMCC_SAL_3558","DMCC_SAL_3560","DMCC_SAL_3561","DMCC_SAL_3564","DMCC_SAL_3565","DMCC_SAL_3566","DMCC_SAL_3597","DMCC_SAL_3600","DMCC_SAL_3604","DMCC_SAL_3779"),
                                     c("DMCC_BAR_3649","DMCC_BAR_3654","DMCC_BAR_3657","DMCC_BAR_3658","DMCC_BAR_3660","DMCC_BAR_3665","DMCC_BAR_3666","DMCC_BAR_3670","DMCC_BAR_3726","DMCC_BAR_3773","DMCC_PBH_3546","DMCC_PBH_3581","DMCC_PBH_3584","DMCC_PBH_3752","DMCC_TAR_3549","DMCC_TAR_3588","DMCC_TAR_3591","DMCC_TAR_3595"),
                                     c("DMCC_A_3502","DMCC_A_3507","DMCC_A_3510","DMCC_A_3512","DMCC_A_3514","DMCC_A_3517","DMCC_A_3518","DMCC_A_3525","DMCC_A_3530","DMCC_PBH_3494","DMCC_PBH_3495","DMCC_PBH_3499","DMCC_PBH_3500","DMCC_PBH_3577","DMCC_PBH_3579","DMCC_PBH_3580","DMCC_PBH_3700","DMCC_PBH_3723","DMCC_PSJ_3468","DMCC_PSJ_3475","DMCC_PSJ_3486","DMCC_PSJ_3487","DMCC_PSJ_3682","DMCC_PSJ_3706","DMCC_PSJ_3708","DMCC_PSJ_3711","DMCC_PSJ_3713"),
                                     c("DMCC_AM_3613","DMCC_AM_3614","DMCC_AM_3615","DMCC_AM_3616","DMCC_AM_3621","DMCC_AM_3622","DMCC_AM_3615","DMCC_AM_3632","DMCC_AM_3788","DMCC_AM_3789"),
                                     c("DMCC_CM_3533","DMCC_CM_3538","DMCC_CM_3550","DMCC_CM_3551","DMCC_CM_3556","DMCC_CM_3568","DMCC_CM_3569","DMCC_CM_3573","DMCC_CM_3777","DMCC_CM_3794"),
                                     c("DMCC_KM_3767","DMCC_TN_3995","DMCC_TN_3996","DMCC_TN_4003","DMCC_TN_4006","DMCC_TN_4011","DMCC_TN_4018","DMCC_TN_4021","DMCC_TN_4023","DMCC_TN_4026","DMCC_TN_4028"),
                                     c("DMCC_BTX_3607","DMCC_BTX_3608","DMCC_BTX_3612","DMCC_BTX_3639","DMCC_BTX_3640","DMCC_BTX_3641","DMCC_BTX_3642","DMCC_BTX_3646","DMCC_BTX_3676","DMCC_BTX_3791")
))
#population
#alabama:  1 
#arkansas  2 
#LA  3 
#Missi  4
#Missour  5
#Tenne  6 
#Texas  7



ernesto.popgenome <- concatenate.regions(ernesto.popgenome)

#Fst
ernesto.popgenome <- F_ST.stats(ernesto.popgenome)

get.F_ST(ernesto.popgenome)
ernesto.popgenome@nuc.F_ST.pairwise

#Tajima.D
ernesto.popgenome <- neutrality.stats(ernesto.popgenome)
ernesto.popgenome@Tajima.D


#tree
ernesto.popgenome_tree <- readData("VCF_Ernesto", format="VCF", include.unknown=TRUE)
ernesto.popgenome_tree <- set.populations(ernesto.popgenome_tree, list(
  c("DMCC_SAL_3560","DMCC_SAL_3564","DMCC_SAL_3565","DMCC_SAL_3597","DMCC_SAL_3779","DMCC_BAR_3649","DMCC_BAR_3657","DMCC_BAR_3658","DMCC_BAR_3660","DMCC_BAR_3665", "DMCC_BAR_3666","DMCC_BAR_3670","DMCC_BAR_3726","DMCC_BAR_3773","DMCC_PBH_3546","DMCC_PBH_3584","DMCC_PBH_3752","DMCC_TAR_3549","DMCC_TAR_3588","DMCC_A_3502","DMCC_A_3507","DMCC_A_3510","DMCC_A_3517","DMCC_A_3530","DMCC_PBH_3499","DMCC_PBH_3500","DMCC_PBH_3580","DMCC_PSJ_3468","DMCC_PSJ_3475","DMCC_PSJ_3487","DMCC_PSJ_3682","DMCC_PSJ_3708","DMCC_PSJ_3713","DMCC_AM_3613","DMCC_AM_3614","DMCC_AM_3621","DMCC_AM_3788","DMCC_AM_3789","DMCC_CM_3550","DMCC_CM_3551","DMCC_CM_3568","DMCC_CM_3569","DMCC_CM_3777","DMCC_KM_3767","DMCC_TN_4006","DMCC_TN_4011","DMCC_TN_4018","DMCC_TN_4021","DMCC_TN_4023","DMCC_TN_4028","DMCC_BTX_3639","DMCC_BTX_3640","DMCC_BTX_3646"),
  c("DMCC_SAL_3558","DMCC_SAL_3561","DMCC_SAL_3566","DMCC_SAL_3600","DMCC_SAL_3604","DMCC_PBH_3581","DMCC_TAR_3591","DMCC_A_3512","DMCC_A_3514","DMCC_A_3518","DMCC_A_3525","DMCC_PBH_3494","DMCC_PBH_3495","DMCC_PBH_3577","DMCC_PBH_3579","DMCC_PBH_3700","DMCC_PBH_3723","DMCC_PSJ_3486","DMCC_PSJ_3706","DMCC_PSJ_3711","DMCC_AM_3615","DMCC_AM_3616","DMCC_AM_3622","DMCC_AM_3615","DMCC_AM_3632","DMCC_CM_3533","DMCC_CM_3538","DMCC_CM_3556","DMCC_CM_3573","DMCC_CM_3794","DMCC_TN_3995","DMCC_TN_3996","DMCC_TN_4026","DMCC_BTX_3607","DMCC_BTX_3612","DMCC_BTX_3641","DMCC_BTX_3642","DMCC_BTX_3676","DMCC_BTX_3791"),
  c("DMCC_BAR_3654","DMCC_TAR_3595","DMCC_TN_4003","DMCC_BTX_3608")
))

#Fst
ernesto.popgenome_tree <- concatenate.regions(ernesto.popgenome_tree)
ernesto.popgenome_tree <- F_ST.stats(ernesto.popgenome_tree)


#Tajima.D
ernesto.popgenome_tree <- neutrality.stats(ernesto.popgenome)
ernesto.popgenome_tree@Tajima.






