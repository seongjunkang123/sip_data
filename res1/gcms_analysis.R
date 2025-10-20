library(erah)
library(rlist)
# import MS library
nist_replib <- importMSP(filename = "C:/Users/there/Desktop/lib_correct/Nist_GC_replib.txt", DB.name = "Nist_replib", 
                         DB.version = "replib", DB.info ='NIST_GC_replib')
nist_mainlib <- importMSP(filename = "C:/Users/there/Desktop/NIST_mainlib/Nist_GC_mainlib.txt", DB.name = "Nist_mainlib", 
                          DB.version = "mainlib", DB.info ='NIST_GC_mainlib')
# separate NIST 2020 mainlib to two files due to the memory limit
nist_main1 <- importMSP(filename = "C:/Users/there/Desktop/lib_correct/main.1.txt", DB.name = "Nist_main1", 
                        DB.version = "main1", DB.info ='NIST_GC_mainpart1')
nist_main2 <- importMSP(filename = "C:/Users/there/Desktop/lib_correct/main.2.txt", DB.name = "Nist_main2", 
                        DB.version = "main2", DB.info ='NIST_GC_mainpart2')

MoNA <- importMSP(filename = "C:/Users/there/Desktop/MoNA-export-GC-MS_Spectra-msp/MoNA-export-GC-MS_Spectra.msp", DB.name = "MoNA_GC", 
                  DB.version = "GC", DB.info ='exportGC')
# set up, tou can jump to load_deconv.rda, the inst.csv could use eRah function: createdt('your_experiment_files__path/') see erah.pdf for more information
createdt("C:/Users/there/Desktop/retest_A")
retest_A <- newExp(instrumental="C:/Users/there/Desktop/retest_A/retest_A_inst.csv",phenotype = "C:/Users/there/Desktop/retest_A/retest_A_pheno.csv"
                    ,info= "B196_retest_A")

createdt("C:/Users/there/Desktop/retest_B")
retest_B <- newExp(instrumental="C:/Users/there/Desktop/retest_B/retest_B_inst.csv",phenotype = "C:/Users/there/Desktop/retest_B/retest_B_pheno.csv"
                   ,info= "B196_retest_B")

createdt("C:/Users/there/Desktop/retest_C")
retest_C <- newExp(instrumental="C:/Users/there/Desktop/retest_C/retest_C_inst.csv",phenotype = "C:/Users/there/Desktop/retest_C/retest_C_pheno.csv"
                   ,info= "B196_retest_C")
#check exp
metaData(retest_A)
metaData(retest_B)
metaData(retest_C)
#deconvolution, using parameter settings in erah.pdf
re_ex1.dec.par <- setDecPar(min.peak.width = 1, 
                            avoid.processing.mz = c(35:69,73:75,147:149))
re_exA <- deconvolveComp(retest_A, re_ex1.dec.par )
re_exB <- deconvolveComp(retest_B, re_ex1.dec.par )
re_exC <- deconvolveComp(retest_C, re_ex1.dec.par )

#setwd("C:/Users/there/Desktop/RE_ABC")

save(re_exA, file = "retest_A_deconv.rda")
save(re_exB, file = "retest_B_deconv.rda")
save(re_exC, file = "retest_C_deconv.rda")

#load deconvolution files, can save time of deconvolution steps
#load("retest_A_deconv.rda")
#load("retest_B_deconv.rda")
#load("retest_C_deconv.rda")
gc()
## Alignment, aliment parameter mz.range is from erah.pdf
re_ex1.al.par <- setAlPar(min.spectra.cor = 0.90, max.time.dist = 3, 
                          mz.range = 70:600)
re_exA2 <- alignComp(re_exA, alParameters = re_ex1.al.par)
re_exB2 <- alignComp(re_exB, alParameters = re_ex1.al.par)
re_exC2 <- alignComp(re_exC, alParameters = re_ex1.al.par)
gc()
re_exA3 <- recMissComp(re_exA2,min.samples = 27)
re_exB3 <- recMissComp(re_exB2,min.samples = 18)
re_exC3 <- recMissComp(re_exC2,min.samples = 17)
gc()
save(re_exA3,file= "retest_A_rec.rda")
save(re_exB3,file= "retest_B_rec.rda")
save(re_exC3,file= "retest_C_rec.rda")
#load("retest_A_rec.rda")
#load("retest_B_rec.rda")
#load("retest_C_rec.rda")
## Identification
re_exA4.1 <- identifyComp(re_exA3,id.database = nist_main1)
re_exA4.1.2 <- identifyComp(re_exA3,id.database = nist_main2)
re_exA4.2 <- identifyComp(re_exA3,id.database = nist_replib)
gc()
re_exB4.1 <- identifyComp(re_exB3,id.database = nist_main1)
re_exB4.1.2 <- identifyComp(re_exB3,id.database = nist_main2)
re_exB4.2 <- identifyComp(re_exB3,id.database = nist_replib)
gc()
re_exC4.1 <- identifyComp(re_exC3,id.database = nist_main1)
re_exC4.1.2 <- identifyComp(re_exC3,id.database = nist_main2)
re_exC4.2 <- identifyComp(re_exC3,id.database = nist_replib)
gc()

save(re_exA4.1,file = 'reid_A_main1.rda')
save(re_exA4.1.2,file = 'reid_A_main2.rda')
save(re_exA4.2,file = 'reid_A_rep.rda')

save(re_exB4.1,file = 'reid_B_main1.rda')
save(re_exB4.1.2,file = 'reid_B_main2.rda')
save(re_exB4.2,file = 'reid_B_rep.rda')

save(re_exC4.1,file = 'reid_C_main1.rda')
save(re_exC4.1.2,file = 'reid_C_main2.rda')
save(re_exC4.2,file = 'reid_C_rep.rda')
# identification result
re_idA1.list <- idList(re_exA4.1,id.database = nist_main1)
re_idA1.2.list <- idList(re_exA4.1.2,id.database = nist_main2)
re_idA2.list <- idList(re_exA4.2,id.database = nist_replib)

re_idB1.list <- idList(re_exB4.1,id.database = nist_main1)
re_idB1.2.list <- idList(re_exB4.1.2,id.database = nist_main2)
re_idB2.list <- idList(re_exB4.2,id.database = nist_replib)

re_idC1.list <- idList(re_exC4.1,id.database = nist_main1)
re_idC1.2.list <- idList(re_exC4.1.2,id.database = nist_main2)
re_idC2.list <- idList(re_exC4.2,id.database = nist_replib)
gc()
#export
re_A_main1_align.list <- alignList(re_exA4.1, by.area=F)
re_A_main2_align.list <- alignList(re_exA4.1.2, by.area=F)
re_A_rep_align.list <- alignList(re_exA4.2, by.area=F)

write.csv(re_A_main1_align.list, file="re_A_alignmain1_list.csv")
write.csv(re_A_main2_align.list, file="re_A_alignmain2_list.csv")
write.csv(re_A_rep_align.list, file="re_A_alignrep_list.csv")
write.csv(re_idA1.list, file="re_A_main1_misscomp_idlist.csv")
write.csv(re_idA1.2.list, file="re_A_main2_misscomp_idlist.csv")
write.csv(re_idA2.list, file="re_A_rep_misscomp_idlist.csv")

re_B_main1_align.list <- alignList(re_exB4.1, by.area=F)
re_B_main2_align.list <- alignList(re_exB4.1.2, by.area=F)
re_B_rep_align.list <- alignList(re_exB4.2, by.area=F)

write.csv(re_B_main1_align.list, file="re_B_alignmain1_list.csv")
write.csv(re_B_main2_align.list, file="re_B_alignmain2_list.csv")
write.csv(re_B_rep_align.list, file="re_B_alignrep_list.csv")
write.csv(re_idB1.list, file="re_B_main1_misscomp_idlist.csv")
write.csv(re_idB1.2.list, file="re_B_main2_misscomp_idlist.csv")
write.csv(re_idB2.list, file="re_B_rep_misscomp_idlist.csv")

re_C_main1_align.list <- alignList(re_exC4.1, by.area=F)
re_C_main2_align.list <- alignList(re_exC4.1.2, by.area=F)
re_C_rep_align.list <- alignList(re_exC4.2, by.area=F)

write.csv(re_C_main1_align.list, file="re_C_alignmain1_list.csv")
write.csv(re_C_main2_align.list, file="re_C_alignmain2_list.csv")
write.csv(re_C_rep_align.list, file="re_C_alignrep_list.csv")
write.csv(re_idC1.list, file="re_C_main1_misscomp_idlist.csv")
write.csv(re_idC1.2.list, file="re_C_main2_misscomp_idlist.csv")
write.csv(re_idC2.list, file="re_C_rep_misscomp_idlist.csv")

gc()
