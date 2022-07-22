#om_assemble -HeadMat fsavg.geom fsavg.cond fsavghm.mat

#om_minverser fsavghm.mat fsavghminv.mat 

om_assemble -SSM fsavg.geom fsavg.cond fsavgSource.tri fsavgcortex.mat

om_assemble -H2MM fsavg.geom fsavg.cond subject00MEG.patches s00cortexH2MM.mat

om_assemble -SS2MM fsavgSource.tri subject00MEG.patches s00cortexMeg.mat

om_gain -MEG fsavghminv.mat fsavgcortex.mat s00cortexH2MM.mat s00cortexMeg.mat s0_MEGleadfield_cortex.mat 
