#om_assemble -HeadMat fsavg.geom fsavg.cond fsavghm.mat

#om_minverser fsavghm.mat fsavghminv.mat 

om_assemble -DSM fsavg.geom fsavg.cond meg114sources.dip fsavgcortex.mat 

om_assemble -H2MM fsavg.geom fsavg.cond subject00MEG.patches s00cortexH2MM.mat

om_assemble -DS2MM meg114sources.dip subject00MEG.patches s00cortexMeg.mat

om_gain -MEG fsavghminv.mat fsavgcortex.mat s00cortexH2MM.mat s00cortexMeg.mat s00MEGleadfield_dip.mat 
