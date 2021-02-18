library(ballgown)
setwd("/data/zhendi/protocol/ballgown")
bg = ballgown(dataDir="/data/zhendi/protocol/ballgown/homo_test_dir", samplePattern='SRR', meas='all')
save(bg, file='homo_bg.rda')

bg = ballgown(dataDir="/data/zhendi/protocol/ballgown/mm10_test_dir", samplePattern='SRR', meas='all')
save(bg, file='mm10_bg.rda')
