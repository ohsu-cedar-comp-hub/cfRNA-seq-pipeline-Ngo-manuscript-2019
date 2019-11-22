library("data.table")
library("qvalue")

stats_table <- snakemake@input[['stats_table']]
cat(sprintf(c('stats table: ', stats_table, '\n')))


qplot <- snakemake@output[['qplot']]
cat(sprintf(c('Qvalue Output: ', qplot, '\n')))

qhist <- snakemake@output[['qhist']]
cat(sprintf(c('Qvalue hist Output: ', qhist, '\n')))

out_table = snakemake@output[['table']]
cat(sprintf(c('Summary results table', out_table,'\n')))

stats_frame = read.table(stats_table, row.names=1, sep='\t', check.names=F)

qobj = qvalue(p=stats_frame$pvalue, fdr.level=T)

stats_frame$qvalues = qobj$qvalues
stats_frame$lfdr = qobj$lfdr
write.table(as.data.frame(stats_frame), file=out_table, quote=FALSE, sep='\t')

pdf(qplot)
plot(qobj)
dev.off()

pdf(qhist)
hist(qobj)
dev.off()
