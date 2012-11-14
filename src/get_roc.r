
args<-commandArgs(TRUE)

input_file = args[1]

data = read.table(input_file,header=F,quote='', sep='\t')
scores = data[,1]
classes = data[,2]

library(ROCR)

pred = prediction(scores, classes)
perf = performance(pred,"tpr","fpr")
pdf("roc_curve.pdf")
plot(perf)
dev.off()