
args<-commandArgs(TRUE)

input_file = args[1]
output_file = args[2]

data = read.table(input_file,header=F,quote='', sep='\t')
scores = data[,1]
classes = data[,2]

library(ROCR)

pred = prediction(scores, classes)
perf = performance(pred,"tpr","fpr")
pdf(output_file)
plot(perf)
dev.off()