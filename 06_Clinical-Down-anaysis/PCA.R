######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("scatterplot3d")


#引用包
library(limma)
library(scatterplot3d)
setwd("C:\\biowolf\\cuproptosis\\25.PCA")        #设置工作目录

#定义PCA分析函数
myPCA=function(input=null,output=null){
	#读取表达数据文件
	rt=read.table(input, header=T, sep="\t", check.names=F)
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	data=data[rowMeans(data)>0.5,]
	
	#删除正常样品
	type=sapply(strsplit(colnames(data),"\\-"),"[",4)
	type=sapply(strsplit(type,""),"[",1)
	type=gsub("2","1",type)
	data=t(data[,type==0])
	rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
		
	#读取risk风险文件
	risk=read.table("risk.all.txt", header=T, sep="\t", row.names=1, check.names=F)
	sameSample=intersect(rownames(data),rownames(risk))
	data=data[sameSample,]
	risk=risk[sameSample,]
	group=as.vector(risk[,"risk"])
		
	#PCA分析
	data.class <- rownames(data)
	data.pca <- prcomp(data, scale. = TRUE)

	#绘制PCA图形
	color=ifelse(group=="low",4,2)
	pcaPredict=predict(data.pca)
	pdf(file=output, width=7, height=7)
	par(oma=c(1,1,2.5,1))
	s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=35)
	legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
	dev.off()
}

######绘制所有基因的PCA图，将04节课symbol.txt复制到当前目录
myPCA(input="symbol.txt",output="PCA.allGene.pdf")
######绘制铜死亡相关基因的PCA图，将09节课cuproptosisExp.txt复制到当前目录
myPCA(input="cuproptosisExp.txt",output="PCA.cuproptosisGene.pdf")
######绘制铜死亡相关lncRNA的PCA图，将09节课cuproptosisLncExp.txt复制到当前目录
myPCA(input="cuproptosisLncExp.txt",output="PCA.cuproptosisLncRNA.pdf")


######读取风险文件,绘制模型lncRNA的PCA图，将14节课risk.all.txt复制到当前目录
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])
		
#PCA分析
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)

#可视化
color=ifelse(group=="low",4,2)
pcaPredict=predict(data.pca)
pdf(file="PCA.riskLnc.pdf", width=6.5, height=6)
par(oma=c(1,1,2.5,1))
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=35)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c(4,2))
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

