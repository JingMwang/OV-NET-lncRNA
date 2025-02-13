######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)
corFilter=0.4            #相关系数的过滤标准
pvalueFilter=0.001       #相关性检验p值的过滤标准
setwd("C:\\Users\\27191\\Desktop\\NETs\\cuproptosis\\09.cuproptosisLnc")     #设置工作目录

#读取lncRNA的表达文件,并对数据进行处理
rt=read.table("lncRNA.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目
sampleType=c(rep(1,conNum), rep(2,treatNum))

#读取铜死亡相关基因的表达文件,并对数据进行处理
rt1=read.table("cuproptosisExp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
cuproptosis=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
cuproptosis=avereps(cuproptosis)
cuproptosis=cuproptosis[rowMeans(cuproptosis)>0.1,]

#删掉正常样品
group=sapply(strsplit(colnames(cuproptosis),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
cuproptosis=cuproptosis[,group==0]

#相关性检验
outTab=data.frame()
for(i in row.names(lncRNA)){
	if(sd(lncRNA[i,])>0.1){
		test=wilcox.test(data[i,] ~ sampleType)
		if(test$p.value<0.05){
			for(j in row.names(cuproptosis)){
				x=as.numeric(lncRNA[i,])
				y=as.numeric(cuproptosis[j,])
				corT=cor.test(x,y)
				cor=corT$estimate
				pvalue=corT$p.value
				if((cor>corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(Cuproptosis=j,lncRNA=i,cor,pvalue,Regulation="postive"))
				}
				if((cor< -corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(Cuproptosis=j,lncRNA=i,cor,pvalue,Regulation="negative"))
				}
			}
		}
	}
}

#输出相关性的结果
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)

#提取铜死亡相关lncRNA的表达量
cuproptosisLncRNA=unique(as.vector(outTab[,"lncRNA"]))
cuproptosisLncRNAexp=data[cuproptosisLncRNA,]
cuproptosisLncRNAexp=rbind(ID=colnames(cuproptosisLncRNAexp), cuproptosisLncRNAexp)
write.table(cuproptosisLncRNAexp,file="cuproptosisLncExp.txt",sep="\t",quote=F,col.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

