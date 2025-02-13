######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("pheatmap")


library(pheatmap)         #引用包
setwd("C:\\biowolf\\cuproptosis\\19.riskPlot")      #设置工作目录

#定义风险曲线的函数
bioRiskPlot=function(inputFile=null, project=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
	rt=rt[order(rt$riskScore),]      #按照病人风险打分对样品进行排序
		
	#绘制风险曲线
	riskClass=rt[,"risk"]
	lowLength=length(riskClass[riskClass=="low"])
	highLength=length(riskClass[riskClass=="high"])
	lowMax=max(rt$riskScore[riskClass=="low"])
	line=rt[,"riskScore"]
	line[line>10]=10
	pdf(file=paste0(project, ".riskScore.pdf"), width=7, height=4)
	plot(line, type="p", pch=20,
		 xlab="Patients (increasing risk socre)",
		 ylab="Risk score",
		 col=c(rep("blue",lowLength),rep("red",highLength)) )
	abline(h=lowMax,v=lowLength,lty=2)
	legend("topleft", c("High risk","Low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	dev.off()
		
	#绘制生存状态图
	color=as.vector(rt$fustat)
	color[color==1]="red"
	color[color==0]="blue"
	pdf(file=paste0(project, ".survStat.pdf"), width=7, height=4)
	plot(rt$futime, pch=19,
		 xlab="Patients (increasing risk socre)",
		 ylab="Survival time (years)",
		 col=color)
	legend("topleft", c("Dead","Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	abline(v=lowLength,lty=2)
	dev.off()
	
	#定义热图注释的颜色
	ann_colors=list()
	bioCol=c("blue", "red")
	names(bioCol)=c("low", "high")
	ann_colors[["Risk"]]=bioCol

	#绘制风险热图
	rt1=rt[c(3:(ncol(rt)-2))]
	rt1=t(rt1)
	annotation=data.frame(Risk=rt[,ncol(rt)])
	rownames(annotation)=rownames(rt)
	pdf(file=paste0(project, ".heatmap.pdf"), width=7, height=4)
	pheatmap(rt1, 
		     annotation=annotation,
		     annotation_colors = ann_colors, 
		     cluster_cols = FALSE,
		     cluster_rows = FALSE,
		     show_colnames = F,
		     scale="row",
		     color = colorRampPalette(c(rep("blue",3.5), "white", rep("red",3.5)))(50),
		     fontsize_col=3,
		     fontsize=7,
		     fontsize_row=8)
	dev.off()
}

#tarin组风险曲线
bioRiskPlot(inputFile="risk.train.txt", project="train")
#test组风险曲线
bioRiskPlot(inputFile="risk.test.txt", project="test")
#所有样品风险曲线
bioRiskPlot(inputFile="risk.all.txt", project="all")


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

