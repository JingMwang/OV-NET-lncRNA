######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)

riskFile="risk.all.txt"      #风险文件
cliFile="clinical.txt"       #临床数据文件
setwd("C:\\biowolf\\cuproptosis\\24.cliGroupSur")      #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

#数据合并
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

#对临床性状的每个分组进行循环
for(j in names(tab)){
	rt1=rt[(rt[,"clinical"]==j),]
	tab1=table(rt1[,"Risk"])
	tab1=tab1[tab1!=0]
	labels=names(tab1)
	if(length(labels)!=2){next}
	if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
		titleName=paste0("age",j)
	}
	
	#计算高低风险组差异pvalue
	diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	
	#绘制生存曲线
	fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
	surPlot=ggsurvplot(fit, 
			           data=rt1,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           title=paste0("Patients with ",j),
			           legend.title="Risk",
			           legend.labs=labels,
			           font.legend=12,
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette=c("red", "blue"),
			           risk.table=F,
			       	   risk.table.title="",
			           risk.table.col = "strata",
			           risk.table.height=.25)
	
	#输出图片
	j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
	pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
			       width = 6,        #图片的宽度
			       height =5)        #图片的高度
	print(surPlot)
	dev.off()
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

