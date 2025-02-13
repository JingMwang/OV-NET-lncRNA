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
cliFile="Survival_SupplementalTable_S1_20171025_xena_sp"    #临床数据文件
setwd("C:\\biowolf\\cuproptosis\\18.PFS")      #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,c("PFI.time", "PFI")]
cli=na.omit(cli)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365
cli=as.matrix(cli)
row.names(cli)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(cli))

#数据合并
sameSample=intersect(row.names(risk), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], risk[sameSample,"risk",drop=F])

#比较高低风险组的PFS差异,得到显著的pvalue
length=length(levels(factor(rt$risk)))
diff=survdiff(Surv(futime, fustat) ~ risk, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit=survfit(Surv(futime, fustat) ~ risk, data = rt)
#print(surv_median(fit))

#绘制生存曲线
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           font.legend=10,
		           xlab="Time(years)",
		           ylab="Progression free survival",
		           break.time.by = 1,
		           palette = c("red", "blue"),
		           risk.table=TRUE,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)

#输出图形
pdf(file="PFS.pdf", width=6.5, height=5.5, onefile=FALSE)
print(surPlot)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

