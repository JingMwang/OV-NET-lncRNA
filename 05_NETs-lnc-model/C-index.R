######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("rms")
#install.packages("pec")


#引用包
library(dplyr)
library(survival)
library(rms)
library(pec)

riskFile="risk.all.txt"     #风险文件
cliFile="clinical.txt"      #临床数据文件
setwd("C:\\biowolf\\cuproptosis\\22.C-index")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)

#计算C-index值
riskScore=cph(Surv(futime,fustat)~riskScore, data=rt, surv=TRUE)
Age=cph(Surv(futime,fustat)~Age, data=rt, surv=TRUE)
Gender=cph(Surv(futime,fustat)~Gender, data=rt, surv=TRUE)
Grade=cph(Surv(futime,fustat)~Grade, data=rt, surv=TRUE)
Stage=cph(Surv(futime,fustat)~Stage, data=rt, surv=TRUE)
c_index  <- cindex(list("Risk score"=riskScore, 
                        "Age"=Age,
                        "Gender"=Gender,
                        "Grade"=Grade,
                        "Stage"=Stage),
                    formula=Surv(futime,fustat)~ .,
                    data=rt,
                    eval.times=seq(0,10,1),
                    splitMethod="bootcv",
                    B=1000
                    )
#输出图形
pdf(file="C-index.pdf", width=5.5, height=5)
plot(c_index, 
     xlim=c(0,10), ylim=c(0.4,0.8), 
     col=bioCol, xlab="Time (years)",
     legend.x=6, legend.y=0.82, legend.cex=1)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

