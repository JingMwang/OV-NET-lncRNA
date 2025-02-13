rm(list = ls())
geneList = read.table("./riskDiff.txt", sep= "\t", header = T)
##----------------------基于基因预测miRNA--------------------

hubGenes = geneList$gene
write.table(hubGenes,"./geneList.txt", sep= "\t", row.names = F, col.names = F, quote = F)

library(multiMiR)

gene2mir <- get_multimir(org     = 'hsa',
                         target  = hubGenes,
                         table   = 'validated',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)

ez = gene2mir@data;dim(ez)

table(ez$support_type)

m_miR = ez[,c(3,4)]
m_miR = m_miR[!duplicated(m_miR),]
colnames(m_miR) = c('miRNAname','geneName')
m_miR$type= 'miRNA'
length(unique(m_miR$miRNAname))
#73


##----------------------基于miRNA预测lncRNA--------------------
miRNAs = unique(ez$mature_mirna_id)

# 下载starBase
# curl 'http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=lncRNA&miRNA=all&clipExpNum=0&degraExpNum=0&pancancerNum=0&programNum=0&program=None&target=all&cellType=all' > starBaseV3_hg19_CLIP-seq_lncRNA_all.txt &

starbase = data.table::fread("d:/PreviousWork/project/annotation_files/starBaseV3_hg19_CLIP-seq_lncRNA_all.txt")
dim(starbase)


# gtf = rtracklayer::import('~/project/annotation_files/gencode.v35.annotation.gtf')
# save(gtf,file = '~/project/annotation_files/annotation_genecodev35.Rdata')
load('d:/PreviousWork/project/annotation_files/annotation_genecodev35.Rdata')
gtf = as.data.frame(gtf)
anno = gtf[,10:13]
lnc_anno = anno[anno$gene_type == "lncRNA",]
lnc_anno = lnc_anno[!duplicated(lnc_anno$gene_name),]

p1 = starbase$geneName %in% lnc_anno$gene_name;table(p1)
starbase = starbase[p1,];dim(starbase)

lnc_mi = starbase[starbase$miRNAname %in% miRNAs,]
length(unique(lnc_mi$miRNAname))


p2 = lnc_mi$pancancerNum >10  &lnc_mi$clipExpNum>4;table(p2)

p3 = lnc_mi$pancancerNum >10  &lnc_mi$clipExpNum>4 & lnc_mi$degraExpNum >0;table(p3)

result1 = lnc_mi[p2,c(2,4)]

mi_result = result1[!duplicated(result1),]
mi_result$type = 'lncRNA'

m_miR = m_miR[m_miR$miRNAname %in% mi_result$miRNAname,]
cyto_input_ceRNA = rbind(m_miR,mi_result)

write.table(cyto_input_ceRNA,file = './mRNA_miRNA_lncRNA_network_20221228.txt',row.names = F,quote = F)


##----------------------基于基因预测TFs转录因子--------------------

# https://maayanlab.cloud/chea3/

# 利用CHEA3数据库预测基因的转录因子
# 选择meanRank top10的转录因子
options(stringsAsFactors = F)
tf_df = read.csv('./TFs_targetGenes.tsv',header = T,sep = '\t')
tfs = tf_df[c(1:10),c(3,6)]

library(stringr)

genes = str_split(tfs$Overlapping_Genes,',')

ll = c(rep(NA,10))
for(i in 1:nrow(tfs)){
  ll[i] = length(genes[[i]])
}

tf_result = data.frame(TFs = c(rep(tfs$TF[1],ll[1]),rep(tfs$TF[2],ll[2]),
                               rep(tfs$TF[3],ll[3]),rep(tfs$TF[4],ll[4]),
                               rep(tfs$TF[5],ll[5]),rep(tfs$TF[6],ll[6]),
                               rep(tfs$TF[7],ll[7]),rep(tfs$TF[8],ll[8]),
                               rep(tfs$TF[9],ll[9]),rep(tfs$TF[10],ll[10])),
                       geneName = unlist(genes))

tf_result$type = 'TF'
write.table(tf_result,file = './Rdata/TF_mRNA_network.txt',quote = F,row.names = F)


##---------------------------

tf_mRNA = tf_result
colnames(tf_mRNA) = colnames(cyto_input_ceRNA)

cyto_input_ceRNA = rbind(cyto_input_ceRNA,tf_mRNA)

write.csv(cyto_input_ceRNA,file = './results/cyto_input_ceRNA_network.csv',row.names = F,quote=F)


mirs = unique(m_miR$miRNAname)
lncs = unique(mi_result$geneName)
tfs = unique(tf_result$TFs)

node_info = data.frame(name = c(hubGenes, mirs,lncs, tfs), 
                       type = c(rep("mRNA",length(hubGenes)),
                                rep("miRNA",length(mirs)),
                                rep("lncRNA",length(lncs)),
                                rep("TF",length(tfs))))
write.csv(node_info, file = "./results/node_info.csv", row.names = F, quote=F)
