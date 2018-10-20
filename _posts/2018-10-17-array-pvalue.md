---
layout:     post
title:      "甲基化芯片P值cutoff的选择"
subtitle:   ""
date:       2018-10-17
author:     "Lee"
header-img: "img/post-bg-unix-linux.jpg"
tags:
    - 甲基化
    - 芯片
    - 数据处理
---

### *p*值如何得到

想确定芯片处理时*p*值的cutoff,首先需要知道*p*值是如何得到。下面是minfi包中关于detectionP函数的源码。


{% raw %}

```r
function (rgSet, type = "m+u")
{
  locusNames <- getManifestInfo(rgSet, "locusNames")
  detP <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
    dimnames = list(locusNames, sampleNames(rgSet)))
  controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")

  #以下八行代码用来提取空白对照探针的红绿信号，计算平均值和方差
  r <- getRed(rgSet)
  rBg <- r[controlIdx, ]
  rMu <- matrixStats::colMedians(rBg)
  rSd <- matrixStats::colMads(rBg)
  g <- getGreen(rgSet)
  gBg <- g[controlIdx, ]
  gMu <- matrixStats::colMedians(gBg)
  gSd <- matrixStats::colMads(gBg)

  #获得一类或二类探针的荧光强度
  TypeII <- getProbeInfo(rgSet, type = "II")
  TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
  TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")

  #以背景探针的荧光强度作为0假设，计算p值
  for (i in 1:ncol(rgSet)) {
    intensity <- r[TypeI.Red$AddressA, i] + r[TypeI.Red$AddressB,
      i]
    detP[TypeI.Red$Name, i] <- 1 - pnorm(intensity, mean = rMu[i] *
      2, sd = rSd[i] * 2)
    intensity <- g[TypeI.Green$AddressA, i] + g[TypeI.Green$AddressB,
      i]
    detP[TypeI.Green$Name, i] <- 1 - pnorm(intensity, mean = gMu[i] *
      2, sd = gSd[i] * 2)
    intensity <- r[TypeII$AddressA, i] + g[TypeII$AddressA,
      i]
    detP[TypeII$Name, i] <- 1 - pnorm(intensity, mean = rMu[i] +
      gMu[i], sd = rSd[i] + gSd[i])
  }
  detP
}
```
{% endraw %}

从代码中可以看出计算探针*p*值时，会先对阴性对照探针的信号强度拟合一个正态分布，然后根据这个正态分布计算其它探针信号的*p*值。（这里需要说明的是，这种计算方式是minfi包中采用的，而illumina官方软件BeadStudio是否采用该方法计算*p*值，我们不得而知。）

### cutoff选择0.01合理吗？
很多文章会选择0.01作为甲基化芯片*p*值的cutoff，但是选择0.01作为cutoff真的合理吗？
为了解答这个问题，我们下载[GSE72773](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1870133)数据来进行相关实验。该数据集，没有提供原始idat文件，我们直接下载信号强度和*p*值文件——[GSE72773_datSignal.csv.gz](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72773/suppl/GSE72773%5FdatSignal%2Ecsv%2Egz)。从文件中提取*p*值和beta值。选择前四个点，以p值为横坐标，beta值为纵坐标做散点图，如下：
![plot](/img/in-post/2018-10-17-pvalue.jpg)
可以看到，图中*p*值为0的点聚集度很高，而*p*值小于0.01并大于0的点明显离散.可以认为即使*p*值小于0.01,得出的beta值还是有很大误差。

### 为什么*p*值小于0.01得到的beta值还是不准？
要回答这个问题，我们首先需要明白*p*值小于0.01的意义。从*p*值的计算过程中，我们可以得出*p*值小于0.01意味着该探针很大概率发出了荧光。但是，发出了荧光就意味着准确测量出甲基化水平了吗？大家注意观察上面的点图，红色的线是*p*值大于0的点的平均值。我们发现这些平均值都在0.5上下徘徊。有经验的同学就会想到，这些点的beta值没有饱和，导致在0.5上下徘徊。至此，我们认为*p*值小于0.01时，甲基化探针确实发出信号，但是并不意味着该信号达到饱和状态。下面就让我们一起来探索合理*p*值cutoff的选择。

### *p*值cutoff的选择
为了得到合理的cutoff，我们设计以下实验。首先，我们选取跨样本甲基化水平非常稳定的点，具体做法就是计算每个位点*p*值等于0的样本beta值的方差，然后选择方差较小（sd<0.05，此时我们假定该点的beta值在样本中呈方差为0.05的正态分布，如果beta值与均值的差距大于0.1288，则为离群点）的点。其次，我们计算当*p*值达到多少的时候，beta值会开始发生偏离。具体代码如下：
{% raw %}
```r
options(stringsAsFactors = F)

#读取p值和beta值文件，每个探针一行，每个样本一列
pvalue = read.table("p_value.txt", header = F, sep = "\t")
beta = read.table("beta_no_na.txt", header = F, sep = "\t")

#将探针号作为行名
row.names(beta) = beta[,1]
beta = beta[,-1]
row.names(pvalue) = pvalue[,1]
pvalue = pvalue[,-1]

#计算每一行中p值介于0和1e-6之间的个数
count_bad_row = function(p) {
  return(length(which(p>0 & p<1e-6)))
}
p_gt_0_row = apply(pvalue, 1, count_bad_row)

#计算每一行p值为0的样本beta值的方差
beta_var_row = c()
for (i in 1:dim(beta)[1]){
  beta_var_row=append(beta_var_row,
                      sd(as.numeric(beta[i,which(as.numeric(pvalue[i,]) == 0)])),
                      length(beta_var_row))
}

#选择方差较小且存在p值大于0的行作为测试
row_selected = which(p_gt_0_row > 1 & beta_var_row < 0.05)

#计算p值大于0时，beta值的离均差
pvalue_plot = c()
beta_dif_plot = c()
for (i in row_selected){
  beta_mean = mean(as.numeric(beta[i,which(pvalue[i,] == 0)]))
  p_gt_0 = which(pvalue[i,] > 0 & pvalue[i,]<1e-6)
  for (v in p_gt_0) {
    pvalue_plot = append(pvalue_plot,
                         pvalue[i,v],
                         length(pvalue_plot))
    beta_dif_plot = append(beta_dif_plot,
                           abs(beta[i,v] - beta_mean),
                           length(beta_dif_plot))
  }
}


#以p值为横坐标，beta值离均差为横坐标做散点图
pdf("p-value_vs_beta-dif.pdf")
plot(pvalue_plot, beta_dif_plot, xlim = c(0,1e-15),
     xlab = "p-value", ylab = "beta-dif")
abline(h=qnorm(p = 0.005, sd = 0.05, mean = 0,lower.tail = F),
       col="red")
dev.off()
```
{% endraw %}


这里要说明一下，用minfi计算*p*值时大于0的最小值为1.110223e-16（返回值类型为double.neg.eps，精度查看命令.Machine）。

采用上述方法，以beta值与均值的差值为纵坐标，*p*值为横坐标，我们得到下图。
![betadif](/img/in-post/2018-10-17-beta-dif.jpg)

从图中，可以看出当p值为1.110223e-16时，仍有很多点离群。我们再在另外几个数据集上做相同的测试。
