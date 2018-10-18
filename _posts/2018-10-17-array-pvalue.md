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

## P值如何得到

想确定芯片处理时*p*值的cutoff,首先需要知道*p*值是如何得到。

下面是minfi包中关于detectionP函数的源码。

{% highlight r linenos %}

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

{% endhighlight %}
