---
title: 'Stata graph edit'
permalink: /posts/2025/04/blog-post-1/
date: 2025-04-11
tags:
  - stata
---

微信公众号推文链接： [Stata修图转代码](https://mp.weixin.qq.com/s/2E8R-I_6FLn7uTsh5L9Nag)

## 适用场景
今天用stata跑代码写论文遇到了一个问题：我在进行a安慰剂检验的时候，使用了陈强老师（2025）《管理世界》上开发的新命令didplacebo，代码大概是这样的：

```
foreach i of varlist $main_dependent{
reghdfe `i' did $control, absorb($FE) $se_vce
est store `i'
didplacebo `i', treatvar(did) pbotime(1(1)1000)
graph save $main/result/Placebo_in-time_`i', replace
graph export $main/result/Placebo_in-time_`i'.png, as(png) name("pbotime") replace
}
*ssc install is needed here.
grc1leg2 $main/result/Placebo_in-time_001.gph $main/result/Placebo_in-time_002.gph , col(2) iscale(1) xsize(100) ysize(50) labsize(medsmall) lrows(1) mainto
```

我希望修改两个子图的标题为变量标签。然而，didplacebo似乎无法DIY subtitle，grc1leg2也没法解决这个问题。只能通过图形编辑器来手动DIY了。

## 怎么办？

由于对于stata来说，我们的手动DIY也相当于是给它派命令，所以手动DIY并不可怕，只要**【do文件可复现】** 就好了！

通过bing搜索，我找到了这篇文章~
[Stata修图与操作记录](https://czxb.github.io/ar/Stata修图与操作记录.html)

解决方案如下：
- 首先运行
```
grc1leg2 $main/result/Placebo_in-time_001.gph $main/result/Placebo_in-time_002.gph , col(2) iscale(1) xsize(100) ysize(50) labsize(medsmall) lrows(1) mainto
```

- 点击图形编辑器，点击那个看起来就很像录制的红色图标。

- 在可视化界面中DIY执行操作

- 再次点击那个看起来就很像录制的红色图标。把录制文件保存到桌面（可以起个名字方便后续好找）

- 用vscode打开那个录制文件，大概是这样的

<img src="https://files.mdnice.com/user/77043/1e5de925-0e78-4365-9de1-f8238198d75c.png" style="zoom:50%;" />

- 因此，只需要在现有代码后面加入新行
```
gr_edit .plotregion1.graph1.subtitle.text.Arrpush a. 0001

gr_edit .plotregion1.graph2.subtitle.text.Arrpush b. 0002
```

这就ok了。现在，不需要我们手动DIY，stata就会执行我们DIY的操作，这个过程可复现，放进do file对应位置就ok啦~

**参考链接**
https://czxb.github.io/ar/Stata修图与操作记录.html

