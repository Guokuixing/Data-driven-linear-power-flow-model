## Data-driven-linear-power-flow-model

本算法为毕业设计《数据驱动的线性化潮流方程建模方法及开源实现》的研究内容，作者郭魁星，指导教师康重庆教授、张宁副教授，单位为清华大学电机系。

使用Caculation.m文件可以实现对IEEE任意节点系统进行潮流分布计算，输入case及节点电压幅值和相角数据即可获得系统潮流分布。

StartUp.m文件为本研究提出的线性化潮流方程建模并应用于N-1准则检验的算法实现过程。

经本研究验证，该算法在IEEE-5、30、33bw、57、118节点系统中均实现了较高的计算精度和较快的计算速度。

开源此代码希望对您的研究和个人使用有所帮助。任何商业用途都应得到上述论文的作者的同意。

此代码需要 Matpower平台和 Matlab python引擎。