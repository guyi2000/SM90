# SM90 结构力学求解器

## 说明

本程序在结构力学求解器教学版及基础上编制了核心计算程序，完成了超静定结构的内力计算以及自由振动的频率计算功能。程序中使用变列宽分解法求解线性方程组，荷载由节点荷载、与单元荷载，其中单元荷载仅有 `4` 类，以及除了斜连接外的结点连接及除斜支座、弹性支座之外的支座约束条件。自由振动采用 `Wittrick-Williams` 算法，方法简单方便，可以达到任意精度并且不会丢根、漏根，且单根、重根均可处理。

使用 `NumKind` 模块进行数据精度的定义，`TypeDef` 模块定义超静定结构求解使用的各种自定义类型，同时实现了单元属性求解及变换矩阵求解的子程序。`BandMat` 模块实现了变带宽矩阵的存贮、释放与求解。`DispMethod` 模块为主要编写的模块，其中实现了单元刚度矩阵求解、整体刚度矩阵集成、固端力求解，整体荷载向量求解，以及之后的单元位移与单元力的求解。`FreqMethod` 为主要频率求解模块，其中实现了固端频率计算、频率数计算，以及单元动力刚度阵计算以及整体动力刚度矩阵集成。
