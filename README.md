# fastFEM
这个程序主要实现有限元的求解，目前用于静磁场的求解。

## 1.开发环境

- visual studio 2013 官网上可以下载到免费的社区版
- qt 需要去http://download.qt.io/ 下载vsaddin 和 qt库，可以选择http://download.qt.io/official_releases/vsaddin/qt-vs-addin-1.2.4-opensource.exe 和 http://download.qt.io/development_releases/qt/5.5/5.5.0-rc/qt-opensource-windows-x86-msvc2013_64-5.5.0-rc.exe 如果qt版本过高，建议安装最新的qtvs插件。

## 2.注意事项

- [ ] 新建工程的时候，如果需要Qt，尽量先用Qtcreator建立一个空的工程，然后导入到visual studio，这样会在工程当中用的是QTDIR这个宏，而不是全路径：
- [ ] VS的.user文件，如果不存在的话，每一次新打开都会自动生成一个，所以理论上不用将它添加追踪，但是似乎是有什么bug，导致生成的文件QTDIR的定义在使用之后，使得它的定义就失去了意义，系统就找不到路径了，所以还是加上吧；

## 3.目录说明

├── armadillo	**<u>armadillo库文件</u>**
├── fastFEM		**<u>主工程</u>**
├── LICENSE		**<u>LICENSE</u>**
├── matlab		**<u>MATLAB测试程序</u>**
├── model		**<u>分网文件</u>**
├── openblas	**<u>openblas在Linux下库文件</u>**
├── OpenBLAS-v0.2.15-Win64-int32	**<u>win64下openblas库文件</u>**
├── qcustomplot		**<u>qcustomplot库文件</u>**
├── README.md		**<u>本文档</u>**
├── spral			**<u>spral库文件</u>**
├── SuperLU_MT_3.1		**<u>superluMT源代码</u>**
├── SuperLU_MT_test	**<u>superLUMT测试项目</u>**
├── tree.txt				**<u>文件树</u>**
├── triangle				**<u>triangle库文件</u>**
└── x64					**<u>生成目录</u>**