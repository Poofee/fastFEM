# fastFEM
这个程序主要实现有限元的求解，目前用于静磁场的求解。

1. ## 开发环境

   - visual studio 2013 官网上可以下载到免费的社区版
   - qt 需要去http://download.qt.io/ 下载vsaddin 和 qt库，可以选择http://download.qt.io/official_releases/vsaddin/qt-vs-addin-1.2.4-opensource.exe 和 http://download.qt.io/development_releases/qt/5.5/5.5.0-rc/qt-opensource-windows-x86-msvc2013_64-5.5.0-rc.exe 如果qt版本过高，建议安装最新的qtvs插件。

   ## 2.注意事项

   - [ ] 新建工程的时候，如果需要Qt，尽量先用Qtcreator建立一个空的工程，然后导入到visual studio，这样会在工程当中用的是QTDIR这个宏，而不是全路径：
test

test2

test3last