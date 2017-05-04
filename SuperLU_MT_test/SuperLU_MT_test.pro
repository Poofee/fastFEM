#-------------------------------------------------
#
# Project created by QtCreator 2017-03-27T09:40:41
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport#necessary for qcustomplot

#output excution file name
TARGET = SuperLU_MT_test

#project type:app lib vcapp vclib
TEMPLATE = app

CONFIG += debug_and_release

#source file lists
SOURCES += \
    ../qcustomplot/qcustomplot.cpp \
    main.cpp \
    mainwindow.cpp


#.ui file lists
FORMS    += mainwindow.ui

#header file lists
HEADERS += \
    ../qcustomplot/qcustomplot.h \
    mainwindow.h


unix {

CONFIG(debug,debug|release){
#output directory
DESTDIR = ../x64/Linux/debug
}else{
#output directory
DESTDIR = ../x64/Linux/release
}
#librarys
LIBS += \
    -L../SuperLU_MT_3.1 -lsuperlu_mt_OPENMP \
    -lgomp\
    -lpthread \
    ../armadillo//lib/libarmadillo.so \
    ../openblas/lib/libopenblas.so

#additional include path
INCLUDEPATH += \
    ../armadillo/include \
    ../qcustomplot \
    ../SuperLU_MT_3.1/SRC \
}

win32 {
CONFIG(debug,debug|release){
#output directory
DESTDIR = ../x64/Debug
}else{
#output directory
DESTDIR = ../x64/Release
}
#librarys
LIBS += \
    ..\SuperLU_MT_3.1\SuperLU_MT_3.1.lib \
    ..\OpenBLAS-v0.2.15-Win64-int32\lib\libopenblas.dll.a

#additional include path
INCLUDEPATH += \
    ..\armadillo\include \
    ..\qcustomplot \
    ..\SuperLU_MT_3.1\SRC \

}
