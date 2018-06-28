#-------------------------------------------------
#
# Project created by QtCreator 2017-03-27T09:40:41
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport#necessary for qcustomplot

#output excution file name
TARGET = fastFEM

#project type:app lib vcapp vclib
TEMPLATE = app

QMAKE_CXXFLAGS += /openmp

CONFIG += debug_and_release

#source file lists
SOURCES += \
    datatype.cpp \
    FastFEMcore.cpp \
    main.cpp \
    material.cpp \
    plot.cpp \
    SuperLU_MT.cpp \
    ../qcustomplot/qcustomplot.cpp \
    mesh/meshtype.cpp


#.ui file lists
FORMS    += mainwindow.ui

#header file lists
HEADERS += \
    datatype.h \
    FastFEMcore.h \
    plot.h \
    spline.h \
    SuperLU_MT.h \
    ../qcustomplot/qcustomplot.h \
    mesh/meshtype.h


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
	..\OpenBLAS-v0.2.15-Win64-int32\include \

}
