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

CONFIG += debug_and_release

#source file lists
SOURCES += \
    datatype.cpp \
    FastFEMcore.cpp \
    main.cpp \
    material.cpp \
    plot.cpp \
    SuperLU_MT.cpp \
    ../qcustomplot/qcustomplot.cpp


#.ui file lists
FORMS    += mainwindow.ui

#header file lists
HEADERS += \
    datatype.h \
    FastFEMcore.h \
    plot.h \
    spline.h \
    SuperLU_MT.h \
    ../qcustomplot/qcustomplot.h


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
    ../armadillo//lib/libarmadillo.so \
    -L../spral/lib -lspralgpu \
    -L../metis -lmetis \
    -L../openblas/lib -lopenblas \
    -L/usr/local/cuda/lib64 -lcudart \
    -L/usr/local/cuda/lib64  -lcublas\
    -lgfortran -lm -lquadmath  \
    -lpthread \
    -lgomp

#additional include path
INCLUDEPATH += \
    ../armadillo/include \
    ../qcustomplot \
    ../SuperLU_MT_3.1/SRC \
    ../spral/include \
}

win32 {
#librarys
LIBS +=

}
