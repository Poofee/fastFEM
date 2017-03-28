#-------------------------------------------------
#
# Project created by QtCreator 2017-03-27T09:40:41
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

#output excution file name
TARGET = fastFEM

#project type:app lib vcapp vclib
TEMPLATE = app

#output directory
DESTDIR = ../x64/Linux

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
    SuperLU_MT.h


unix {
#librarys
LIBS +=

#additional include path
INCLUDEPATH += \
    ../armadillo/include \
    ../qcustomplot \
    ../SuperLU_MT_3.1/SRC \
}

win32 {
#librarys
LIBS +=

}
