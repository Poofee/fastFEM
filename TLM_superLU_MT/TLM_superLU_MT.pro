#-------------------------------------------------
#
# Project created by QtCreator 2016-04-11T14:20:28
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets  printsupport

TARGET = TLM_superLU_MT
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow.cpp \
    datatype.cpp\
    qcustomplot.cpp

HEADERS  += mainwindow.h \
    datatype.h\
    qcustomplot.h

FORMS    += mainwindow.ui

INCLUDEPATH += .\
    E:\Projects\cplusplus\SuperLU_MT_3.1\SRC \
    E:\Projects\cplusplus\OpenBLAS-v0.2.15-Win64-int32\include \
    E:\Projects\cplusplus\armadillo-6.500.4\include

LIBS += E:\Projects\cplusplus\SuperLU_MT_3.1\SuperLU_MT\x64\Release \
    E:\Projects\cplusplus\OpenBLAS-v0.2.15-Win64-int32\lib \
