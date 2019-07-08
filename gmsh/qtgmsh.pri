
win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../bin/ -lgmsh
     else:LIBS += -L$$PWD/../bin/ -lgmsh

INCLUDEPATH += $$PWD/
DEPENDPATH += $$PWD/


