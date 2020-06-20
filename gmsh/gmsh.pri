DEFINES += _USE_MATH_DEFINES NOMINMAX WIN32

win32-msvc* {
    #ignore warning C4819
    QMAKE_CXXFLAGS += /wd"4819"
    QMAKE_CXXFLAGS *=  /wd"4100"
    contains (QMAKE_CXXFLAGS_WARN_ON, -w34100) : QMAKE_CXXFLAGS_WARN_ON -= -w34100
    DEFINES += _CRT_SECURE_NO_WARNINGS
}

include($$PWD/Geo/Geo.pri)
include($$PWD/Common/Common.pri)
include($$PWD/Mesh/Mesh.pri)
include($$PWD/Numeric/Numeric.pri)
include($$PWD/Parser/Parser.pri)
include($$PWD/Post/Post.pri)
