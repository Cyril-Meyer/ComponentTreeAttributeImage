TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        preprocess_nenist.cpp

HEADERS += \
    Algorithms/ComponentTree.h \
    Algorithms/ComponentTree.hxx \
    Algorithms/Morphology.h \
    Algorithms/Morphology.hxx \
    Common/FlatSE.h \
    Common/FlatSE.hxx \
    Common/Image.h \
    Common/Image.hxx \
    Common/ImageIO.hxx \
    Common/ImageIterators.h \
    Common/Point.h \
    Common/Types.h
