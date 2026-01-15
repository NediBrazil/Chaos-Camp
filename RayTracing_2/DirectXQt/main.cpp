#include <QApplication>
#include "MainWindow.h"

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    MainWindow w;
    w.showFullScreen();
    return app.exec();
}
