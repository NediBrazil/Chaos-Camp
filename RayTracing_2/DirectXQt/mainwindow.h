#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#pragma once
#include <QMainWindow>
#include <QTimer>
#include <QLabel>
#include "DXRenderer.h"
#include "ViewportWidget.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    MainWindow(QWidget* parent = nullptr);

private:
    DXRenderer renderer;
    ViewportWidget* viewport;
    QLabel* fpsLabel;
    bool initialized = false;

    QTimer renderTimer;
    QTimer fpsTimer;

    int frameCounter = 0;

    void renderToQt();

protected:
    void keyPressEvent(QKeyEvent* event) override;
};

#endif // MAINWINDOW_H
