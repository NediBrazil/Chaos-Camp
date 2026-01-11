#include "mainwindow.h"
#include "MainWindow.h"
#include <QVBoxLayout>
#include <QKeyEvent>

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    resize(1000, 600);
    QWidget* c = new QWidget;
    setCentralWidget(c);

    viewport = new ViewportWidget;
    setFocusPolicy(Qt::StrongFocus);
    viewport->setFocusPolicy(Qt::StrongFocus);

    viewport->setAttribute(Qt::WA_NativeWindow);
    viewport->winId();

    fpsLabel = new QLabel("FPS: 0");

    QVBoxLayout* layout = new QVBoxLayout(c);
    layout->addWidget(viewport);
    layout->addWidget(fpsLabel);

    connect(&renderTimer, &QTimer::timeout, this, &MainWindow::renderToQt);
    renderTimer.start(0);

    connect(&fpsTimer, &QTimer::timeout, [&] {
        fpsLabel->setText(QString("FPS: %1").arg(frameCounter));
        frameCounter = 0;
    });
    fpsTimer.start(1000);
}

void MainWindow::keyPressEvent(QKeyEvent* event)
{
    if (event->key() == Qt::Key_R)
        renderer.setRenderMode(RenderMode::Raster);

    if (event->key() == Qt::Key_T)
        renderer.setRenderMode(RenderMode::RayTracing);
}

void MainWindow::renderToQt()
{
    if (!initialized)
    {
        HWND hwnd = (HWND)viewport->winId();
        if (IsWindow(hwnd))
        {
            initialized = renderer.initialize(
                hwnd,
                viewport->width(),
                viewport->height()
                );
        }
        return;
    }
    renderer.renderFrame();

    unsigned char* p = renderer.getBackBufferCPU();
    if (!p) return;

    QImage img(
        p,
        renderer.getWidth(),
        renderer.getHeight(),
        renderer.getRowPitch(),
        QImage::Format_RGBA8888);

    viewport->updateFrame(img);
    frameCounter++;
}
