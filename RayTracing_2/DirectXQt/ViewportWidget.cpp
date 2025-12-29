#include "ViewportWidget.h"
#include <QPainter>

ViewportWidget::ViewportWidget(QWidget* parent)
    : QWidget(parent)
{
    setAttribute(Qt::WA_NativeWindow);
}

void ViewportWidget::updateFrame(const QImage& img)
{
    frame = img;
    update();
}

void ViewportWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.setRenderHint(QPainter::SmoothPixmapTransform, false);
    p.drawImage(rect(), frame, frame.rect());
}

