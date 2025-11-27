#include "ViewportWidget.h"
#include <QPainter>

ViewportWidget::ViewportWidget(QWidget* parent)
    : QWidget(parent)
{
}

void ViewportWidget::updateFrame(const QImage& img)
{
    frame = img;
    update();
}

void ViewportWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.drawImage(rect(), frame);
}
