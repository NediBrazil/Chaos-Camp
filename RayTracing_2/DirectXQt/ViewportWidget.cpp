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

    QSize targetSize = frame.size();
    targetSize.scale(size(), Qt::KeepAspectRatio);

    QRect targetRect(
        QPoint(0, 0),
        targetSize
        );
    targetRect.moveCenter(rect().center());

    p.drawImage(targetRect, frame);
}
