#ifndef VIEWPORTWIDGET_H
#define VIEWPORTWIDGET_H
#pragma once
#include <QWidget>
#include <QImage>

class ViewportWidget : public QWidget
{
    Q_OBJECT

public:
    ViewportWidget(QWidget* parent = nullptr);
    void updateFrame(const QImage& img);

protected:
    void paintEvent(QPaintEvent* event) override;

private:
    QImage frame;
};

#endif
