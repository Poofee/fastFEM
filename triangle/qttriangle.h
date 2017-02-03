#ifndef QTTRIANGLE_H
#define QTTRIANGLE_H

#include <QMainWindow>

namespace Ui {
class Qttriangle;
}

class Qttriangle : public QMainWindow
{
    Q_OBJECT

public:
    explicit Qttriangle(QWidget *parent = 0);
    ~Qttriangle();

private:
    Ui::Qttriangle *ui;
};

#endif // QTTRIANGLE_H
