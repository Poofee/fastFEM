#include "qttriangle.h"
#include "ui_qttriangle.h"

Qttriangle::Qttriangle(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Qttriangle)
{
    ui->setupUi(this);
}

Qttriangle::~Qttriangle()
{
    delete ui;
}
