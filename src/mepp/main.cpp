///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#include <QApplication>
#include <QtGui>

#include "mainwindow.hxx"

void loadStyleSheet()
{
    // Let's use QFile and point to a resource...
    QFile data("./mepp.qss");
    QString style;

    // ...to open the file
    if (data.open(QFile::ReadOnly))
	{
        // QTextStream...
        QTextStream styleIn(&data);
        // ...read file to a string.
        style = styleIn.readAll();
        data.close();

        // We'll use qApp macro to get the QApplication pointer and set the style sheet application wide.
        qApp->setStyleSheet(style);
    }
}

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    mainwindow window;
	loadStyleSheet();
	window.show();

    return app.exec();
}