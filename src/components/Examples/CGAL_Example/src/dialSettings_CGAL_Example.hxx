///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef HEADER_MEPP_COMPONENT_CGAL_EXAMPLE_PLUGIN_SETTINGS_H
#define HEADER_MEPP_COMPONENT_CGAL_EXAMPLE_PLUGIN_SETTINGS_H

#include <mepp_config.h>
#ifdef BUILD_component_CGAL_Example

#ifndef _MSC_VER
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
#include <QtGui>
#ifndef _MSC_VER
#pragma GCC diagnostic warning "-Wuninitialized"
#endif

#if (QT_VERSION >= QT_VERSION_CHECK(5, 0, 0))
#include <QtWidgets>
#endif

#include "ui_dialSettings_CGAL_Example.h"

class SettingsDialog_CGAL_Example : public QDialog, public Ui_Settings
{
    Q_OBJECT

public:
    SettingsDialog_CGAL_Example(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_CGAL_EXAMPLE_PLUGIN_SETTINGS_H