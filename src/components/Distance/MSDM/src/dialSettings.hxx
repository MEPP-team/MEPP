///////////////////////////////////////////////////////////////////////////
// Author: Guillaume Lavoué
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef HEADER_MEPP_COMPONENT_MSDM_PLUGIN_SETTINGS_H
#define HEADER_MEPP_COMPONENT_MSDM_PLUGIN_SETTINGS_H

#include <mepp_config.h>
#ifdef BUILD_component_MSDM

#ifndef _MSC_VER
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
#include <QtGui>
#ifndef _MSC_VER
#pragma GCC diagnostic warning "-Wuninitialized"
#endif

#if QT_VERSION >= 0x050000
#include <QtWidgets>
#endif

#include "ui_dialSettings.h"

class SettingsDialog : public QDialog, public Ui_Settings
{
    Q_OBJECT

public:
    SettingsDialog(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_MSDM_PLUGIN_SETTINGS_H