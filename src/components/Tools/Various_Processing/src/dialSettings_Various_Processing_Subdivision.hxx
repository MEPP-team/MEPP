#ifndef HEADER_MEPP_COMPONENT_Various_Processing_Subdivision_PLUGIN_SETTINGS_H
#define HEADER_MEPP_COMPONENT_Various_Processing_Subdivision_PLUGIN_SETTINGS_H

#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

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

#include "ui_dialSettings_Various_Processing_Subdivision.h"

class SettingsDialog_Various_Processing_Subdivision : public QDialog, public Ui_Settings_Subdivision
{
    Q_OBJECT

public:
    SettingsDialog_Various_Processing_Subdivision(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_Various_Processing_Subdivision_PLUGIN_SETTINGS_H