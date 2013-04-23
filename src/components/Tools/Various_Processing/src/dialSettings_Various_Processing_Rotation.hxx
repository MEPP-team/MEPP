#ifndef HEADER_MEPP_COMPONENT_Various_Processing_Rotation_PLUGIN_SETTINGS_H
#define HEADER_MEPP_COMPONENT_Various_Processing_Rotation_PLUGIN_SETTINGS_H

#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#ifndef _MSC_VER
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
#include <QtGui/QDialog>
#ifndef _MSC_VER
#pragma GCC diagnostic warning "-Wuninitialized"
#endif

#include "ui_dialSettings_Various_Processing_Rotation.h"

class SettingsDialog_Various_Processing_Rotation : public QDialog, public Ui_Settings_Rotation
{
    Q_OBJECT

public:
    SettingsDialog_Various_Processing_Rotation(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_Various_Processing_Rotation_PLUGIN_SETTINGS_H