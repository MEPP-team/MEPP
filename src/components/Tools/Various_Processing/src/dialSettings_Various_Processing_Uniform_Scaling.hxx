#ifndef HEADER_MEPP_COMPONENT_Various_Processing_Uniform_Scaling_PLUGIN_SETTINGS_H
#define HEADER_MEPP_COMPONENT_Various_Processing_Uniform_Scaling_PLUGIN_SETTINGS_H

#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include <QtGui/QDialog>

#include "ui_dialSettings_Various_Processing_Uniform_Scaling.h"

class SettingsDialog_Various_Processing_Uniform_Scaling : public QDialog, public Ui_Settings_Uniform_Scaling
{
    Q_OBJECT

public:
    SettingsDialog_Various_Processing_Uniform_Scaling(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_Various_Processing_Uniform_Scaling_PLUGIN_SETTINGS_H