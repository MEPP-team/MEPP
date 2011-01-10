#ifndef HEADER_MEPP_COMPONENT_Various_Processing_Simplification_PLUGIN_SETTINGS_H
#define HEADER_MEPP_COMPONENT_Various_Processing_Simplification_PLUGIN_SETTINGS_H

#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include <QtGui/QDialog>

#include "ui_dialSettings_Various_Processing_Simplification.h"

class SettingsDialog_Various_Processing_Simplification : public QDialog, public Ui_Settings_Simplification
{
    Q_OBJECT

public:
    SettingsDialog_Various_Processing_Simplification(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_Various_Processing_Simplification_PLUGIN_SETTINGS_H