#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "dialSettings_Various_Processing_Subdivision.hxx"

SettingsDialog_Various_Processing_Subdivision::SettingsDialog_Various_Processing_Subdivision(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog_Various_Processing_Subdivision::loadDefaults()
{
}

void SettingsDialog_Various_Processing_Subdivision::loadFromSettings()
{
}

void SettingsDialog_Various_Processing_Subdivision::saveToSettings()
{
}

void SettingsDialog_Various_Processing_Subdivision::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif