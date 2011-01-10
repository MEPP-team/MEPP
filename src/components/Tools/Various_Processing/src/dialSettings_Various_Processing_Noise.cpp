#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "dialSettings_Various_Processing_Noise.hxx"

SettingsDialog_Various_Processing_Noise::SettingsDialog_Various_Processing_Noise(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog_Various_Processing_Noise::loadDefaults()
{
}

void SettingsDialog_Various_Processing_Noise::loadFromSettings()
{
}

void SettingsDialog_Various_Processing_Noise::saveToSettings()
{
}

void SettingsDialog_Various_Processing_Noise::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif