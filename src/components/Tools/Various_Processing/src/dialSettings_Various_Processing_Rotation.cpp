#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "dialSettings_Various_Processing_Rotation.hxx"

SettingsDialog_Various_Processing_Rotation::SettingsDialog_Various_Processing_Rotation(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog_Various_Processing_Rotation::loadDefaults()
{
}

void SettingsDialog_Various_Processing_Rotation::loadFromSettings()
{
}

void SettingsDialog_Various_Processing_Rotation::saveToSettings()
{
}

void SettingsDialog_Various_Processing_Rotation::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif