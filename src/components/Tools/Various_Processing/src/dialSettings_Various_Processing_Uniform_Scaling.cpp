#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "dialSettings_Various_Processing_Uniform_Scaling.hxx"

SettingsDialog_Various_Processing_Uniform_Scaling::SettingsDialog_Various_Processing_Uniform_Scaling(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog_Various_Processing_Uniform_Scaling::loadDefaults()
{
}

void SettingsDialog_Various_Processing_Uniform_Scaling::loadFromSettings()
{
}

void SettingsDialog_Various_Processing_Uniform_Scaling::saveToSettings()
{
}

void SettingsDialog_Various_Processing_Uniform_Scaling::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif