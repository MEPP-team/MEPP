#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "dialSettings_Various_Processing_Smoothing.hxx"

SettingsDialog_Various_Processing_Smoothing::SettingsDialog_Various_Processing_Smoothing(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog_Various_Processing_Smoothing::loadDefaults()
{
}

void SettingsDialog_Various_Processing_Smoothing::loadFromSettings()
{
}

void SettingsDialog_Various_Processing_Smoothing::saveToSettings()
{
}

void SettingsDialog_Various_Processing_Smoothing::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif