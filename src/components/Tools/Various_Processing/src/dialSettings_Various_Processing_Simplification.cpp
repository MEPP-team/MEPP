#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "dialSettings_Various_Processing_Simplification.hxx"

SettingsDialog_Various_Processing_Simplification::SettingsDialog_Various_Processing_Simplification(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog_Various_Processing_Simplification::loadDefaults()
{
}

void SettingsDialog_Various_Processing_Simplification::loadFromSettings()
{
}

void SettingsDialog_Various_Processing_Simplification::saveToSettings()
{
}

void SettingsDialog_Various_Processing_Simplification::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif