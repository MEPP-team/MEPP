#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "dialSettings_Various_Processing_Translation.hxx"

SettingsDialog_Various_Processing_Translation::SettingsDialog_Various_Processing_Translation(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog_Various_Processing_Translation::loadDefaults()
{
}

void SettingsDialog_Various_Processing_Translation::loadFromSettings()
{
}

void SettingsDialog_Various_Processing_Translation::saveToSettings()
{
}

void SettingsDialog_Various_Processing_Translation::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif