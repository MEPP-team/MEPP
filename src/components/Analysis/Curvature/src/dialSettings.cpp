#include <mepp_config.h>
#ifdef BUILD_component_Curvature

#include "dialSettings.hxx"

SettingsDialog::SettingsDialog(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog::loadDefaults()
{
}

void SettingsDialog::loadFromSettings()
{
}

void SettingsDialog::saveToSettings()
{
}

void SettingsDialog::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif