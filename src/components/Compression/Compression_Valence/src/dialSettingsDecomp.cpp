#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "dialSettingsDecomp.hxx"

SettingsDialogDecomp::SettingsDialogDecomp(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialogDecomp::loadDefaults()
{
}

void SettingsDialogDecomp::loadFromSettings()
{
}

void SettingsDialogDecomp::saveToSettings()
{
}

void SettingsDialogDecomp::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif