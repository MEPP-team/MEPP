#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#include "dialSettings_Various_Processing_Quantization.hxx"

SettingsDialog_Various_Processing_Quantization::SettingsDialog_Various_Processing_Quantization(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog_Various_Processing_Quantization::loadDefaults()
{
}

void SettingsDialog_Various_Processing_Quantization::loadFromSettings()
{
}

void SettingsDialog_Various_Processing_Quantization::saveToSettings()
{
}

void SettingsDialog_Various_Processing_Quantization::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif