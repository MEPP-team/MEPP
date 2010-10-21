///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#include <mepp_config.h>
#ifdef BUILD_component_CGAL_Example

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