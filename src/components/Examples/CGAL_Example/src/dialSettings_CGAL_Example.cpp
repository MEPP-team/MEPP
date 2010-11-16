///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#include <mepp_config.h>
#ifdef BUILD_component_CGAL_Example

#include "dialSettings_CGAL_Example.hxx"

SettingsDialog_CGAL_Example::SettingsDialog_CGAL_Example(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    loadDefaults();
    loadFromSettings();
}

void SettingsDialog_CGAL_Example::loadDefaults()
{
}

void SettingsDialog_CGAL_Example::loadFromSettings()
{
}

void SettingsDialog_CGAL_Example::saveToSettings()
{
}

void SettingsDialog_CGAL_Example::accept()
{
    saveToSettings();

    QDialog::accept();
}

#endif