#ifndef HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H
#define HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include <QtGui/QDialog>

#include "ui_dialSettingsComp.h"

class SettingsDialogComp : public QDialog, public Ui_SettingsComp
{
    Q_OBJECT

public:
    SettingsDialogComp(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

	void setFilename();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H