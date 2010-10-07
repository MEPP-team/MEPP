#ifndef HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_DECOMP_H
#define HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_DECOMP_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include <QtGui/QDialog>

#include "ui_dialSettingsDecomp.h"

class SettingsDialogDecomp : public QDialog, public Ui_SettingsDecomp
{
    Q_OBJECT

public:
    SettingsDialogDecomp(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_DECOMP_H