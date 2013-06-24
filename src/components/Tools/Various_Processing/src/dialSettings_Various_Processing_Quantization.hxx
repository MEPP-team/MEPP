#ifndef HEADER_MEPP_COMPONENT_Various_Processing_Quantization_PLUGIN_SETTINGS_H
#define HEADER_MEPP_COMPONENT_Various_Processing_Quantization_PLUGIN_SETTINGS_H

#include <mepp_config.h>
#ifdef BUILD_component_Various_Processing

#ifndef _MSC_VER
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
#include <QtGui>
#ifndef _MSC_VER
#pragma GCC diagnostic warning "-Wuninitialized"
#endif

#if QT_VERSION >= 0x050000
#include <QtWidgets>
#endif

#include "ui_dialSettings_Various_Processing_Quantization.h"

class SettingsDialog_Various_Processing_Quantization : public QDialog, public Ui_Settings_Quantization
{
    Q_OBJECT

public:
    SettingsDialog_Various_Processing_Quantization(QWidget *parent = 0);
    void accept();

private slots:
    void loadDefaults();
    void loadFromSettings();
    void saveToSettings();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_Various_Processing_Quantization_PLUGIN_SETTINGS_H