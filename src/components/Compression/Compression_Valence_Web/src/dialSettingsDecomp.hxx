#ifndef HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_DECOMP_H
#define HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_DECOMP_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence_Web

#include <QtGui/QDialog>

#include "ui_dialSettingsDecomp.h"

/**
 \class	SettingsDialogDecomp

 \brief	Settings dialog decomp. 

 */

class SettingsDialogDecomp : public QDialog, public Ui_SettingsDecomp
{
    Q_OBJECT

public:

    /**
     \fn	SettingsDialogDecomp::SettingsDialogDecomp(QWidget *parent = 0);
    
     \brief	Constructor.
    
    
     \param [in,out]	parent	If non-null, the parent.
     */

    SettingsDialogDecomp(QWidget *parent = 0);

    /**
     \fn	void SettingsDialogDecomp::accept();
    
     \brief	Accepts this object.
    
     */

    void accept();

private slots:

    /**
     \fn	void SettingsDialogDecomp::loadDefaults();
    
     \brief	Loads the defaults.

     */

    void loadDefaults();

    /**
     \fn	void SettingsDialogDecomp::loadFromSettings();
    
     \brief	Loads from settings.
    
     */

    void loadFromSettings();

    /**
     \fn	void SettingsDialogDecomp::saveToSettings();
    
     \brief	Saves to settings.

     */

    void saveToSettings();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_DECOMP_H