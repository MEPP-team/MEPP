#ifndef HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H
#define HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence_Web

#include <QtGui/QDialog>

#include "ui_dialSettingsComp.h"

/**
 \class	SettingsDialogComp

 \brief	Settings dialog comp. 

 */

class SettingsDialogComp : public QDialog, public Ui_SettingsComp
{
    Q_OBJECT

public:

    /**
     \fn	SettingsDialogComp::SettingsDialogComp(QWidget *parent = 0);
    
     \brief	Constructor.
    
    
     \param [in,out]	parent	If non-null, the parent.
     */

    SettingsDialogComp(QWidget *parent = 0);

    /**
     \fn	void SettingsDialogComp::accept();
    
     \brief	Accepts this object.
    
     */

    void accept();

private slots:

    /**
     \fn	void SettingsDialogComp::loadDefaults();
    
     \brief	Loads the defaults.
    
     */

    void loadDefaults();

    /**
     \fn	void SettingsDialogComp::loadFromSettings();
    
     \brief	Loads from settings.
    
     */

    void loadFromSettings();

    /**
     \fn	void SettingsDialogComp::saveToSettings();
    
     \brief	Saves to settings.
    
     */

    void saveToSettings();

	/**
	 \fn	void SettingsDialogComp::setFilename();
	
	 \brief	Sets the filename.

	 */

	void setFilename();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H