#ifndef HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_JCW_H
#define HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_JCW_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include <QtGui/QDialog>

#include "ui_dialSettingsJCW.h"

/**
 \class	SettingsDialogJCW

 \brief	Settings dialog JCW. 

 */

class SettingsDialogJCW : public QDialog, public Ui_SettingsJCW
{
    Q_OBJECT

public:

    /**
     \fn	SettingsDialogJCW::SettingsDialogJCW(QWidget *parent = 0);
    
     \brief	Constructor.
    
    
     \param [in,out]	parent	If non-null, the parent.
     */

    SettingsDialogJCW(QWidget *parent = 0);

    /**
     \fn	void SettingsDialogJCW::accept();
    
     \brief	Accepts this object.
    
     */

    void accept();

private slots:

    /**
     \fn	void SettingsDialogJCW::loadDefaults();
    
     \brief	Loads the defaults.
    
     */

    void loadDefaults();

    /**
     \fn	void SettingsDialogJCW::loadFromSettings();
    
     \brief	Loads from settings.
    
     */

    void loadFromSettings();

    /**
     \fn	void SettingsDialogJCW::saveToSettings();
    
     \brief	Saves to settings.
    
     */

    void saveToSettings();

	/**
	 \fn	void SettingsDialogJCW::setFilename();
	
	 \brief	Sets the filename.

	 */

	void setFilename();

private:
};

#endif

#endif // HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H