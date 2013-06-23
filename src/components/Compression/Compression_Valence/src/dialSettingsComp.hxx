#ifndef HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H
#define HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#ifndef _MSC_VER
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
#include <QtGui>
#ifndef _MSC_VER
#pragma GCC diagnostic warning "-Wuninitialized"
#endif

#if (QT_VERSION >= QT_VERSION_CHECK(5, 0, 0))
#include <QtWidgets>
#endif

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

    SettingsDialogComp(QWidget *parent, QString &saveLocation);

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
	QString &saveLocation_;
};

#endif

#endif // HEADER_MEPP_COMPONENT_VALENCE_PLUGIN_SETTINGS_COMP_H