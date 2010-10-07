#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include <QDir>
#include <QFileDialog>
#include <QMessageBox>

#include "dialSettingsComp.hxx"

SettingsDialogComp::SettingsDialogComp(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

	connect(pushButtonBrowse, SIGNAL(clicked()), this, SLOT(setFilename()));

    loadDefaults();
    loadFromSettings();
}

void SettingsDialogComp::loadDefaults()
{
}

void SettingsDialogComp::loadFromSettings()
{
}

void SettingsDialogComp::saveToSettings()
{
}

void SettingsDialogComp::accept()
{
    saveToSettings();

	if (!file_name->text().isEmpty())
		QDialog::accept();
	else
		QMessageBox::information(this, APPLICATION, tr("You must specify a filename."));
}

void SettingsDialogComp::setFilename()
{
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save P3D File - from Valence"),
											 QDir::currentPath(), //QString()
											 tr("P3D files (*.p3d)"));
	if (!fileName.isEmpty())
		file_name->setText(fileName);
}

#endif