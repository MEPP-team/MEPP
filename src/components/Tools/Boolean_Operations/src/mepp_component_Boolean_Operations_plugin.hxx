#ifndef HEADER_MEPP_COMPONENT_BOOLEAN_OPERATIONS_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_BOOLEAN_OPERATIONS_PLUGIN_INTERFACE_H

/*!
 * \file mepp_component_Boolean_Operations_plugin.hxx
 */

#include <QtGlobal> // important, for QT_VERSION

#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_Boolean_Operations

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

class mepp_component_Boolean_Operations_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);
#if QT_VERSION >= 0x050000
	Q_PLUGIN_METADATA(IID "mepp_component_Boolean_Operations_plugin")
#endif

public:
	void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
	{
		this->mw = mainWindow;
		this->lwindow = lw;
		this->mPluginName = this->metaObject()->className();
		
		// choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
		// menuCompression, menuWatermaking, menuExamples
		mParentMenu = mainWindow->menuTools;

		// début --- actions ---
		actionNew_Position = new QAction(tr("New_Position"), this);
		if (actionNew_Position)
			connect(actionNew_Position, SIGNAL(triggered()), this, SLOT(New_Position()));

		actionSubdiviser = new QAction(tr("Subdiviser"), this);
		if (actionSubdiviser)
			connect(actionSubdiviser, SIGNAL(triggered()), this, SLOT(Subdiviser()));

		actionUnion = new QAction(tr("Union"), this);
		if (actionUnion)
			connect(actionUnion, SIGNAL(triggered()), this, SLOT(Union()));

		actionInter = new QAction(tr("Intersection"), this);
		if (actionInter)
			connect(actionInter, SIGNAL(triggered()), this, SLOT(Inter()));

		actionMinus = new QAction(tr("Subtraction"), this);
		if (actionMinus)
			connect(actionMinus, SIGNAL(triggered()), this, SLOT(Minus()));
		// fin --- actions ---
	}

	QList<QAction*> actions() const
	{
		return QList<QAction*>()	<< actionNew_Position
									<< NULL
									<< actionSubdiviser
									<< NULL
									<< actionUnion
									<< actionInter
									<< actionMinus;
	}

public slots:
	void New_Position();
	void Subdiviser();
	void Union();
	void Inter();
	void Minus();

private:
	QAction *actionNew_Position;
	QAction *actionSubdiviser;
	QAction *actionUnion;
	QAction *actionInter;
	QAction *actionMinus;
};

#endif // BUILD_component_Boolean_Operations

#endif // HEADER_MEPP_COMPONENT_BOOLEAN_OPERATIONS_PLUGIN_INTERFACE_H
