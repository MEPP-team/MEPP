#ifndef HEADER_MEPP_COMPONENT_COMPRESSION_VALENCE_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_COMPRESSION_VALENCE_PLUGIN_INTERFACE_H

#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

class mepp_component_Compression_Valence_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{
	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);

	public:
		mepp_component_Compression_Valence_plugin() : mepp_component_plugin_interface() {}
		~mepp_component_Compression_Valence_plugin()
		{
			delete actionCompress; delete actionOpen_P3D_file; delete actionDecompress_all;
			delete actionDecompress_one_level; delete actionDecompress_precedent_level; delete actionDecompress_go_to_specific_level;
			delete actionDecompress_mesh_sequence_on_off;
			delete actionJCW;
			delete actionJCWdecompres;
			delete actionJCWdecompress_without_extraction;
		}

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();

			// choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
			// menuCompression, menuWatermaking, menuExamples
			mParentMenu = mainWindow->menuCompression;

			// debut --- actions ---
			actionCompress = new QAction(tr("Compress"), this);
			if (actionCompress)
				connect(actionCompress, SIGNAL(triggered()), this, SLOT(OnCompress()));

			actionOpen_P3D_file = new QAction(tr("Open P3D file"), this);
			if (actionOpen_P3D_file)
				connect(actionOpen_P3D_file, SIGNAL(triggered()), this, SLOT(load_P3D_file()));

			actionDecompress_all = new QAction(tr("Decompress: all"), this);
			if (actionDecompress_all)
				connect(actionDecompress_all, SIGNAL(triggered()), this, SLOT(OnDecompress_all()));

			actionDecompress_one_level = new QAction(tr("Decompress: one level"), this);
			if (actionDecompress_one_level)
				connect(actionDecompress_one_level, SIGNAL(triggered()), this, SLOT(OnDecompress_one_level()));

			actionDecompress_precedent_level = new QAction(tr("Decompress: precedent level"), this);
			if (actionDecompress_precedent_level)
				connect(actionDecompress_precedent_level, SIGNAL(triggered()), this, SLOT(OnDecompress_precedent_level()));

			actionDecompress_go_to_specific_level = new QAction(tr("Decompress: go to specific level"), this);
			if (actionDecompress_go_to_specific_level)
				connect(actionDecompress_go_to_specific_level, SIGNAL(triggered()), this, SLOT(OnDecompress_go_to_specific_level()));

			actionDecompress_mesh_sequence_on_off = new QAction(tr("Decompress: mesh sequence on/off"), this);
			if (actionDecompress_mesh_sequence_on_off)
				connect(actionDecompress_mesh_sequence_on_off, SIGNAL(triggered()), this, SLOT(OnDecompress_mesh_sequence_on_off()));
			
			actionJCW = new QAction(tr("JCW: Insertion"), this);
			if (actionJCW)
				connect(actionJCW, SIGNAL(triggered()), this, SLOT(OnJCW()));
			
			actionJCWdecompres = new QAction(tr("JCW: Decompress one level"), this);
			if (actionJCWdecompres)
				connect(actionJCWdecompres, SIGNAL(triggered()), this, SLOT(OnJCWdecompres()));
			
			actionJCWdecompress_without_extraction = new QAction(tr("JCW: Decompress one level without extraction"), this);
			if (actionJCWdecompress_without_extraction)
				connect(actionJCWdecompress_without_extraction, SIGNAL(triggered()), this, SLOT(OnJCWdecompress_without_extraction()));

			// fin --- actions ---
		}

		QList<QAction*> actions() const
		{
			return QList<QAction*>()	<< actionCompress
										<< NULL
										<< actionOpen_P3D_file
										<< NULL
										<< actionDecompress_all
										<< actionDecompress_one_level
										<< actionDecompress_precedent_level
										<< actionDecompress_go_to_specific_level
										<< NULL
										<< actionDecompress_mesh_sequence_on_off
										<< NULL
										<< actionJCW
										<< actionJCWdecompres
										<< actionJCWdecompress_without_extraction;
		}
		
		virtual void pre_draw() {}
		virtual void post_draw() {}
		virtual void pre_draw_all_scene() {}
		virtual void post_draw_all_scene() {}

		virtual void OnMouseLeftDown(QMouseEvent *event) {}
		virtual void OnMouseLeftUp(QMouseEvent *event);
		virtual void OnMouseRightDown(QMouseEvent *event) {}
		virtual void OnMouseRightUp(QMouseEvent *event);
		virtual void OnMouseMotion(QMouseEvent *event) {}
		virtual void OnMouseWheel(QWheelEvent *event);
		virtual void OnKeyPress(QKeyEvent *event) {}
		virtual void OnKeyRelease(QKeyEvent *event) {}

	public slots:
		void OnCompress();
		void load_P3D_file();

		void OnDecompress_all();
		void OnDecompress_one_level();
		void OnDecompress_precedent_level();
		void OnDecompress_go_to_specific_level();

		void OnDecompress_mesh_sequence_on_off();

		void OnJCW();

		void OnJCWdecompres();
		void OnJCWdecompress_without_extraction();

	protected:
		void ShowText(void);
		void WriteInfo(void);

	private:
		QAction *actionCompress;
		QAction *actionOpen_P3D_file;

		QAction *actionDecompress_all, *actionDecompress_one_level, *actionDecompress_precedent_level, *actionDecompress_go_to_specific_level;
		QAction *actionDecompress_mesh_sequence_on_off;
		QAction *actionJCW;
		QAction *actionJCWdecompres, * actionJCWdecompress_without_extraction;
};

#endif

#endif
