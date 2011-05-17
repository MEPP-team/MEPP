#ifndef HEADER_MEPP_COMPONENT_COMPRESSION_VALENCE_PLUGIN_INTERFACE_H
#define HEADER_MEPP_COMPONENT_COMPRESSION_VALENCE_PLUGIN_INTERFACE_H

#include <QObject>

#include <mepp_config.h>
#ifdef BUILD_component_Compression_Valence

#include "../../../../mepp/mepp_component_plugin_interface.h"

#include <QAction>
#include <QtPlugin>

/**
 \class	mepp_component_Compression_Valence_plugin

 \brief	Mepp component compression valence plugin. 

 */

class mepp_component_Compression_Valence_plugin : 
  public QObject,
  public mepp_component_plugin_interface
{

	Q_OBJECT
	Q_INTERFACES(mepp_component_plugin_interface);

	public:

		/**
		 \fn	mepp_component_Compression_Valence_plugin()
		
		 \brief	Default constructor.

		 */

		mepp_component_Compression_Valence_plugin() : mepp_component_plugin_interface() {}

		/**
		 \fn	~mepp_component_Compression_Valence_plugin()
		
		 \brief	Finaliser.
		

		 */

		~mepp_component_Compression_Valence_plugin()
		{
			delete actionCompress; delete actionOpen_P3D_file; delete actionDecompress_all;
			delete actionDecompress_one_level; delete actionDecompress_precedent_level; delete actionDecompress_go_to_specific_level;
			delete actionDecompress_mesh_sequence_on_off;
			delete actionJCW;
			delete actionJCWdecompress;
			delete actionJCWdecompress_without_extraction;
		}

		/**
		 \fn	void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		
		 \brief	Initialises this object.

		
		 \param [in,out]	mainWindow	If non-null, the main window.
		 \param [in,out]	lw		  	If non-null, the lw.
		 */

		void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
		{
			this->mw = mainWindow;
			this->lwindow = lw;
			this->mPluginName = this->metaObject()->className();

			// choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
			// menuCompression, menuWatermaking, menuExamples
			mParentMenu = mainWindow->menuCompression;

			// debut --- actions ---
			actionCompress = new QAction(tr("Compression"), this);
			if (actionCompress)
				connect(actionCompress, SIGNAL(triggered()), this, SLOT(OnCompress()));

			actionOpen_P3D_file = new QAction(tr("Open P3D file"), this);
			if (actionOpen_P3D_file)
				connect(actionOpen_P3D_file, SIGNAL(triggered()), this, SLOT(load_P3D_file()));

			actionDecompress_all = new QAction(tr("Decompression : all LoDs"), this);
			if (actionDecompress_all)
				connect(actionDecompress_all, SIGNAL(triggered()), this, SLOT(OnDecompress_all()));

			actionDecompress_one_level = new QAction(tr("Decompression : next LoD"), this);
			if (actionDecompress_one_level)
				connect(actionDecompress_one_level, SIGNAL(triggered()), this, SLOT(OnDecompress_one_level()));

			actionDecompress_precedent_level = new QAction(tr("Decompression : previous LoD"), this);
			if (actionDecompress_precedent_level)
				connect(actionDecompress_precedent_level, SIGNAL(triggered()), this, SLOT(OnDecompress_precedent_level()));

			actionDecompress_go_to_specific_level = new QAction(tr("Decompression : go to specific LoD"), this);
			if (actionDecompress_go_to_specific_level)
				connect(actionDecompress_go_to_specific_level, SIGNAL(triggered()), this, SLOT(OnDecompress_go_to_specific_level()));	
			
			actionJCW = new QAction(tr("JCW - Compression and Watermark embedding"), this);
			if (actionJCW)
				connect(actionJCW, SIGNAL(triggered()), this, SLOT(OnJCW()));
			
			actionJCWdecompress = new QAction(tr("JCW - Decompression and Watermark extraction : next LoD"), this);
			if (actionJCWdecompress)
				connect(actionJCWdecompress, SIGNAL(triggered()), this, SLOT(OnJCWdecompress()));
			
			actionJCWdecompress_without_extraction = new QAction(tr("JCW - Decompression without extraction : next LoD"), this);
			if (actionJCWdecompress_without_extraction)
				connect(actionJCWdecompress_without_extraction, SIGNAL(triggered()), this, SLOT(OnJCWdecompress_without_extraction()));

			actionDecompress_mesh_sequence_on_off = new QAction(tr("Activate/Desactivate mesh sequence generation"), this);
			if (actionDecompress_mesh_sequence_on_off)
				connect(actionDecompress_mesh_sequence_on_off, SIGNAL(triggered()), this, SLOT(OnDecompress_mesh_sequence_on_off()));

			// fin --- actions ---
		}

		/**
		 \fn	QList<QAction*> actions() const
		
		 \brief	Gets the actions.
		
		
		 \return	null if it fails, else a list of.
		 */

		QList<QAction*> actions() const
		{
			return QList<QAction*>()	<< actionCompress
										<< actionJCW
										<< NULL
										<< actionOpen_P3D_file
										<< NULL
										<< actionDecompress_all
										<< actionDecompress_one_level
										<< actionDecompress_precedent_level
										<< actionDecompress_go_to_specific_level
										<< NULL																														
										<< actionJCWdecompress
										<< actionJCWdecompress_without_extraction
										<< NULL
										<< actionDecompress_mesh_sequence_on_off;
		}

		/**
		 \fn	virtual void pre_draw()
		
		 \brief	Pre draw.

		 */

		virtual void pre_draw() {}

		/**
		 \fn	virtual void post_draw()
		
		 \brief	Posts the draw.

		 */

		virtual void post_draw() {}

		/**
		 \fn	virtual void pre_draw_all_scene()
		
		 \brief	Pre draw all scene.

		 */

		virtual void pre_draw_all_scene() {}

		/**
		 \fn	virtual void post_draw_all_scene()
		
		 \brief	Posts the draw all scene.
		
		 */

		virtual void post_draw_all_scene() {}

		/**
		 \fn	virtual void OnMouseLeftDown(QMouseEvent *event)
		
		 \brief	Executes the mouse left down action.

		
		 \param [in,out]	event	If non-null, the event.
		 */

		virtual void OnMouseLeftDown(QMouseEvent *event) {}
		virtual void OnMouseLeftUp(QMouseEvent *event);

		/**
		 \fn	virtual void OnMouseRightDown(QMouseEvent *event)
		
		 \brief	Executes the mouse right down action.

		
		 \param [in,out]	event	If non-null, the event.
		 */

		virtual void OnMouseRightDown(QMouseEvent *event) {}
		virtual void OnMouseRightUp(QMouseEvent *event);

		/**
		 \fn	virtual void OnMouseMotion(QMouseEvent *event)
		
		 \brief	Executes the mouse motion action.
		
		
		 \param [in,out]	event	If non-null, the event.
		 */

		virtual void OnMouseMotion(QMouseEvent *event) {}
		virtual void OnMouseWheel(QWheelEvent *event);

		/**
		 \fn	virtual void OnKeyPress(QKeyEvent *event)
		
		 \brief	Executes the key press action.
		
		
		 \param [in,out]	event	If non-null, the event.
		 */

		virtual void OnKeyPress(QKeyEvent *event) {}

		/**
		 \fn	virtual void OnKeyRelease(QKeyEvent *event)
		
		 \brief	Executes the key release action.
		
		
		 \param [in,out]	event	If non-null, the event.
		 */

		virtual void OnKeyRelease(QKeyEvent *event) {}

		public slots:

		/**
		 \fn	void mepp_component_Compression_Valence_plugin::OnCompress();
		
		 \brief	Executes the compress action.
		
		 */

		void OnCompress();

		/**
		 \fn	void load_P3D_file();
		
		 \brief	Loads the p3d file.
		

		 */

		void load_P3D_file();

		/**
		 \fn	void OnDecompress_all();
		
		 \brief	Executes the decompress all action.
		
		 */

		void OnDecompress_all();

		/**
		 \fn	void OnDecompress_one_level();
		
		 \brief	Executes the decompress one level action.

		 */

		void OnDecompress_one_level();

		/**
		 \fn	void OnDecompress_precedent_level();
		
		 \brief	Executes the decompress precedent level action.
		

		 */

		void OnDecompress_precedent_level();

		/**
		 \fn	void OnDecompress_go_to_specific_level();
		
		 \brief	Executes the decompress go to specific level action.

		 */

		void OnDecompress_go_to_specific_level();

		/**
		 \fn	void OnDecompress_mesh_sequence_on_off();
		
		 \brief	Executes the decompress mesh sequence on off action.

		 */

		void OnDecompress_mesh_sequence_on_off();

		/**
		 \fn	void OnJCW();
		
		 \brief	Executes the Joint Compression and Watermark (JCW) action.

		 */

		void OnJCW();

		/**
		 \fn	void OnJCWdecompres();
		
		 \brief	Executes the JCW decompres action.

		 */

		void OnJCWdecompress();

		/**
		 \fn	void OnJCWdecompress_without_extraction();
		
		 \brief	Executes the JCW decompress without extraction action.
		
		 */

		void OnJCWdecompress_without_extraction();

	protected:
		void ShowText(void);
		//void WriteInfo(void);

	private:
		QAction *actionCompress;
		QAction *actionOpen_P3D_file;

		QAction *actionDecompress_all, *actionDecompress_one_level, *actionDecompress_precedent_level, *actionDecompress_go_to_specific_level;
		QAction *actionDecompress_mesh_sequence_on_off;
		QAction *actionJCW;
		QAction *actionJCWdecompress, * actionJCWdecompress_without_extraction;
};

#endif

#endif
