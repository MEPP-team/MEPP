/*!
 * \file mainwindow.hxx
 * \brief mainwindow file.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2010
 */
#ifndef HEADER_MAINWINDOW
#define HEADER_MAINWINDOW

#ifndef _MSC_VER
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
#include <QtGui>
#ifndef _MSC_VER
#pragma GCC diagnostic warning "-Wuninitialized"
#endif

#include "ui_mainwindow.h"

#include "viewer.hxx"
#include "mepp_action.hxx"

class mepp_component_plugin_interface;

/*! 
 * \class mainwindow
 * \brief mainwindow class.
 */
class mainwindow : public QMainWindow, /*private*/public Ui::mainwindow
{
    Q_OBJECT

    public:
		/*!
		 * \brief Constructor.
		 *
		 * \param parent : parent window.
		 */
        mainwindow(QMainWindow *parent = 0);
		/*!
		 * \brief Destructor.
		 */
		~mainwindow();

		/*!
		 * \fn QAction *get_actionNewEmpty()
		 * \brief Get a QAction for emit a request to create a viewer with an empty mesh.
		 *
		 * \return QAction*.
		 */
		QAction *get_actionNewEmpty() { return actionNewEmpty; }
		/*!
		 * \fn QAction *get_actionAddEmpty()
		 * \brief Get a QAction for emit a request to add an empty mesh in the current active viewer.
		 *
		 * \return QAction*.
		 */
		QAction *get_actionAddEmpty() { return actionAddEmpty; }

		/*!
		 * \fn QAction *get_mainwindowActionOpen()
		 * \brief Get a QAction for emit a request to open mesh file(s) in a new viewer in normal mode.
		 *
		 * \return mepp_action*.
		 */
		mepp_action *get_mainwindowActionOpen() { return mainwindowActionOpen; }

		/*!
		 * \fn QAction *get_mainwindowActionOpen_space()
		 * \brief Get a QAction for emit a request to open mesh file(s) in a new viewer in space mode.
		 *
		 * \return mepp_action*.
		 */
		mepp_action *get_mainwindowActionOpen_space() { return mainwindowActionOpen_space; }
		/*!
		 * \fn QAction *get_mainwindowActionOpen_time()
		 * \brief Get a QAction for emit a request to open mesh file(s) in a new viewer in time mode.
		 *
		 * \return mepp_action*.
		 */
		mepp_action *get_mainwindowActionOpen_time() { return mainwindowActionOpen_time; }

		/*!
		 * \fn QAction *get_mainwindowActionOpen_and_Add_space()
		 * \brief Get a QAction for emit a request to open and add mesh file(s) in the current active viewer in space mode.
		 *
		 * \return mepp_action*.
		 */
		mepp_action *get_mainwindowActionOpen_and_Add_space() { return mainwindowActionOpen_and_Add_space; }
		/*!
		 * \fn QAction *get_mainwindowActionOpen_and_Add_time()
		 * \brief Get a QAction for emit a request to open and add mesh file(s) in the current active viewer in time mode.
		 *
		 * \return mepp_action*.
		 */
		mepp_action *get_mainwindowActionOpen_and_Add_time() { return mainwindowActionOpen_and_Add_time; }

		/*!
		 * \fn QAction *get_mainwindowActionSave_As()
		 * \brief Get a QAction for emit a request to save polyhedron(s) in file(s).
		 *
		 * \return mepp_action*.
		 */
		mepp_action *get_mainwindowActionSave_As() { return mainwindowActionSave_As; }

		/*!
		 * \fn int loadFile(const QString &fileName, int loadType, typeFuncOpenSave f)
		 * \brief Open and load a mesh in a new viewer.
		 *
		 * \param fileName the mesh filename.
		 * \param loadType the type of load (Normal, Space or Time).
		 * \param f a pointer to a typeFuncOpenSave function to load the file (default function if NULL).
		 * \return 0 if loading is ok.
		 */
		int loadFile(const QString &fileName, int loadType, typeFuncOpenSave f);
		/*!
		 * \fn int addFile(Viewer *viewer, const QString &fileName, int loadType, typeFuncOpenSave f)
		 * \brief Open and load a mesh in a viewer.
		 *
		 * \param viewer the viewer.
		 * \param fileName the mesh filename.
		 * \param loadType the type of load (Normal, Space or Time).
		 * \param f a pointer to a typeFuncOpenSave function to load the file (default function if NULL).
		 * \return 0 if loading is ok.
		 */
		int addFile(Viewer *viewer, const QString &fileName, int loadType, typeFuncOpenSave f);

		/*!
		 * \fn QWidget *activeMdiChild()
		 * \brief Return the current active viewer.
		 *
		 * \return the current active viewer or 0 if not.
		 */
		QWidget *activeMdiChild()
		{
			if (QMdiSubWindow *activeSubWindow = mdiArea->activeSubWindow())
				return qobject_cast<QWidget *>(activeSubWindow->widget());

			return 0;
		}

		/*!
		 * \fn void update_mesh_properties(bool update_component=false, bool update_boundary=false)
		 * \brief Update current active mesh properties (component(s), vertice(s), facet(s), edge(s), boundary(ies), genus) in status bar.
		 *
		 * \param update_component update component(s) information (default is none).
		 * \param update_boundary update boundary(ies) information (default is none).
		 */
		void update_mesh_properties(bool update_component=false, bool update_boundary=false);

		/*!
		 * \fn readSettings()
		 * \brief Read Mepp settings: mainwindow pos and size, recent files list, tree location (init).
		 */
		void readSettings();

		/*!
		 * \fn writeSettings()
		 * \brief Write Mepp settings: mainwindow pos and size, recent files list, tree location (exit).
		 */
		void writeSettings();

		// copy/paste viewpoint
		bool viewpointCopied;							//!< true if a camera viewpoint is copied
		qglviewer::Vec copyCameraPosition;				//!< the copied camera position
		qglviewer::Quaternion copyCameraOrientation;	//!< the copied camera orientation

		QString savePNGLocation, saveAVILocation;		//!< default paths to save PNG snapshots or AVI video

		Viewer *lastViewerCreated;						//!< last viewer created pointeur

		QString saveLocation;							//!< default or last save location

	protected:
		/*!
		 * \brief Close all viewers (event).
		 *
		 * \param event the event.
		 */
		void closeEvent(QCloseEvent *event);

		/*!
		 * \fn void clearMenu(QMenu*)
		 * \brief Clear a menu (recursive).
		 *
		 * \param menu a menu.
		 */
		void clearMenu(QMenu*);

		/*!
		 * \fn void loadPlugins()
		 * \brief Load all plugins (components).
		 */
		void loadPlugins();
		/*!
		 * \fn bool initPlugin(QObject*)
		 * \brief Init all plugins (components).
		 *
		 * \param obj a plugin.
		 * \param nu_plugin plugin number.
		 * \return true if ok.
		 */
		bool initPlugin(QObject* obj, int nu_plugin);

    private slots:
		/*!
		 * \fn void updateMenus()
		 * \brief Update all menus and toolbars.
		 */
		void updateMenus();
		/*!
		 * \fn void updateWindowMenu()
		 * \brief Update list of viewers in Window menu.
		 */
		void updateWindowMenu();

		/*!
		 * \fn on_actionNew_triggered
		 * \brief Create a viewer with an intermal mesh: a simple tetrahedron (from menu).
		 */
		void on_actionNew_triggered();

		/*!
		 * \fn on_actionOpen_triggered
		 * \brief Create a viewer and load mesh file(s) (from menu).
		 */
		void on_actionOpen_triggered();
		/*!
		 * \fn on_actionOpen_space_triggered
		 * \brief Create a viewer and load mesh file(s) in space mode (from menu).
		 */
		void on_actionOpen_space_triggered();
		/*!
		 * \fn on_actionOpen_time_triggered
		 * \brief Create a viewer and load mesh file(s) in time mode (from menu).
		 */
		void on_actionOpen_time_triggered();
		/*!
		 * \fn on_actionOpen_and_Add_space_triggered
		 * \brief Load mesh file(s) in space mode in the current active viewer (from menu).
		 */
		void on_actionOpen_and_Add_space_triggered();
		/*!
		 * \fn on_actionOpen_and_Add_time_triggered
		 * \brief Load mesh file(s) in time mode in the current active viewer (from menu).
		 */
		void on_actionOpen_and_Add_time_triggered();

		/*!
		 * \fn on_actionSave_As_triggered
		 * \brief Save polyhedron(s) from the current viewer in file(s) (from menu).
		 */
		void on_actionSave_As_triggered();

		/*!
		 * \fn on_actionClose_triggered
		 * \brief Close the current viewer (from menu).
		 */
		void on_actionClose_triggered();

		/*!
		 * \fn on_actionDelete_triggered()
		 * \brief Delete current mesh - only in Space or Time mode (from menu).
		 */
		void on_actionDelete_triggered();

		/*!
		 * \fn on_actionClone_triggered
		 * \brief Clone the current viewer (from menu).
		 */
		void on_actionClone_triggered();

		void on_actionOpen_texture_triggered();
		void on_actionTexture_settings_triggered();
		void on_actionTexture_to_Vertex_Color_triggered();

		/*!
		 * \brief Close the current viewer (from window menu).
		 */
		void on_actionClose_Window_triggered();
		/*!
		 * \fn on_actionClose_All_triggered
		 * \brief Close all viewers (from window menu).
		 */
		void on_actionClose_All_triggered();

		/*!
		 * \fn on_actionChange_MDI_View_Mode_triggered
		 * \brief Change type of MDI view mode (tabbed/subwindow view).
		 */
		void on_actionChange_MDI_View_Mode_triggered();
		/*!
		 * \fn on_actionChange_Viewer_Mode_Space_Time_triggered
		 * \brief Change type of viewer mode (space/time mode).
		 */
		void on_actionChange_Viewer_Mode_Space_Time_triggered();

		/*!
		 * \fn on_actionExit_triggered
		 * \brief Quit Mepp.
		 */
		void on_actionExit_triggered();

		/*!
		 * \fn actionNewEmpty_slot()
		 * \brief Create a viewer with an empty mesh (from action).
		 */
		void actionNewEmpty_slot();
		/*!
		 * \fn actionAddEmpty_slot()
		 * \brief Add an empty mesh in the current active viewer (from action).
		 */
		void actionAddEmpty_slot();

		/*!
		 * \fn actionOpen_slot(QString, QString, int, typeFuncOpenSave)
		 * \brief Open mesh file(s) in a new viewer in normal mode (from action).
		 *
		 * \param title dialog box title.
		 * \param typeFiles filter for file(s), for example: "OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply);;X3D files (*.x3d);;ALL files (*.*)".
		 * \param loadType the type of load (Normal, Space or Time).
		 * \param f a pointer to a typeFuncOpenSave function to load the file(s) (default function if NULL).
		 */
		void actionOpen_slot(QString, QString, int, typeFuncOpenSave);
		/*!
		 * \fn actionOpen_space_slot(QString, QString, typeFuncOpenSave)
		 * \brief Open mesh file(s) in a new viewer in space mode (from action).
		 *
		 * \param title dialog box title.
		 * \param typeFiles filter for file(s), for example: "OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply);;X3D files (*.x3d);;ALL files (*.*)".
		 * \param f a pointer to a typeFuncOpenSave function to load the file(s) (default function if NULL).
		 */
		void actionOpen_space_slot(QString, QString, typeFuncOpenSave);
		/*!
		 * \fn actionOpen_time_slot(QString, QString, typeFuncOpenSave)
		 * \brief Open mesh file(s) in a new viewer in time mode (from action).
		 *
		 * \param title dialog box title.
		 * \param typeFiles filter for file(s), for example: "OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply);;X3D files (*.x3d);;ALL files (*.*)".
		 * \param f a pointer to a typeFuncOpenSave function to load the file(s) (default function if NULL).
		 */
		void actionOpen_time_slot(QString, QString, typeFuncOpenSave);
		/*!
		 * \brief Open and add mesh file(s) in the current active viewer in space mode (from action).
		 *
		 * \param title dialog box title.
		 * \param typeFiles filter for file(s), for example: "OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply);;X3D files (*.x3d);;ALL files (*.*)".
		 * \param f a pointer to a typeFuncOpenSave function to load the file(s) (default function if NULL).
		 */
		void actionOpen_and_Add_space_slot(QString, QString, typeFuncOpenSave);
		/*!
		 * \fn actionOpen_and_Add_time_slot(QString, QString, typeFuncOpenSave)
		 * \brief Open and add mesh file(s) in the current active viewer in time mode (from action).
		 *
		 * \param title dialog box title.
		 * \param typeFiles filter for file(s), for example: "OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply);;X3D files (*.x3d);;ALL files (*.*)".
		 * \param f a pointer to a typeFuncOpenSave function to load the file(s) (default function if NULL).
		 */
		void actionOpen_and_Add_time_slot(QString, QString, typeFuncOpenSave);

		/*!
		 * \fn actionSave_As_slot(QString, QString, typeFuncOpenSave)
		 * \brief Save polyhedron(s) in file(s) (from action).
		 *
		 * \param title dialog box title.
		 * \param typeFiles filter for file(s), for example: "OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply);;X3D files (*.x3d);;ALL files (*.*)".
		 * \param f a pointer to a typeFuncOpenSave function to save the file(s) (default function if NULL).
		 */
		void actionSave_As_slot(QString, QString, typeFuncOpenSave);

		/*!
		 * \fn on_actionAbout_triggered()
		 * \brief Open about Mepp dialog box (from menu).
		 */
		void on_actionAbout_triggered();

		/*!
		 * \fn on_actionAbout_QGLViewer_triggered()
		 * \brief Open about QGLViewer dialog box (from menu).
		 */
		void on_actionAbout_QGLViewer_triggered();

		/*!
		 * \fn popupAboutBox(QString title, QString html_resource_name)
		 * \brief Popup an about dialog box.
		 *
		 * \param title dialog box title.
		 * \param html_resource_name an html file.
		 */
		void popupAboutBox(QString title, QString html_resource_name);
		/*!
		 * \fn on_actionAbout_CGAL_triggered()
		 * \brief Open about CGAL dialog box (from menu).
		 */
		void on_actionAbout_CGAL_triggered();

		/*!
		 * \fn openRecentFile()
		 * \brief Open a recent file (from menu).
		 */
		void openRecentFile();
		/*!
		 * \fn setActiveSubWindow(QWidget *window)
		 * \brief Set the active viewer.
		 *
		 * \param window a viewer.
		 */
		void setActiveSubWindow(QWidget *window);

		// rendering options
		/*!
		 * \fn on_actionRender_Point_triggered()
		 * \brief Render mesh in point mode (from menu).
		 */
		void on_actionRender_Point_triggered();
		/*!
		 * \fn on_actionRender_Line_triggered()
		 * \brief Render mesh in line mode (from menu).
		 */
		void on_actionRender_Line_triggered();
		/*!
		 * \fn on_actionRender_Fill_triggered()
		 * \brief Render mesh in fill mode (from menu).
		 */
		void on_actionRender_Fill_triggered();

		/*!
		 * \fn on_actionSuperimpose_Edges_triggered()
		 * \brief Superimpose edges of mesh(es) (from menu).
		 */
		void on_actionSuperimpose_Edges_triggered();
		/*!
		 * \fn on_actionSuperimpose_Vertices_triggered()
		 * \brief Superimpose vertices of mesh(es) (from menu).
		 */
		void on_actionSuperimpose_Vertices_triggered();
		/*!
		 * \fn on_actionSuperimpose_Vertices_big_triggered()
		 * \brief Superimpose vertices (bigger mode) of mesh(es) (from menu).
		 */
		void on_actionSuperimpose_Vertices_big_triggered();

		/*!
		 * \fn on_actionVertex_Color_triggered()
		 * \brief View mesh(es) in vertex color mode (from menu).
		 */
		void on_actionVertex_Color_triggered();
		/*!
		 * \fn on_actionFace_Color_triggered()
		 * \brief View mesh(es) in face color mode (from menu).
		 */
		void on_actionFace_Color_triggered();

		void on_actionTexture_Mode_triggered();

		/*!
		 * \fn on_actionLighting_triggered()
		 * \brief Toggle lighting (on/off) (from menu).
		 */
		void on_actionLighting_triggered();
		/*!
		 * \fn on_actionSmooth_Shading_triggered()
		 * \brief Toggle smooth shading (on/off) (from menu).
		 */
		void on_actionSmooth_Shading_triggered();

		/*!
		 * \fn on_actionAntialiasing_triggered()
		 * \brief Toggle antialiasing (on/off) (from menu).
		 */
		void on_actionAntialiasing_triggered();
		/*!
		 * \fn on_actionCulling_triggered()
		 * \brief Toggle culling (on/off) (from menu).
		 */
		void on_actionCulling_triggered();
		
		// color options
		/*!
		 * \fn on_actionBackground_color_triggered()
		 * \brief Choose and set background color (from menu).
		 */
		void on_actionBackground_color_triggered();
		/*!
		 * \fn on_actionVertex_color_triggered()
		 * \brief Choose and set vertex color (from menu).
		 */
		void on_actionVertex_color_triggered();
		/*!
		 * \fn on_actionEdge_color_triggered()
		 * \brief Choose and set edge color (from menu).
		 */
		void on_actionEdge_color_triggered();
		/*!
		 * \fn on_actionFace_color_triggered()
		 * \brief Choose and set face color (from menu).
		 */
		void on_actionFace_color_triggered();
		/*!
		 * \fn on_actionMaterial_triggered()
		 * \brief Choose and set material (from menu).
		 */
		void on_actionMaterial_triggered();
		// color options

		// show options
		/*!
		 * \fn on_actionShow_FPS_triggered()
		 * \brief Toggle fps display information (on/off) (from menu).
		 */
		void on_actionShow_FPS_triggered();
		/*!
		 * \fn on_actionShow_axis_triggered()
		 * \brief Toggle axis display (on/off) (from menu).
		 */
		void on_actionShow_axis_triggered();
		/*!
		 * \fn on_actionShow_grid_triggered()
		 * \brief Toggle grid display (on/off) (from menu).
		 */
		void on_actionShow_grid_triggered();

		/*!
		 * \fn on_actionShow_normals_triggered()
		 * \brief Toggle normals display (on/off) (from menu).
		 */
		void on_actionShow_normals_triggered();

		/*!
		 * \fn on_actionBounding_box_triggered()
		 * \brief Toggle bounding box display (on/off) (from menu).
		 */
		void on_actionBounding_box_triggered();
		/*!
		 * \fn on_actionBounding_box_when_moving_triggered()
		 * \brief Toggle bounding box display when moving (on/off) (from menu).
		 */
		void on_actionBounding_box_when_moving_triggered();
		// show options

		// view options
		/*!
		 * \fn on_actionReset_viewpoint_triggered()
		 * \brief Reset view point (from menu).
		 */
		void on_actionReset_viewpoint_triggered();
		/*!
		 * \fn on_actionCopy_viewpoint_triggered()
		 * \brief Copy view point (from menu).
		 */
		void on_actionCopy_viewpoint_triggered();
		/*!
		 * \fn on_actionPaste_viewpoint_triggered()
		 * \brief Past view point (from menu).
		 */
		void on_actionPaste_viewpoint_triggered();

		/*!
		 * \fn on_actionCenter_all_objects_triggered()
		 * \brief Center all meshes in space mode (from menu).
		 */
		void on_actionCenter_all_objects_triggered();

		/*!
		 * \fn on_actionCouplingRotations_triggered()
		 * \brief Coupling rotation of meshes (from menu).
		 */
		void on_actionCouplingRotations_triggered();

		/*!
		 * \fn on_actionVBO_triggered()
		 * \brief Toggle display lists mode (on/off) (from menu).
		 */
		void on_actionVBO_triggered();
		// view options

		// capture options
		/*!
		 * \fn on_actionScreenshot_triggered()
		 * \brief Take a bitmap screenshot (from menu).
		 */
		void on_actionScreenshot_triggered();
		/*!
		 * \fn on_actionScreenshot_sequence_triggered()
		 * \brief Take a bitmaps sequence screenshot (from menu).
		 */
		void on_actionScreenshot_sequence_triggered();

		/*!
		 * \fn on_actionClipboard_screenshot_triggered()
		 * \brief Copy a bitmap screenshot in clipboard (from menu).
		 */
		void on_actionClipboard_screenshot_triggered();
		// capture options

		// dynamic options
		/*!
		 * \fn on_actionParams_triggered()
		 * \brief Set params for dynamic time sequence (from menu).
		 */
		void on_actionParams_triggered();

		/*!
		 * \fn on_actionReverse_start_triggered()
		 * \brief Start dynamic time sequence in reverse order (from menu).
		 */
		void on_actionReverse_start_triggered();
		/*!
		 * \fn on_actionReverse_start_loop_triggered()
		 * \brief Start dynamic time sequence in reverse order with loop (from menu).
		 */
		void on_actionReverse_start_loop_triggered();
		/*!
		 * \fn on_actionStart_triggered()
		 * \brief Start dynamic time sequence (from menu).
		 */
		void on_actionStart_triggered();
		/*!
		 * \fn on_actionStart_loop_triggered()
		 * \brief Start dynamic time sequence with loop (from menu).
		 */
		void on_actionStart_loop_triggered();
		/*!
		 * \fn on_actionStop_triggered()
		 * \brief Stop dynamic time sequence (from menu).
		 */
		void on_actionStop_triggered();

		/*!
		 * \fn on_actionDynFirst_triggered()
		 * \brief Go to first position of dynamic time sequence (from menu).
		 */
		void on_actionDynFirst_triggered();
		/*!
		 * \fn on_actionDynPrevious_triggered()
		 * \brief Go to previous position of dynamic time sequence (from menu).
		 */
		void on_actionDynPrevious_triggered();
		/*!
		 * \fn on_actionDynNext_triggered()
		 * \brief Go to next position of dynamic time sequence (from menu).
		 */
		void on_actionDynNext_triggered();
		/*!
		 * \fn on_actionDynLast_triggered()
		 * \brief Go to last position of dynamic time sequence (from menu).
		 */
		void on_actionDynLast_triggered();
		// dynamic options

	private:
		Viewer* the_viewer;									//!< last new/load viewer
		QGLViewer* aboutQGLViewer;							//!< QGLViewer for about dialog box
		QList<mepp_component_plugin_interface *> lplugin;	//!< list of plugins/components

		QSignalMapper *windowMapper;

		/*!
		 * \fn createActions()
		 * \brief Create all Mepp actions (init).
		 */
		void createActions();
		/*!
		 * \fn createMenus()
		 * \brief Create all Mepp menus (init).
		 */
		void createMenus();
		/*!
		 * \fn createStatusBar()
		 * \brief Create Mepp status bar (init).
		 */
		void createStatusBar();

		QLabel *STATUSBAR_COMPONENTS;	//!< label for status bar
		QLabel *STATUSBAR_VERTICES;		//!< label for status bar
		QLabel *STATUSBAR_FACETS;		//!< label for status bar
		QLabel *STATUSBAR_EDGES;		//!< label for status bar
		QLabel *STATUSBAR_BOUNDARIES;	//!< label for status bar
		QLabel *STATUSBAR_GENUS;		//!< label for status bar

		/*!
		 * \fn setCurrentFile()
		 * \brief Set a mesh filename in recent files list.
		 *
		 * \param fileName filename of the mesh.
		 */
		void setCurrentFile(const QString &fileName);
		/*!
		 * \fn strippedName(const QString &fullFileName)
		 * \brief Return the filename of a mesh without path.
		 *
		 * \param fullFileName filename of the mesh.
		 * \return QString.
		 */
		QString strippedName(const QString &fullFileName);

		QAction *separatorAct_menuFile;				//!< separator in file menu
		QAction *separatorAct_menuWindow;			//!< separator in window menu

		/*!
		 * \fn updateRecentFileActions()
		 * \brief Update Mepp recent files list.
		 */
		void updateRecentFileActions();
		enum { MaxRecentFiles = 5 };				//!< only 5 mesh files in recent list
		QAction *recentFileActs[MaxRecentFiles];	//!< tab of mesh files recent list

		QAction *actionNewEmpty;
		QAction *actionAddEmpty;

		mepp_action *mainwindowActionOpen, *mainwindowActionOpen_space, *mainwindowActionOpen_time;
		mepp_action *mainwindowActionOpen_and_Add_space, *mainwindowActionOpen_and_Add_time;
		mepp_action *mainwindowActionSave_As;

		QDockWidget *dockComponents, *dockDirView;
		QMainWindow *inner;

		QFileSystemModel *model;
		QSortFilterProxyModel *proxyModel;
		QTreeView *tree;

		QString treeLocation;
		QString openLocation;

		int m_dockComponents_MinimumWidth, m_dockDirView_MinimumWidth;
};

#endif