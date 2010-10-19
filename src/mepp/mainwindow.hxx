///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef HEADER_MAINWINDOW
#define HEADER_MAINWINDOW

#include <QtGui>

#include "ui_mainwindow.h"

#include "viewer.hxx"
#include "mepp_action.hxx"

class mepp_component_plugin_interface;

class mainwindow : public QMainWindow, /*private*/public Ui::mainwindow
{
    Q_OBJECT

    public:
        mainwindow(QMainWindow *parent = 0);
		~mainwindow();

		QAction *get_actionNewEmpty() { return actionNewEmpty; }
		QAction *get_actionAddEmpty() { return actionAddEmpty; }

		mepp_action *get_mainwindowActionOpen() { return mainwindowActionOpen; }

		mepp_action *get_mainwindowActionOpen_space() { return mainwindowActionOpen_space; }
		mepp_action *get_mainwindowActionOpen_time() { return mainwindowActionOpen_time; }

		mepp_action *get_mainwindowActionOpen_and_Add_space() { return mainwindowActionOpen_and_Add_space; }
		mepp_action *get_mainwindowActionOpen_and_Add_time() { return mainwindowActionOpen_and_Add_time; }

		mepp_action *get_mainwindowActionSave_As() { return mainwindowActionSave_As; }

		int loadFile(const QString &fileName, int loadType, typeFuncOpenSave f);
		int addFile(Viewer *viewer, const QString &fileName, int loadType, typeFuncOpenSave f);

		QWidget *activeMdiChild()
		{
			if (QMdiSubWindow *activeSubWindow = mdiArea->activeSubWindow())
				return qobject_cast<QWidget *>(activeSubWindow->widget());

			return 0;
		}

		void update_mesh_properties(bool update_component=false, bool update_boundary=false);

		// copy/paste viewpoint
		bool viewpointCopied;
		qglviewer::Vec copyCameraPosition;
		qglviewer::Quaternion copyCameraOrientation;

	protected:
		void closeEvent(QCloseEvent *event);

		void clearMenu(QMenu*);

		void loadPlugins();
		bool initPlugin(QObject*);

    private slots:
		void updateMenus();
		void updateWindowMenu();

		void on_actionNew_triggered();

		void on_actionOpen_triggered();
		void on_actionOpen_space_triggered();
		void on_actionOpen_time_triggered();
		void on_actionOpen_and_Add_space_triggered();
		void on_actionOpen_and_Add_time_triggered();

		void on_actionSave_As_triggered();

		void on_actionClose_triggered();

		void on_actionClone_triggered();

		void on_actionClose_Window_triggered();
		void on_actionClose_All_triggered();

		void on_actionChange_MDI_View_Mode_triggered();
		void on_actionChange_Viewer_Mode_Space_Time_triggered();

		void on_actionExit_triggered();

		void actionNewEmpty_slot();
		void actionAddEmpty_slot();

		void actionOpen_slot(QString, QString, int, typeFuncOpenSave);
		void actionOpen_space_slot(QString, QString, typeFuncOpenSave);
		void actionOpen_time_slot(QString, QString, typeFuncOpenSave);
		void actionOpen_and_Add_space_slot(QString, QString, typeFuncOpenSave);
		void actionOpen_and_Add_time_slot(QString, QString, typeFuncOpenSave);

		void actionSave_As_slot(QString, QString, typeFuncOpenSave);

		void on_actionAbout_triggered();

		void on_actionAbout_QGLViewer_triggered();

		void popupAboutBox(QString title, QString html_resource_name);
		void on_actionAbout_CGAL_triggered();

		void openRecentFile();
		void setActiveSubWindow(QWidget *window);

		// rendering options
		void on_actionRender_Point_triggered();
		void on_actionRender_Line_triggered();
		void on_actionRender_Fill_triggered();

		void on_actionSuperimpose_Edges_triggered();
		void on_actionSuperimpose_Vertices_triggered();
		void on_actionSuperimpose_Vertices_big_triggered();

		void on_actionVertex_Color_triggered();
		void on_actionFace_Color_triggered();

		void on_actionLighting_triggered();
		void on_actionSmooth_Shading_triggered();

		void on_actionAntialiasing_triggered();
		void on_actionCulling_triggered();
		
		// color options
		void on_actionBackground_color_triggered();
		void on_actionVertex_color_triggered();
		void on_actionEdge_color_triggered();
		void on_actionFace_color_triggered();
		void on_actionMaterial_triggered();
		// color options

		// show options
		void on_actionShow_FPS_triggered();
		void on_actionShow_axis_triggered();
		void on_actionShow_grid_triggered();

		void on_actionShow_normals_triggered();

		void on_actionBounding_box_triggered();
		void on_actionBounding_box_when_moving_triggered();
		// show options

		// view options
		void on_actionReset_viewpoint_triggered();
		void on_actionCopy_viewpoint_triggered();
		void on_actionPaste_viewpoint_triggered();

		void on_actionCenter_all_objects_triggered();

		void on_actionCouplingRotations_triggered();

		void on_actionVBO_triggered();
		// view options

		// capture options
		void on_actionScreenshot_triggered();
		void on_actionScreenshot_sequence_triggered();

		void on_actionClipboard_screenshot_triggered();
		// capture options

		// dynamic options
		void on_actionParams_triggered();

		void on_actionReverse_start_triggered();
		void on_actionReverse_start_loop_triggered();
		void on_actionStart_triggered();
		void on_actionStart_loop_triggered();
		void on_actionStop_triggered();

		void on_actionDynFirst_triggered();
		void on_actionDynPrevious_triggered();
		void on_actionDynNext_triggered();
		void on_actionDynLast_triggered();

		void on_actionDynDelete_triggered();
		// dynamic options

	private:
		Viewer* the_viewer;
		QGLViewer* aboutQGLViewer;
		QList<mepp_component_plugin_interface *> lplugin;

		QSignalMapper *windowMapper;

		void createActions();
		void createMenus();
		void createStatusBar();

		QLabel *STATUSBAR_COMPONENTS, *STATUSBAR_VERTICES, *STATUSBAR_FACETS, *STATUSBAR_EDGES, *STATUSBAR_BOUNDARIES, *STATUSBAR_GENUS;

		void readSettings();
		void writeSettings();

		void setCurrentFile(const QString &fileName);
		QString strippedName(const QString &fullFileName);

		QAction *separatorAct_menuFile;
		QAction *separatorAct_menuWindow;

		void updateRecentFileActions();
		enum { MaxRecentFiles = 5 };
		QAction *recentFileActs[MaxRecentFiles];

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
};

#endif