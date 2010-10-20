///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#include "mainwindow.hxx"

#define MEPP_VERSION "v0.36d - 20/10/2010"

#include "mepp_component_plugin_interface.h"

typedef boost::shared_ptr<QAction> QActionPtr;

mainwindow::mainwindow(QMainWindow *parent) : QMainWindow(parent)
{
	setupUi(this);

	connect(actionTile, SIGNAL(triggered()), mdiArea, SLOT(tileSubWindows()));
	connect(actionCascade, SIGNAL(triggered()), mdiArea, SLOT(cascadeSubWindows()));

	connect(actionNext, SIGNAL(triggered()), mdiArea, SLOT(activateNextSubWindow()));
	connect(actionPrevious, SIGNAL(triggered()), mdiArea, SLOT(activatePreviousSubWindow()));

	connect(mdiArea, SIGNAL(subWindowActivated(QMdiSubWindow *)), this, SLOT(updateMenus()));

	connect(actionAbout_Qt, SIGNAL(triggered()), qApp, SLOT(aboutQt()));

	windowMapper = new QSignalMapper(this);	// to see
    connect(windowMapper, SIGNAL(mapped(QWidget *)), this, SLOT(setActiveSubWindow(QWidget *)));

	createActions();
	createMenus();
	createStatusBar();
	updateMenus();

	clearMenu(menuComponents);

	readSettings();

	// dock

	// Components
	dockComponents = new QDockWidget(tr(" Components (right-click)"), this);
	dockComponents->setMinimumWidth(220);
	inner = new QMainWindow(dockComponents);
	inner->setWindowFlags(Qt::Widget); // <---------
	dockComponents->setWidget(inner);
	this->addDockWidget(Qt::LeftDockWidgetArea, dockComponents);

	// ---

	// DirView
	dockDirView = new QDockWidget(tr(" Directory view"), this);
	dockDirView->setMinimumWidth(260);
	this->addDockWidget(Qt::RightDockWidgetArea, dockDirView);

		model = new QFileSystemModel;
		QStringList filters;
		filters << "*.off" << "*.obj" << "*.smf" << "*.ply";  // extensions
		model->setNameFilters(filters);
		model->setNameFilterDisables(false);

		proxyModel = new QSortFilterProxyModel;
		proxyModel->setSourceModel(model);

		tree = new QTreeView();
		tree->setModel(proxyModel);	//model
		tree->setSortingEnabled(true);

		tree->setColumnHidden(2, true);
		tree->setColumnWidth(0, 260);
		tree->sortByColumn(0, Qt::AscendingOrder);
		
		//QString location("C:\\_prj_\\MEPP2\\SVN\\data\\");
		QModelIndex index = model->setRootPath(treeLocation);
		QModelIndex proxyIndex = proxyModel->mapFromSource(index);

		tree->scrollTo(proxyIndex);
		tree->setExpanded(proxyIndex, true);
		tree->setCurrentIndex(proxyIndex);

		tree->setSelectionMode(QAbstractItemView::ExtendedSelection);	//QAbstractItemView::MultiSelection
		tree->setDragEnabled(true);

    dockDirView->setWidget(tree);
	// dock
	
	loadPlugins();

	this->setWindowTitle(tr("%1 - %2").arg(MAINWINDOW_TITLE).arg(MEPP_VERSION));
	aboutQGLViewer = new QGLViewer(); // for aboutQGLViewer

	viewpointCopied = false;
	qglviewer::Vec copyCameraPosition(0.f, 0.f, 0.f);
	qglviewer::Quaternion copyCameraOrientation(0, 0, 0, 1); // identity Quaternion

	//setUnifiedTitleAndToolBarOnMac(true);

	setlocale(LC_ALL, "C");

	mdiArea->setMainWindow(this);
	mdiArea->setAcceptDrops(true);
}

mainwindow::~mainwindow()
{
	for (int i=0; i<lplugin.size(); i++)
	{
		delete lplugin[i]->mMenu;
		delete lplugin[i]->mToolBar;

		delete lplugin[i];
	}

	// dock

	// Components
	delete inner;
	delete dockComponents;

	// DirView
	delete tree;
	delete proxyModel;
	delete model;
	delete dockDirView;
	// dock

	delete aboutQGLViewer; // for aboutQGLViewer

	delete windowMapper;

	delete actionNewEmpty;
	delete actionAddEmpty;

	delete mainwindowActionOpen;
	delete mainwindowActionOpen_space;
	delete mainwindowActionOpen_time;
	delete mainwindowActionOpen_and_Add_space;
	delete mainwindowActionOpen_and_Add_time;

	delete mainwindowActionSave_As;

	for (int i = 0; i < MaxRecentFiles; ++i)
	{
        delete recentFileActs[i];
    }

	delete STATUSBAR_COMPONENTS;
	delete STATUSBAR_VERTICES;
	delete STATUSBAR_FACETS;
	delete STATUSBAR_EDGES;
	delete STATUSBAR_BOUNDARIES;
	delete STATUSBAR_GENUS;
}

void mainwindow::clearMenu(QMenu* menu)
{
	Q_FOREACH (QAction* action, menu->actions())
	{
		QMenu* menu = action->menu();
		if (menu)
			clearMenu(menu);
	}
	menu->menuAction()->setEnabled(false);
}

void mainwindow::loadPlugins()
{
  QDir pluginsDir(qApp->applicationDirPath());
  Q_FOREACH (QString fileName, pluginsDir.entryList(QDir::Files))
  {
    if (fileName.contains("mepp_") && QLibrary::isLibrary(fileName))
	{
      QPluginLoader loader;
      loader.setFileName(pluginsDir.absoluteFilePath(fileName));
      QObject *obj = loader.instance();
      if (obj)
	  {
        if (initPlugin(obj))
		{
			qDebug("### Loading \"%s\"...", fileName.toUtf8().data());
			menuComponents->menuAction()->setEnabled(true);
		}
		else
		{
			qDebug("### Duplicate plugin \"%s\". Plugin not loaded.", fileName.toUtf8().data());
			delete obj;
		}
      }
      else
	  {
        qDebug("Error loading \"%s\": %s", qPrintable(fileName), qPrintable(loader.errorString()));

		// Editeur de liens, Système : mettre "Console (/SUBSYSTEM:CONSOLE)" ou "Windows (/SUBSYSTEM:WINDOWS)"
		QMessageBox::critical(this,
                            tr("Cannot load plugin"),
                            tr("Error loading \"%1\": %2")
                            .arg(fileName)
							.arg(loader.errorString())); // pour plugin en release, mettre QT_NO_DEBUG dans C/C++, Préprocesseur
      }
    }
  }
}

bool mainwindow::initPlugin(QObject* obj)
{
	mepp_component_plugin_interface* plugin = qobject_cast<mepp_component_plugin_interface*>(obj);

	// call plugin's init() method
	plugin->init(this, mdiArea->subWindowList());

	for (int i=0; i<lplugin.size(); i++)
		if (plugin->getPluginName() == lplugin[i]->getPluginName())
			return false;

	// menu
	if (plugin->mParentMenu)
	{
		plugin->mMenu = new QMenu(plugin->getPluginName(), plugin->mParentMenu);
		plugin->mParentMenu->menuAction()->setEnabled(true);
		plugin->mParentMenu->addMenu(plugin->mMenu);
	}
	// menu

	// dock
	plugin->mToolBar = new QToolBar(inner);
	plugin->mToolBar->setOrientation(Qt::Vertical);
	plugin->mToolBar->setWindowTitle("Toolbar_" + plugin->getPluginName());
	plugin->mToolBar->setVisible(false);
	inner->addToolBar(Qt::LeftToolBarArea, plugin->mToolBar);
	// dock

	// put plugin into list
	lplugin << plugin;

	QLabel pLabel(plugin->getPluginName().toUpper() + ": ");
	pLabel.setStyleSheet("QLabel { color: black; font: bold 10px; }");
	plugin->mToolBar->addWidget(&pLabel);
	QLabel pLabelEmpty("");
	plugin->mToolBar->addWidget(&pLabelEmpty);

	Q_FOREACH (QAction* action, plugin->actions())
	{
		// menu
		if (plugin->mMenu)
		{
			if (action)
				plugin->mMenu->addAction(action);
			else
				plugin->mMenu->addSeparator();
		}
		// menu

		// dock
		if (action)
			plugin->mToolBar->addAction(action);
		else
			plugin->mToolBar->addSeparator();
		// dock
	}
	return true;
}

void mainwindow::createActions()
{
	actionNewEmpty = new QAction(this);
	connect(actionNewEmpty, SIGNAL(triggered()), this, SLOT(actionNewEmpty_slot()), Qt::DirectConnection);

	actionAddEmpty = new QAction(this);
	connect(actionAddEmpty, SIGNAL(triggered()), this, SLOT(actionAddEmpty_slot()), Qt::DirectConnection);

	mainwindowActionOpen = new mepp_action();
	connect(mainwindowActionOpen, SIGNAL(sendParamsOpen(QString, QString, int, typeFuncOpenSave)), this, SLOT(actionOpen_slot(QString, QString, int, typeFuncOpenSave)), Qt::DirectConnection);

	mainwindowActionOpen_space = new mepp_action();
	connect(mainwindowActionOpen_space, SIGNAL(sendParams(QString, QString, typeFuncOpenSave)), this, SLOT(actionOpen_space_slot(QString, QString, typeFuncOpenSave)), Qt::DirectConnection);

	mainwindowActionOpen_time = new mepp_action();
	connect(mainwindowActionOpen_time, SIGNAL(sendParams(QString, QString, typeFuncOpenSave)), this, SLOT(actionOpen_time_slot(QString, QString, typeFuncOpenSave)), Qt::DirectConnection);

	mainwindowActionOpen_and_Add_space = new mepp_action();
	connect(mainwindowActionOpen_and_Add_space, SIGNAL(sendParams(QString, QString, typeFuncOpenSave)), this, SLOT(actionOpen_and_Add_space_slot(QString, QString, typeFuncOpenSave)), Qt::DirectConnection);

	mainwindowActionOpen_and_Add_time = new mepp_action();
	connect(mainwindowActionOpen_and_Add_time, SIGNAL(sendParams(QString, QString, typeFuncOpenSave)), this, SLOT(actionOpen_and_Add_time_slot(QString, QString, typeFuncOpenSave)), Qt::DirectConnection);

	mainwindowActionSave_As = new mepp_action();
	connect(mainwindowActionSave_As, SIGNAL(sendParams(QString, QString, typeFuncOpenSave)), this, SLOT(actionSave_As_slot(QString, QString, typeFuncOpenSave)), Qt::DirectConnection);

    for (int i = 0; i < MaxRecentFiles; ++i)
	{
        recentFileActs[i] = new QAction(this);
        recentFileActs[i]->setVisible(false);
        connect(recentFileActs[i], SIGNAL(triggered()), this, SLOT(openRecentFile()));
    }
}

void mainwindow::createMenus()
{
    separatorAct_menuFile = menuFile->addSeparator();
    for (int i = 0; i < MaxRecentFiles; ++i)
        menuFile->addAction(recentFileActs[i]);

    updateRecentFileActions();

	connect(menuWindow, SIGNAL(aboutToShow()), this, SLOT(updateWindowMenu()));
}

void mainwindow::createStatusBar()
{
	STATUSBAR_COMPONENTS = new QLabel(this);
	STATUSBAR_VERTICES = new QLabel(this);
	STATUSBAR_FACETS = new QLabel(this);
	STATUSBAR_EDGES = new QLabel(this);
	STATUSBAR_BOUNDARIES = new QLabel(this);
	STATUSBAR_GENUS = new QLabel(this);

	statusBar()->addPermanentWidget(STATUSBAR_COMPONENTS);
	statusBar()->addPermanentWidget(STATUSBAR_VERTICES);
	statusBar()->addPermanentWidget(STATUSBAR_FACETS);
	statusBar()->addPermanentWidget(STATUSBAR_EDGES);
	statusBar()->addPermanentWidget(STATUSBAR_BOUNDARIES);
	statusBar()->addPermanentWidget(STATUSBAR_GENUS);

	statusBar()->showMessage(tr("Ready"));
}

void mainwindow::update_mesh_properties(bool update_component, bool update_boundary)
{
	bool hasMdiChild = (activeMdiChild() != 0);

	QString components = QString(tr(" 0 component "));
	QString vertices = QString(tr(" 0 vertice "));
	QString facets = QString(tr(" 0 facet "));
	QString edges = QString(tr(" 0 edge "));
	QString boundaries = QString(tr(" no boundary "));
	QString genus = QString(tr(" genus 0 "));

	if (hasMdiChild)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		unsigned int m_nb_components = viewer->getScenePtr()->get_polyhedron()->nb_components();
		unsigned int m_nb_boundaries = viewer->getScenePtr()->get_polyhedron()->nb_boundaries();

		if (update_component)
		{
			m_nb_components = viewer->getScenePtr()->get_polyhedron()->calc_nb_components();
		}

		if (update_boundary)
		{
			m_nb_boundaries = viewer->getScenePtr()->get_polyhedron()->calc_nb_boundaries();
		}

		unsigned int v = viewer->getScenePtr()->get_polyhedron()->size_of_vertices();
		unsigned int e = viewer->getScenePtr()->get_polyhedron()->size_of_halfedges()/2;
		unsigned int f = viewer->getScenePtr()->get_polyhedron()->size_of_facets();
		unsigned int g = viewer->getScenePtr()->get_polyhedron()->genus(m_nb_components,v,f,e,m_nb_boundaries);

		if (m_nb_components > 1)
			components = QString(tr(" %1 components ")).arg(m_nb_components);
		else
			components = QString(tr(" %1 component ")).arg(m_nb_components);

		vertices = QString(tr(" %1 vertices ")).arg(v);

		facets = QString(tr(" %1 facets ")).arg(f);

		edges = QString(tr(" %1 edges ")).arg(e);

		if (m_nb_components == 0)
			boundaries = QString(tr(" no boundary "));
		else
			if (m_nb_components == 1)
				boundaries = QString(tr(" 1 boundary "));
			else
				boundaries = QString(tr(" %1 boundaries ")).arg(m_nb_components);

		genus = QString(tr(" genus %1 ")).arg(g);
	}

	STATUSBAR_COMPONENTS->setText(components);
	STATUSBAR_VERTICES->setText(vertices);
	STATUSBAR_FACETS->setText(facets);
	STATUSBAR_EDGES->setText(edges);
	STATUSBAR_BOUNDARIES->setText(boundaries);
	STATUSBAR_GENUS->setText(genus);
}

void mainwindow::setCurrentFile(const QString &fileName)
{
    QSettings settings(ORGANIZATION, APPLICATION);
    QStringList files = settings.value("recentFileList").toStringList();
    files.removeAll(fileName);
    files.prepend(fileName);
    while (files.size() > MaxRecentFiles)
        files.removeLast();

    settings.setValue("recentFileList", files);

	updateRecentFileActions();
}

void mainwindow::updateRecentFileActions()
{
    QSettings settings(ORGANIZATION, APPLICATION);
    QStringList files = settings.value("recentFileList").toStringList();

    int numRecentFiles = qMin(files.size(), (int)MaxRecentFiles);

    for (int i = 0; i < numRecentFiles; ++i)
	{
        QString text = tr("&%1 %2").arg(i + 1).arg(strippedName(files[i]));
        recentFileActs[i]->setText(text);
        recentFileActs[i]->setData(files[i]);
        recentFileActs[i]->setVisible(true);
    }
    for (int j = numRecentFiles; j < MaxRecentFiles; ++j)
        recentFileActs[j]->setVisible(false);

    separatorAct_menuFile->setVisible(numRecentFiles > 0);
}

QString mainwindow::strippedName(const QString &fullFileName)
{
    return QFileInfo(fullFileName).fileName();
}

void mainwindow::writeSettings()
{
	QSettings settings(ORGANIZATION, APPLICATION);

	QString path;
	QFileInfo fileInfo = model->fileInfo(proxyModel->mapToSource(tree->currentIndex()));
	if (fileInfo.isFile())
		path = fileInfo.absolutePath();
	else
		path = fileInfo.absoluteFilePath();
	settings.setValue("treeLocation", path);

	settings.beginGroup("MainWindow");
	settings.setValue("size", size());
	settings.setValue("pos", pos());
	settings.endGroup();
}

void mainwindow::readSettings()
{
	QSettings settings(ORGANIZATION, APPLICATION);

	treeLocation = settings.value("treeLocation", QDir::currentPath()).toString();

	settings.beginGroup("MainWindow");
	resize(settings.value("size", QSize(1024, 768)).toSize());
	move(settings.value("pos", QPoint(200, 200)).toPoint());
	settings.endGroup();
}

void mainwindow::updateMenus()
{
    bool hasMdiChild = (activeMdiChild() != 0);

	Viewer *viewer = NULL;
	int loadType = Normal;

	if (hasMdiChild)
	{
		viewer = (Viewer *)activeMdiChild();
		loadType = viewer->getScenePtr()->get_loadType();
    
		if (loadType==Normal || loadType==Space)
			actionOpen_and_Add_space->setEnabled(true);
		else
			actionOpen_and_Add_space->setEnabled(false);

		if (loadType==Normal || loadType==Time)
			actionOpen_and_Add_time->setEnabled(true);
		else
			actionOpen_and_Add_time->setEnabled(false);
	}
	else
	{
		actionOpen_and_Add_space->setEnabled(false);
		actionOpen_and_Add_time->setEnabled(false);
	}

	actionSave_As->setEnabled(hasMdiChild);
    actionClose->setEnabled(hasMdiChild);

	actionClone->setEnabled(hasMdiChild);

	actionClose_Window->setEnabled(hasMdiChild);
	actionClose_All->setEnabled(hasMdiChild);
	actionTile->setEnabled(hasMdiChild && !actionChange_MDI_View_Mode->isChecked());
	actionCascade->setEnabled(hasMdiChild && !actionChange_MDI_View_Mode->isChecked());
	actionNext->setEnabled(hasMdiChild);
	actionPrevious->setEnabled(hasMdiChild);

	actionChange_MDI_View_Mode->setEnabled(hasMdiChild);
	actionChange_Viewer_Mode_Space_Time->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() != Normal);

	this->setWindowTitle(tr("%1 - %2").arg(MAINWINDOW_TITLE).arg(MEPP_VERSION));

	// rendering options
	actionRender_Point->setEnabled(hasMdiChild);
	actionRender_Line->setEnabled(hasMdiChild);
	actionRender_Fill->setEnabled(hasMdiChild);

	actionSuperimpose_Edges->setEnabled(hasMdiChild);
	actionSuperimpose_Vertices->setEnabled(hasMdiChild);
	actionSuperimpose_Vertices_big->setEnabled(hasMdiChild);

	actionVertex_Color->setEnabled(hasMdiChild);
	actionFace_Color->setEnabled(hasMdiChild);

	actionLighting->setEnabled(hasMdiChild);
	actionSmooth_Shading->setEnabled(hasMdiChild);

	actionAntialiasing->setEnabled(hasMdiChild);
	actionCulling->setEnabled(hasMdiChild);

	if (hasMdiChild)
	{
		if (viewer->getRender_Mode() == GL_POINT)
		{
			actionRender_Point->setChecked(true);
			actionRender_Line->setChecked(false);
			actionRender_Fill->setChecked(false);
		}
		else if (viewer->getRender_Mode() == GL_LINE)
		{
			actionRender_Point->setChecked(false);
			actionRender_Line->setChecked(true);
			actionRender_Fill->setChecked(false);
		}
		else if (viewer->getRender_Mode() == GL_FILL)
		{
			actionRender_Point->setChecked(false);
			actionRender_Line->setChecked(false);
			actionRender_Fill->setChecked(true);
		}

		actionSuperimpose_Edges->setChecked(viewer->getSuperimpose_Edges());
		actionSuperimpose_Vertices->setChecked(viewer->getSuperimpose_Vertices());
		actionSuperimpose_Vertices_big->setChecked(viewer->getSuperimpose_Vertices_big());

		actionVertex_Color->setChecked(viewer->getVertex_Color());
		actionFace_Color->setChecked(viewer->getFace_Color());

		actionLighting->setChecked(viewer->getLighting());
		actionSmooth_Shading->setChecked(viewer->getSmooth_Shading());

		actionAntialiasing->setChecked(viewer->getAntialiasing());
		actionCulling->setChecked(viewer->getCulling());

		// show options
		actionShow_FPS->setChecked(viewer->FPSIsDisplayed());
		actionShow_axis->setChecked(viewer->axisIsDrawn());
		actionShow_grid->setChecked(viewer->gridIsDrawn());

		actionShow_normals->setChecked(viewer->getShowNormals());

		actionBounding_box->setChecked(viewer->getBounding_box());
		actionBounding_box_when_moving->setChecked(viewer->getBounding_box_when_moving());
		// show options

		// view options
		actionCouplingRotations->setChecked(viewer->getCouplingRotations());

		actionVBO->setChecked(viewer->getVBO_mode());
		// view options

		// capture options
		actionScreenshot_sequence->setChecked(viewer->getSave_animation());
		// capture options
	}
	// rendering options

	for (int i=0; i<lplugin.size(); i++)
		lplugin[i]->setListWindow(mdiArea->subWindowList());

	// color options
	actionBackground_color->setEnabled(hasMdiChild);
	actionVertex_color->setEnabled(hasMdiChild);
	actionEdge_color->setEnabled(hasMdiChild);
	actionFace_color->setEnabled(hasMdiChild);

	actionMaterial->setEnabled(hasMdiChild);
	// color options

	// show options
	actionShow_FPS->setEnabled(hasMdiChild);
	actionShow_axis->setEnabled(hasMdiChild);
	actionShow_grid->setEnabled(hasMdiChild);

	actionShow_normals->setEnabled(hasMdiChild);

	actionBounding_box->setEnabled(hasMdiChild);
	actionBounding_box_when_moving->setEnabled(hasMdiChild);
	// show options

	// view options
	actionReset_viewpoint->setEnabled(hasMdiChild);
	actionCopy_viewpoint->setEnabled(hasMdiChild);
	actionPaste_viewpoint->setEnabled(hasMdiChild && viewpointCopied);

	actionCenter_all_objects->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Space);
	actionCouplingRotations->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Space);

	actionVBO->setEnabled(hasMdiChild);
	// view options

	// capture options
	actionScreenshot->setEnabled(hasMdiChild);
	actionScreenshot_sequence->setEnabled(hasMdiChild);

	actionClipboard_screenshot->setEnabled(hasMdiChild);
	// capture options

	// dynamic options
	actionParams->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);

	actionReverse_start->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);
	actionReverse_start_loop->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);
	actionStart->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);
	actionStart_loop->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);
	actionStop->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);

	actionDynFirst->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);
	actionDynPrevious->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);
	actionDynNext->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);
	actionDynLast->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);

	actionDynDelete->setEnabled(hasMdiChild && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time);
	// dynamic options

	// status bar
	update_mesh_properties(false, false);
}

void mainwindow::updateWindowMenu()
{
    menuWindow->clear();

	menuWindow->addAction(actionClose_Window);
	menuWindow->addAction(actionClose_All);
	menuWindow->addSeparator();

	menuWindow->addAction(actionTile);
	menuWindow->addAction(actionCascade);
	menuWindow->addSeparator();

	menuWindow->addAction(actionNext);
	menuWindow->addAction(actionPrevious);
	menuWindow->addSeparator();

	menuWindow->addAction(actionChange_MDI_View_Mode);
	menuWindow->addSeparator();

	menuWindow->addAction(actionChange_Viewer_Mode_Space_Time);

	//

	separatorAct_menuWindow = menuWindow->addSeparator();

    QList<QMdiSubWindow *> windows = mdiArea->subWindowList();
    separatorAct_menuWindow->setVisible(!windows.isEmpty());

    for (int i = 0; i < windows.size(); ++i)
	{
        Viewer *viewer = qobject_cast<Viewer *>(windows.at(i)->widget());

        QString text;
        if (i < 9)
		{
            text = tr("&%1 %2").arg(i + 1)
                               .arg(viewer->userFriendlyCurrentFile());
        }
		else
		{
            text = tr("%1 %2").arg(i + 1)
                              .arg(viewer->userFriendlyCurrentFile());
        }
        QAction *action  = menuWindow->addAction(text);
        action->setCheckable(true);
        action ->setChecked(viewer == activeMdiChild());

        connect(action, SIGNAL(triggered()), windowMapper, SLOT(map()));	// to see
        windowMapper->setMapping(action, windows.at(i));
    }
}

void mainwindow::on_actionNew_triggered()
{
    the_viewer = new Viewer(this, lplugin);
    the_viewer->setScenePtr(ScenePtr(new Scene()));

	int res = the_viewer->getScenePtr()->add_mesh(INTERNAL_MESH, Normal, NULL, the_viewer);
	if (!res)
	{
		actionChange_Viewer_Mode_Space_Time->setText(tr("Change Viewer Mode (Normal)"));

		mdiArea->addSubWindow(the_viewer);
		the_viewer->setWindowTitle(INTERNAL_MESH);
		the_viewer->show();

		statusBar()->showMessage(tr("Mesh created"), 2000);

		for (int i=0; i<lplugin.size(); i++)
			lplugin[i]->setListWindow(mdiArea->subWindowList());
	}
	else
	{
		delete the_viewer;

		if (res != 1)
			QMessageBox::warning(this, APPLICATION, tr("Error creating mesh."));
	}
}

void mainwindow::actionNewEmpty_slot()
{
    the_viewer = new Viewer(this, lplugin);
    the_viewer->setScenePtr(ScenePtr(new Scene()));

	int res = the_viewer->getScenePtr()->add_mesh(EMPTY_MESH, Normal, NULL, the_viewer);
	if (!res)
	{
		actionChange_Viewer_Mode_Space_Time->setText(tr("Change Viewer Mode (Normal)"));

		mdiArea->addSubWindow(the_viewer);
		the_viewer->setWindowTitle(EMPTY_MESH);
		the_viewer->show();

		statusBar()->showMessage(tr("Empty mesh created"), 2000);

		for (int i=0; i<lplugin.size(); i++)
			lplugin[i]->setListWindow(mdiArea->subWindowList());
	}
	else
	{
		delete the_viewer;
		
		if (res != 1)
			QMessageBox::warning(this, APPLICATION, tr("Error creating empty mesh."));
	}
}

int mainwindow::loadFile(const QString &fileName, int loadType, typeFuncOpenSave f)
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	the_viewer = new Viewer(this, lplugin);
    the_viewer->setScenePtr(ScenePtr(new Scene()));

	int res = the_viewer->getScenePtr()->add_mesh(fileName, loadType, f, the_viewer);
	if (!res)
	{
		actionChange_Viewer_Mode_Space_Time->setText(tr("Change Viewer Mode (%1)").arg(the_viewer->getScenePtr()->get_stringLoadType()));
		if (loadType == Time)
			actionChange_Viewer_Mode_Space_Time->setChecked(true);

		mdiArea->addSubWindow(the_viewer);
		the_viewer->setWindowTitle(strippedName(fileName));
		the_viewer->show();

		if (f==NULL)
		{
			setCurrentFile(fileName);
			statusBar()->showMessage(tr("File loaded: %1").arg(strippedName(fileName)), 2000);
		}

		for (int i=0; i<lplugin.size(); i++)
			lplugin[i]->setListWindow(mdiArea->subWindowList());
	}
	else
	{
		delete the_viewer;

		QApplication::restoreOverrideCursor();

		if (res==-1)
			QMessageBox::warning(this, APPLICATION,
									   tr("Error while loading file: \"%1\".").
									   arg(fileName));
	}

	QApplication::restoreOverrideCursor();

	return res;
}

int mainwindow::addFile(Viewer *viewer, const QString &fileName, int loadType, typeFuncOpenSave f)
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	int res = viewer->getScenePtr()->add_mesh(fileName, loadType, f, viewer);
	if (!res)
	{
		actionChange_Viewer_Mode_Space_Time->setText(tr("Change Viewer Mode (%1)").arg(the_viewer->getScenePtr()->get_stringLoadType()));
		if (loadType == Time)
			actionChange_Viewer_Mode_Space_Time->setChecked(true);

		viewer->setWindowTitle(tr("%1 - (%2: %3/%4)")
								.arg(strippedName(fileName))
								.arg(viewer->getScenePtr()->get_stringLoadType())
								.arg(viewer->getScenePtr()->get_current_polyhedron()+1)
								.arg(viewer->getScenePtr()->get_nb_polyhedrons()));
		if (f==NULL)
		{
			setCurrentFile(fileName);
			statusBar()->showMessage(tr("File added: %1").arg(strippedName(fileName)), 2000);
		}

		updateMenus();
	}
	else
	{
		QApplication::restoreOverrideCursor();

		if (res==-1)
			QMessageBox::warning(this, APPLICATION,
									   tr("Error while loading file: \"%1\".").
									   arg(fileName));
	}

	QApplication::restoreOverrideCursor();

	return res;
}

void mainwindow::actionOpen_slot(QString title, QString typeFiles, int loadType, typeFuncOpenSave f)
{
    QStringList files = QFileDialog::getOpenFileNames(this, title,
                                         QDir::currentPath(),
                                         typeFiles);

	QStringList::Iterator it = files.begin();
    while (it != files.end())
	{
		if (loadFile(*it, loadType, f))
				break;

        ++it;
    }

	if (files.size()>1)
		mdiArea->tileSubWindows();
}
void mainwindow::on_actionOpen_triggered()
{
    actionOpen_slot(tr("Open Mesh File(s)"), tr("OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply);;ALL files (*.*)"), Normal, NULL);
}

void mainwindow::actionOpen_space_slot(QString title, QString typeFiles, typeFuncOpenSave f)
{
    QStringList files = QFileDialog::getOpenFileNames(this, title,
                                         QDir::currentPath(),
                                         typeFiles);

	QStringList::Iterator it = files.begin();
    while (it != files.end())
	{
		if (it == files.begin())
		{
			if (loadFile(*it, Normal, f))
				break;
		}
		else
		{
			Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild());

			if (addFile(viewer, *it, Space, f))
				break;
		}

        ++it;
    }
}
void mainwindow::on_actionOpen_space_triggered()
{
	actionOpen_space_slot(tr("Open Mesh File(s) (space)"), tr("OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply)"), NULL);
}

void mainwindow::actionOpen_time_slot(QString title, QString typeFiles, typeFuncOpenSave f)
{
    QStringList files = QFileDialog::getOpenFileNames(this, title,
                                         QDir::currentPath(),
                                         typeFiles);

	QStringList::Iterator it = files.begin();
    while (it != files.end())
	{
		if (it == files.begin())
		{
			if (loadFile(*it, Normal, f))
				break;
		}
		else
		{
			Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild());

			if (addFile(viewer, *it, Time, f))
				break;
		}

        ++it;
    }
}
void mainwindow::on_actionOpen_time_triggered()
{
	actionOpen_time_slot(tr("Open Mesh File(s) (time)"), tr("OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply)"), NULL);
}

void mainwindow::actionOpen_and_Add_space_slot(QString title, QString typeFiles, typeFuncOpenSave f)
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		QStringList files = QFileDialog::getOpenFileNames(this, title,
                                         QDir::currentPath(),
                                         typeFiles);

		QStringList::Iterator it = files.begin();
		while (it != files.end())
		{
			if (addFile(viewer, *it, Space, f))
				break;

			++it;
		}
		viewer->recreateListsAndUpdateGL();
	}
}
void mainwindow::on_actionOpen_and_Add_space_triggered()
{
	actionOpen_and_Add_space_slot(tr("Open and Add Mesh File(s) (space)"), tr("OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply)"), NULL);
}

void mainwindow::actionOpen_and_Add_time_slot(QString title, QString typeFiles, typeFuncOpenSave f)
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		QStringList files = QFileDialog::getOpenFileNames(this, title,
                                         QDir::currentPath(),
                                         typeFiles);

		QStringList::Iterator it = files.begin();
		while (it != files.end())
		{
			if (addFile(viewer, *it, Time, f))
				break;

			++it;
		}
		viewer->recreateListsAndUpdateGL();
	}
}
void mainwindow::on_actionOpen_and_Add_time_triggered()
{
	actionOpen_and_Add_time_slot(tr("Open and Add Mesh File(s) (time)"), tr("OFF Files (*.off);;OBJ files (*.obj);;SMF files (*.smf);;PLY files (*.ply)"), NULL);
}

void mainwindow::actionAddEmpty_slot()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		int loadType = viewer->getScenePtr()->get_loadType();
		if (loadType == Normal)
			loadType = Space;

		int res = viewer->getScenePtr()->add_mesh(EMPTY_MESH, loadType, NULL, viewer);
		if (!res)
		{
			actionChange_Viewer_Mode_Space_Time->setText(tr("Change Viewer Mode (%1)").arg(the_viewer->getScenePtr()->get_stringLoadType()));
			if (loadType == Time)
				actionChange_Viewer_Mode_Space_Time->setChecked(true);

			viewer->setWindowTitle(tr("%1 - (%2: %3/%4)")
											.arg(EMPTY_MESH)
											.arg(viewer->getScenePtr()->get_stringLoadType())
											.arg(viewer->getScenePtr()->get_current_polyhedron()+1)
											.arg(viewer->getScenePtr()->get_nb_polyhedrons()));
			statusBar()->showMessage(tr("Empty mesh added"), 2000);

			updateMenus();
		}
	}
}

void mainwindow::openRecentFile()
{
    QAction *action = qobject_cast<QAction *>(sender());
    if (action)
        loadFile(action->data().toString(), Normal, NULL);
}

void mainwindow::actionSave_As_slot(QString title, QString typeFiles, typeFuncOpenSave f)
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		QString fileName = QFileDialog::getSaveFileName(this, title,
											 QDir::currentPath(),
											 typeFiles);
		if (!fileName.isEmpty())
		{
			int res = viewer->getScenePtr()->save_file(fileName, f, viewer);

			if (res==-1)
				QMessageBox::warning(this, APPLICATION,
									   tr("Error while writing file: \"%1\".").
									   arg(fileName));
			else
				statusBar()->showMessage(tr("File saved: %1").arg(strippedName(fileName)), 2000);
		}
	}
}
void mainwindow::on_actionSave_As_triggered()
{
	actionSave_As_slot(tr("Save Mesh File"), tr("OFF Files (*.off)"), NULL);
}

void mainwindow::on_actionClose_triggered()
{
	mdiArea->closeActiveSubWindow();
}

void mainwindow::on_actionClone_triggered()
{
	QMessageBox::information(this, APPLICATION, tr("Function not yet implemented."));
}

void mainwindow::on_actionClose_Window_triggered()
{
	on_actionClose_triggered();
}

void mainwindow::on_actionClose_All_triggered()
{
	mdiArea->closeAllSubWindows();
}

void mainwindow::on_actionChange_MDI_View_Mode_triggered()
{
	if (actionChange_MDI_View_Mode->isChecked())
	{
		mdiArea->setViewMode(QMdiArea::TabbedView);
		actionChange_MDI_View_Mode->setText(tr("Change MDI View Mode (Tabbed View)"));
	}
	else
	{
		mdiArea->setViewMode(QMdiArea::SubWindowView);
		actionChange_MDI_View_Mode->setText(tr("Change MDI View Mode (SubWindow View)"));
	}

	updateMenus();
}

void mainwindow::on_actionChange_Viewer_Mode_Space_Time_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		if (viewer->getScenePtr()->get_loadType() == Space)
		{
			viewer->getScenePtr()->set_loadType(Time);
			actionChange_Viewer_Mode_Space_Time->setText(tr("Change Viewer Mode (Time)"));
		}
		else
		{
			viewer->getScenePtr()->set_loadType(Space);
			actionChange_Viewer_Mode_Space_Time->setText(tr("Change Viewer Mode (Space)"));
			viewer->getScenePtr()->todoIfModeSpace(viewer);
		}

		updateMenus();

		viewer->setWindowTitle(tr("%1 - (%2: %3/%4)")
								.arg(viewer->getScenePtr()->userFriendlyCurrentFile())
								.arg(viewer->getScenePtr()->get_stringLoadType())
								.arg(viewer->getScenePtr()->get_current_polyhedron()+1)
								.arg(viewer->getScenePtr()->get_nb_polyhedrons()));

		viewer->recreateListsAndUpdateGL();
	}
}

void mainwindow::on_actionExit_triggered()
{
    close();
}

void mainwindow::closeEvent(QCloseEvent *event)
{
    mdiArea->closeAllSubWindows();
    if (activeMdiChild())
	{
        event->ignore();
    }
	else
	{
        writeSettings();
        event->accept();
    }
}

void mainwindow::setActiveSubWindow(QWidget *window)
{
    if (!window)
        return;
    mdiArea->setActiveSubWindow(qobject_cast<QMdiSubWindow *>(window));
}

void mainwindow::on_actionAbout_triggered()
{
	QMessageBox::about(this, tr("About MEPP"),
		tr("<br>"
			"<b>MEPP</b><br>"
			"<br>"
			"3D MEsh Processing Platform<br>"
			"LIRIS M2DISCO (c) 2010<br>"
			"<br>"
			"Martial TOLA<br>"
			"<br>"
			MEPP_VERSION"<br>"
			"<br>"));
}

void mainwindow::on_actionAbout_QGLViewer_triggered()
{
	aboutQGLViewer->aboutQGLViewer();
}

void mainwindow::popupAboutBox(QString title, QString html_resource_name)
{
	QFile about_CGAL(html_resource_name);
	about_CGAL.open(QIODevice::ReadOnly);
	QMessageBox mb(QMessageBox::NoIcon,
				 title,
				 QTextStream(&about_CGAL).readAll(),
				 QMessageBox::Ok,
				 this);
	mb.exec();
}

void mainwindow::on_actionAbout_CGAL_triggered()
{
	QIcon icon(":/logo/Pictures/cgal.png");
	QMessageBox messageBox(QMessageBox::NoIcon,
		tr("About CGAL"),
		tr("<h2>Computational Geometry Algorithms Library</h2>"
           "<p>CGAL provides efficient and reliable geometric algorithms in the form of a C++ library.</p>"
           "<p>For more information visit <a href=\"http://www.cgal.org/\">www.cgal.org</a>.</p>")/*.arg(CGAL_VERSION)*/, QMessageBox::Ok);
	messageBox.setIconPixmap(icon.pixmap(157, 41));
	messageBox.exec();
}

// rendering options
void mainwindow::on_actionRender_Point_triggered()
{
	if (activeMdiChild() != 0)
	{
		((Viewer *)activeMdiChild())->setRender_Point();		
		actionRender_Line->setChecked(false);
		actionRender_Fill->setChecked(false);
	}
}
void mainwindow::on_actionRender_Line_triggered()
{
	if (activeMdiChild() != 0)
	{
		((Viewer *)activeMdiChild())->setRender_Line();
		actionRender_Point->setChecked(false);
		actionRender_Fill->setChecked(false);
	}
}
void mainwindow::on_actionRender_Fill_triggered()
{
	if (activeMdiChild() != 0)
	{
		((Viewer *)activeMdiChild())->setRender_Fill();
		actionRender_Point->setChecked(false);
		actionRender_Line->setChecked(false);
	}
}

void mainwindow::on_actionSuperimpose_Edges_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setSuperimpose_Edges(actionSuperimpose_Edges->isChecked());
}
void mainwindow::on_actionSuperimpose_Vertices_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setSuperimpose_Vertices(actionSuperimpose_Vertices->isChecked());
}
void mainwindow::on_actionSuperimpose_Vertices_big_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setSuperimpose_Vertices_big(actionSuperimpose_Vertices_big->isChecked());
}

void mainwindow::on_actionVertex_Color_triggered()
{
	if (activeMdiChild() != 0)
	{
		((Viewer *)activeMdiChild())->setVertex_Color(actionVertex_Color->isChecked());
		actionFace_Color->setChecked(((Viewer *)activeMdiChild())->getFace_Color());
	}
}
void mainwindow::on_actionFace_Color_triggered()
{
	if (activeMdiChild() != 0)
	{
		((Viewer *)activeMdiChild())->setFace_Color(actionFace_Color->isChecked());
		actionVertex_Color->setChecked(((Viewer *)activeMdiChild())->getVertex_Color());
	}
}

void mainwindow::on_actionLighting_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setLighting(actionLighting->isChecked());
}
void mainwindow::on_actionSmooth_Shading_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setSmooth_Shading(actionSmooth_Shading->isChecked());
}

void mainwindow::on_actionAntialiasing_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setAntialiasing(actionAntialiasing->isChecked());
}
void mainwindow::on_actionCulling_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setCulling(actionCulling->isChecked());
}
// rendering options

// color options
void mainwindow::on_actionBackground_color_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		QColor new_color = QColorDialog::getColor(viewer->getViewerBackgroundColor(), this);
		if (new_color.isValid())
			viewer->setViewerBackgroundColor(new_color);
	}
}
void mainwindow::on_actionVertex_color_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		QColor new_color = QColorDialog::getColor(viewer->getViewerVertexColor(), this);
		if (new_color.isValid())
			viewer->setViewerVertexColor(new_color);
	}
}
void mainwindow::on_actionEdge_color_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		QColor new_color = QColorDialog::getColor(viewer->getViewerEdgeColor(), this);
		if (new_color.isValid())
			viewer->setViewerEdgeColor(new_color);
	}
}
void mainwindow::on_actionFace_color_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		QColor new_color = QColorDialog::getColor(viewer->getViewerFaceColor(), this);
		if (new_color.isValid())
			viewer->setViewerFaceColor(new_color);
	}
}

void mainwindow::on_actionMaterial_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		QStringList items;
		items << tr("Silver") << tr("Gold") << tr("Jade") << tr("Light blue") << tr("Emerald") << tr("Polished silver") << tr("Chrome") << tr("Copper")
				<< tr("Polished gold") << tr("Pewter") << tr("Obsidian") << tr("Black plastic") << tr("Polished bronze") << tr("Polished copper")
				<< tr("Pearl") << tr("Ruby") << tr("Turquoise") << tr("Brass") << tr(" None ");
		items.sort();

		bool ok;
		QString item = QInputDialog::getItem(this, tr("Select Material"), tr("Material:"), items, items.indexOf(viewer->get_material().c_str()), false, &ok);
		if (ok && !item.isEmpty())
		{
			viewer->change_material(item.toStdString());
			viewer->recreateListsAndUpdateGL();
		}
	}
}
// color options

// show options
void mainwindow::on_actionShow_FPS_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->toggleFPSIsDisplayed();
}
void mainwindow::on_actionShow_axis_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->toggleAxisIsDrawn();
}
void mainwindow::on_actionShow_grid_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->toggleGridIsDrawn();
}

void mainwindow::on_actionShow_normals_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setShowNormals(actionShow_normals->isChecked());
}

void mainwindow::on_actionBounding_box_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setBounding_box(actionBounding_box->isChecked());
}
void mainwindow::on_actionBounding_box_when_moving_triggered()
{
	if (activeMdiChild() != 0)
	{
                Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

                if (viewer->getVBO_mode())
		{
			actionBounding_box_when_moving->setChecked(false);
			QMessageBox::information(this, APPLICATION, tr("Sorry, 'bounding box when moving' is not possible with this mode."));
		}

                viewer->setBounding_box_when_moving(actionBounding_box_when_moving->isChecked());
	}
}
// show options

// view options
void mainwindow::on_actionReset_viewpoint_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->resetView();
}

void mainwindow::on_actionCopy_viewpoint_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild());

		copyCameraPosition = viewer->camera()->position();
		copyCameraOrientation = viewer->camera()->orientation();

		viewpointCopied = true;
		actionPaste_viewpoint->setEnabled(viewpointCopied);
	}
}
void mainwindow::on_actionPaste_viewpoint_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild());

		viewer->camera()->setPosition(copyCameraPosition);
		viewer->camera()->setOrientation(copyCameraOrientation);

		viewer->updateGL();
	}
}

void mainwindow::on_actionCenter_all_objects_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Space)
		((Viewer *)activeMdiChild())->centerAllObjects();
}
void mainwindow::on_actionCouplingRotations_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Space)
		((Viewer *)activeMdiChild())->setCouplingRotations(actionCouplingRotations->isChecked());
}

void mainwindow::on_actionVBO_triggered()
{
	if (activeMdiChild() != 0)
		((Viewer *)activeMdiChild())->setVBO_mode(actionVBO->isChecked());
}
// view options

// capture options
void mainwindow::on_actionScreenshot_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		viewer->setSnapshotFileName("");
		viewer->saveSnapshot(false, false);
	}
}
void mainwindow::on_actionScreenshot_sequence_triggered()
{
	if (activeMdiChild() != 0)
	{
		Viewer *viewer = qobject_cast<Viewer *>(activeMdiChild()); // avoid bug under Linux

		viewer->saveAnimation(actionScreenshot_sequence->isChecked());
		actionScreenshot_sequence->setChecked(viewer->getSave_animation());

		//mencoder "mf://*.png" -mf fps=25 -o video.avi -ovc x264 -x264encopts qp=26:subq=5:8x8dct:frameref=2:bframes=3:weight_b:pass=1
		//org --> mencoder video.avi -forceidx -of lavf -ovc lavc -oac mp3lame -lavcopts vcodec=flv:vbitrate=320:autoaspect:abitrate=32 -vf scale=320:-3 -af resample=22050 -o video.flv
		//mencoder video.avi -forceidx -of lavf -ovc lavc -oac mp3lame -lavcopts vcodec=flv:vbitrate=4096:autoaspect:abitrate=32 -vf scale=800:-3 -af resample=22050 -o video.flv

/*
		Explications de la ligne de commande :
		- forceidx permet de forcer la réindéxation du fichier dans le cas où le fichier est désynchronisé (désynchronisation audio/video), cela limite les erreurs sur les fichiers mal encodés.
		- of lavf libavformat se charge de l'endage vidéo.
		- ovc lavc les codecs de libavcodec sont utilisés pour encoder la vidéo.
		- oac mp3lame lame est utilisé comme codec audio.
		- lavcopts permet de spécifier les options des codecs vidéo libavcodec
		vcodec=flv utilisation du codec flv
		vbitrate=320 les videos seront encodées avec un bitrate de 320kb/s
		autoaspect permet de définir automatiquement l'aspect ratio
		abitrate=32 le bitrate audio est fixé à 32kb/s
		- vf permet de crée un filtre video
		scale=320:-3 Mise à l'échelle en fixant la largeur à 320 et en présérvant l'aspect ratio original
		- af permet de crée un filtre audio
		resample=22050 Resamplage de l'audio à 22050 Hz
*/
	}
}

void mainwindow::on_actionClipboard_screenshot_triggered()
{
	if (activeMdiChild() != 0)
#if (QGLVIEWER_VERSION < 0x020303)
		QMessageBox::information(this, APPLICATION, tr("You need libQGLViewer superior to version 2.3.3 to use this function."));
#else
		((Viewer *)activeMdiChild())->snapshotToClipboard();
#endif
}
// capture options

// dynamic options
void mainwindow::on_actionParams_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynParams();
}

void mainwindow::on_actionReverse_start_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynReverseStart();
}
void mainwindow::on_actionReverse_start_loop_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynReverseStartLoop();
}
void mainwindow::on_actionStart_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynStart();
}
void mainwindow::on_actionStart_loop_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynStartLoop();
}
void mainwindow::on_actionStop_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynStop();
}

void mainwindow::on_actionDynFirst_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynFirst();
}
void mainwindow::on_actionDynPrevious_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynPrevious();
}
void mainwindow::on_actionDynNext_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynNext();
}
void mainwindow::on_actionDynLast_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynLast();
}

void mainwindow::on_actionDynDelete_triggered()
{
	if (activeMdiChild() != 0 && ((Viewer *)activeMdiChild())->getScenePtr()->get_loadType() == Time)
		((Viewer *)activeMdiChild())->setDynDelete();
}
// dynamic options
