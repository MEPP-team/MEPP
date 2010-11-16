///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

A) Développer un nouveau composant:
-----------------------------------

1) s’inspirer de CGAL_Example dans trunk\src\components\Examples\CGAL_Example:
un renommage propre et précis tenant compte de la « case » dans les .h, .hxx et .cpp devrait vous permettre d’obtenir très facilement le squelette de votre composant
(todo: cette étape sera détaillée bientôt de manière plus exaustive)

2) le « noyau » de MEPP ne doit en principe pas être modifié

3) votre composant doit avoir un nom en anglais (validé par Florent) et doit être rattaché à une catégorie ci-dessous:

Analysis
Compression
Distance
Remeshing
Segmentation
Tools
Watermarking

ou éventuellement

Examples

a) Ce nom doit être également le nom du sous-dossier (contenant le code source de votre composant) au sein de cette sous-catégorie sur le svn (cf. composant CGAL_Example).
b) Le nom de la classe C++ de ce composant doit être du type VotreComposant_Component (cf. composant CGAL_Example).
c) Le nom de la classe C++ du plugin (qui va de paire avec la classe ci-dessus) de votre composant doit être du type mepp_component_VotreComposant_plugin (cf. composant CGAL_Example).
Important: attention à la "case sensitive" pour ces 3 points.

A partir de là, le script CMake s'occupe de tout, il n'y a pas à toucher une seule ligne de code du noyau de Mepp (ni le mepp_config.h.in, ni le polyhedron_enriched_polyhedron.h, ni le mainwindow.ui pour le menu).

Votre menu se déclare en effet dans le fichier mepp_component_VotreComposant_plugin.hxx en le rattachant à la catégorie évoquée ci-dessus et en déclarant les actions des menus comme ci-dessous (cf. composant CGAL_Example).

-----

public:
        mepp_component_CGAL_Example_plugin() : mepp_component_plugin_interface() {}
        ~mepp_component_CGAL_Example_plugin()
        {
            delete actionStep_1; delete actionStep_2; delete actionStep_3;
            delete actionStep_4;
        }

        void init(mainwindow* mainWindow, QList<QMdiSubWindow *> lw)
        {
            this->mw = mainWindow;
            this->lwindow = lw;
            this->mPluginName = this->metaObject()->className();
           
            // choice: menuTools, menuDistance_Quality_measure, menuAnalysis_Filtering, menuSegmentation, menuRemeshing_Subdivision,
            // menuCompression, menuWatermaking, menuExamples
            mParentMenu = mainWindow->menuExamples;

            // début --- actions ---
            actionStep_1 = new QAction(tr("Triangulate And Random Color Facets"), this);
            if (actionStep_1)
                connect(actionStep_1, SIGNAL(triggered()), this, SLOT(step1()));

            actionStep_2 = new QAction(tr("Create Center Vertex"), this);
            if (actionStep_2)
                connect(actionStep_2, SIGNAL(triggered()), this, SLOT(step2()));

            actionStep_3 = new QAction(tr("Show Black And White Facets"), this);
            if (actionStep_3)
                connect(actionStep_3, SIGNAL(triggered()), this, SLOT(step3()));

            actionStep_4 = new QAction(tr("Draw Connections"), this);
            if (actionStep_4)
                connect(actionStep_4, SIGNAL(triggered()), this, SLOT(step4()));
            // fin --- actions ---
        }

        QList<QAction*> actions() const
        {
            return QList<QAction*>()    << actionStep_1
                                        << actionStep_2
                                        << actionStep_3
                                        << NULL			// menu separator
                                        << actionStep_4;
        }


B) Utiliser le code d'un composant X au sein de votre composant (sans la gestion des boîtes de dialogue du composant X pour le moment):
---------------------------------------------------------------------------------------------------------------------------------------


Note: dans cette exemple, on utilise le composant Curvature dans le composant CGAL_Example


1) dans trunk\src\components\Examples\CGAL_Example\cmake\use_components.txt, renseigner les fichiers *.c/*.cpp dont vous avez besoin avec la fonction set de CMake (cf. dernière ligne)
	set( use_components ../../Analysis/Curvature/src/Curvature_Component.cpp ../../Analysis/Curvature/src/extract_Vpropres.cpp )

2) dans trunk\src\components\Examples\CGAL_Example\mepp_component_CGAL_Example_plugin.cpp, renseigner en haut du fichier les lignes suivantes:
	// we want to use Curvature component
	#include "../../../Analysis/Curvature/src/Curvature_Component.h"
	typedef boost::shared_ptr<Curvature_Component> Curvature_ComponentPtr;
	// we want to use Curvature component
	
3) ensuite, l'appel et l'utilisation de Curvature dans trunk\src\components\Examples\CGAL_Example\mepp_component_CGAL_Example_plugin.cpp, se fait comme ceci, cf. fonction step9():

void mepp_component_CGAL_Example_plugin::step9()
{
	// active viewer
	if (mw->activeMdiChild() != 0)
	{
		Viewer* viewer = (Viewer *)mw->activeMdiChild();
		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();

		QApplication::setOverrideCursor(Qt::WaitCursor);

		// we use CGAL_Example component here (as usual)
			CGAL_Example_ComponentPtr component_ptr = findOrCreateComponentForViewer<CGAL_Example_ComponentPtr, CGAL_Example_Component>(viewer, polyhedron_ptr);
			component_ptr->TriangulateAndRandomColorFacets(polyhedron_ptr);

			component_ptr->set_init(2);
		// we use CGAL_Example component here (as usual)

		// we use Curvature component here
			bool IsGeo;
			double radius;
			Curvature_ComponentPtr component_ptr_curvature = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);
				
			// params
			radius = 0.001;
			IsGeo=true;

			mw->statusBar()->showMessage(tr("Curvature..."));				
			component_ptr_curvature->principal_curvature(polyhedron_ptr,IsGeo,radius);
			mw->statusBar()->showMessage(tr("Curvature is done"));

			component_ptr_curvature->set_init(2);

			component_ptr_curvature->ConstructColorMap(polyhedron_ptr,1);
			viewer->recreateListsAndUpdateGL();
		// we use Curvature component here
	}

	QApplication::restoreOverrideCursor();
}