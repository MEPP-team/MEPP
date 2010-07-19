///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Marche à suivre pour Mepp sous Windows :
----------------------------------------

1) récupérer et installer Visual Studio Express 2008 SP1:
http://download.gforge.liris.cnrs.fr/meppbin/windows/Visual_C++_2008_Express_with_SP1/VS2008ExpressWithSP1FRAx1504731.iso (869 Mo)

2) récupérer les dépendances du projet MEPP (headers & libs, CMake & CMake-gui):
http://download.gforge.liris.cnrs.fr/meppbin/windows/MEPP/mepp_cmake_dep_prebuilt_binaries_win32_v01.rar (376 Mo)
et décompresser l'archive dans le répertoire de votre choix (exemple: C:\dev32)

3) positionner 3 variables d’environnement (menu « Poste de travail » -> bouton droit -> Propriétés -> onglet « Avancé » -> bouton « Variables d’environnement » (en bas) -> puis dans la partie « Variables système » (en bas):
 - bouton « nouveau » : rajouter la variable QTDIR avec comme valeur :
C:\dev32\Qt\4.6.3 si vous avez décompressez le fichier ci-dessus dans C:\dev32
 - bouton « nouveau » : rajouter la variable CGAL_DIR avec comme valeur :
C:\dev32\CGAL-3.6 si vous avez décompressez le fichier ci-dessus dans C:\dev32
 - bouton « modifier » : rajouter au sein (à la fin par exemple) de la variable Path:
;C:\dev32\Qt\4.6.3\bin (attention au ;)

4) télécharger les sources de Mepp:
Par exemple, dans C:/mepp/SVN,
avec TortoiseSVN (SVN Checkout...): https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp/trunk

5) utiliser CMake-gui (dans C:\dev32\_cmake-2.8.1-win32-x86_\bin)
 - renseigner le champ "Where is the source code" avec par exemple C:/mepp/SVN/trunk
 - renseigner le champ "Where to build the binaries" avec par exemple C:/mepp/SVN/trunk/build
 - cliquer sur Configure (en bas à gauche) et choisir comme Generator "Visual Studio 9 2008"
 - activer/désactiver les composants que vous désirez ou non (premières lignes en haut toujours du type BUILD_component_nomducomposant, exemple: BUILD_component_Curvature)
 - cliquer sur Configure (en bas à gauche) à nouveau
 - cliquer sur Generate (en bas à gauche)
 - ouvrir avec Visual Studio la solution mepp.sln générée dans C:/mepp/SVN/trunk/build puis compiler Mepp
 - se positionner sur le "sous-projet" mepp, faire un "bouton droit" puis cliquer sur "Définir comme projet de démarrage"
