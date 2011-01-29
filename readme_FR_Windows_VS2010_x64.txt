///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Marche à suivre pour Mepp (64 bits) sous Windows avec Visual Studio Express 2010 :
----------------------------------------------------------------------------------


Note: vous devez disposer d'une machine 64 bits (Vista 64 bits ou Seven 64 bits)

1a) récupérer et installer Windows Installer 4.5 Redistributable:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/WindowsXP-KB942288-v3-x86.exe (3 Mo)

1b) récupérer et installer Visual Studio Express 2010:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/VS2010ExpressFRA.iso (728 Mo)

1c) récupérer et installer MSDK for Windows 7 and .NET Framework 4:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/GRMSDKX_EN_DVD_7_and_4.0_(7.1).iso (571 Mo)

1d) récupérer et installer un patch pour Visual Studio Express 2010:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/VS10-KB2280741-x86.exe (8 Mo)

2) récupérer les dépendances du projet MEPP (headers & libs, CMake & CMake-gui):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/MEPP/mepp_prebuilt_binaries_vs2010_x64_v01.rar (646 Mo)
et décompresser l'archive dans le répertoire de votre choix (exemple: C:\dev)

3a) positionner 3 variables d’environnement (menu « Poste de travail » -> bouton droit -> Propriétés -> onglet « Avancé » -> bouton « Variables d’environnement » (en bas) -> puis dans la partie « Variables système » (en bas):
 - bouton « nouveau » : rajouter la variable QTDIR avec comme valeur :
C:\dev\qt-4.7.1 si vous avez décompressé le fichier ci-dessus dans C:\dev
 - bouton « nouveau » : rajouter la variable CGAL_DIR avec comme valeur :
C:\dev\CGAL-3.7 si vous avez décompressé le fichier ci-dessus dans C:\dev
 - bouton « modifier » : rajouter au sein (à la fin par exemple) de la variable Path:
;C:\dev\qt-4.7.1\bin (attention au ;)

3b) installer Graphviz : http://download.gforge.liris.cnrs.fr/meppbin/windows/graphviz-2.26.3.msi

3c) redémarrer la machine pour la prise en compte des variables ci-dessus

4) télécharger les sources de Mepp:
Par exemple, dans C:/mepp/SVN,
avec TortoiseSVN (SVN Checkout...): https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp/trunk

5) utiliser CMake-gui (dans C:\dev\_cmake-2.8.3.20110118_\bin)
 - renseigner le champ "Where is the source code" avec par exemple C:/mepp/SVN/trunk
 - renseigner le champ "Where to build the binaries" avec par exemple C:/mepp/SVN/trunk/build
 - cliquer sur Configure (en bas à gauche) et choisir comme Generator "Visual Studio 10 2010 Win64" (attention, ne pas se tromper avec "Visual Studio 10 2010")
 - activer/désactiver les composants que vous désirez ou non (premières lignes en haut toujours du type BUILD_component_nomducomposant, exemple: BUILD_component_Curvature)
 - cliquer sur Configure (en bas à gauche) à nouveau
 - cliquer sur Generate (en bas à gauche)
 - ouvrir avec Visual Studio la solution mepp.sln générée dans C:/mepp/SVN/trunk/build puis compiler Mepp
 - se positionner sur le "sous-projet" mepp, faire un "bouton droit" puis cliquer sur "Définir comme projet de démarrage"
 
6) la documentation de Mepp (à venir...) et de votre composant au format Doxygen se génère également via Visual Studio



* Au premier lancement de Visual Studio 2010, il est conseillé d'activer cette option:
Menu Outils -> Paramètres -> Paramètres avancés
