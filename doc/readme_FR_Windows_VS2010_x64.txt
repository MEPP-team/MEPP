///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011-2012
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Marche à suivre pour Mepp (64 bits) sous Windows avec Visual Studio Express 2010 :
----------------------------------------------------------------------------------


Note: vous devez disposer d'une machine 64 bits (Vista 64 bits ou Seven 64 bits)

1a) récupérer et installer DAEMON Tools Lite afin de pouvoir installer les images 'iso' des logiciels ci-dessous:
http://download.gforge.liris.cnrs.fr/meppbin/windows/utils/ (DTLitexxxx-yyyy.exe)

1b) récupérer et installer Visual Studio Express 2010:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/VS2010ExpressFRA.iso (728 Mo)

1c) récupérer et installer Windows SDK 7.1 (ie MSDK for Windows 7 and .NET Framework 4):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/GRMSDKX_EN_DVD_7_and_4.0_(7.1).iso (571 Mo)
Note: si vous rencontrez un problème d'installation, vous pouvez utiliser à la place l'installeur web: http://www.microsoft.com/en-us/download/details.aspx?id=8279

1d) récupérer et installer Visual Studio Express 2010 SP1 (IMPORTANT):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/sp1/VS2010SP1dvd1.iso (1.5 Go)

1e) récupérer et installer Visual C++ 2010 SP1 Compiler Update for the Windows SDK 7.1 (IMPORTANT):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/sp1/VC-Compiler-KB2519277.exe (121 Mo)

2a) récupérer les dépendances du projet MEPP (headers & libs, CMake & CMake-gui):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/MEPP/mepp_prebuilt_binaries_vs2010_x64_v01.rar (646 Mo)
et décompresser l'archive dans le répertoire de votre choix (exemple: C:\dev - Attention, ne pas utiliser un chemin 'exotique' du type "Documents and Settings\$USER$\My ...\My ...")

2b) récupérer la mise à jour des dépendances du projet MEPP (support vidéo ffmpeg):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/MEPP/mepp_prebuilt_binaries_vs2010_x64_ffmpeg_addon.rar (8 Mo)
et décompresser l'archive dans le même répertoire que ci-dessus (si besoin, répondre oui pour écraser des fichiers)

3a) positionner 3 variables d’environnement (menu « Poste de travail » -> bouton droit -> Propriétés -> onglet « Avancé » -> bouton « Variables d’environnement » (en bas) -> puis dans la partie « Variables système » (en bas):
 - bouton « nouveau » : rajouter la variable QTDIR avec comme valeur :
C:\dev\qt-4.7.1 si vous avez décompressé le fichier ci-dessus dans C:\dev
 - bouton « nouveau » : rajouter la variable CGAL_DIR avec comme valeur :
C:\dev\CGAL-3.7 si vous avez décompressé le fichier ci-dessus dans C:\dev
 - bouton « modifier » : rajouter au sein (à la fin par exemple) de la variable Path:
;C:\dev\qt-4.7.1\bin (attention au ;)

Note: il se peut que vous disposiez d'un autre logiciel (exemple: MiKTeX) utilisant une autre version de Qt (et donc des dll incompatibles) ce qui provoquera une erreur au lancement de Mepp.
Dans ce cas, il faut donc changer l'ordre de votre variable Path et positionner C:\dev\qt-x.x.x\bin avant le logiciel en question.

3b) installer Graphviz : http://download.gforge.liris.cnrs.fr/meppbin/windows/utils/ (graphviz-x.yy.z.msi)

Note: Graphviz est utilisé par Doxygen pour la génération des images des graphes de dépendances. C'est l'outil 'dot.exe' qui est appelé.
Mais, si vous avez préalablement installé MATLAB, celui-ci utilise lui aussi un outil 'dot.exe' ce qui posera problème et aura pour conséquence d'avoir des images 'vides'.
Dans ce cas, il faut donc changer l'ordre de votre variable Path et positionner Graphviz avant MATLAB.

3c) redémarrer la machine pour la prise en compte des variables d'environnement ci-dessus

4) utiliser CMake-gui (dans C:\dev\_cmake-2.8.3.20110118_\bin)
 - renseigner le champ "Where is the source code" avec C:\MEPP.git\MEPP (attention, pas C:\MEPP.git\MEPP\src !)
 - renseigner le champ "Where to build the binaries" avec C:\MEPP.git\MEPP\build
 - cliquer sur Configure (en bas à gauche) et choisir comme 'Generator': "Visual Studio 10 Win64" (attention, ne pas se tromper avec "Visual Studio 10")
 - activer/désactiver les composants que vous désirez ou non (premières lignes en haut toujours du type BUILD_component_nomducomposant, exemple: BUILD_component_Curvature)
 - cliquer sur Configure (en bas à gauche) à nouveau
 - cliquer sur Generate (en bas à gauche)
 - ouvrir avec Visual Studio la solution mepp.sln générée dans C:\MEPP.git\MEPP\build puis compiler Mepp
 - se positionner sur le "sous-projet" mepp, faire un "bouton droit" puis cliquer sur "Définir comme projet de démarrage"
 
Note: attention, par défaut le projet se compile en Debug, à vous de basculer en Release si vous le souhaitez.
 
5) la documentation de Mepp et de votre composant au format Doxygen se génère également via Visual Studio



* Au premier lancement de Visual Studio 2010, il est conseillé d'activer cette option:
Menu "Outils -> Paramètres -> Paramètres avancés"

* Avec Visual Studio 2010 Ultimate, pour compiler Mepp en debug il peut être nécessaire de désactiver cette option:
Projet « mepp » -> bouton droit -> Propriétés -> Propriétés de configuration -> Outil Manifeste -> Sortie des commentaires -> Non
 