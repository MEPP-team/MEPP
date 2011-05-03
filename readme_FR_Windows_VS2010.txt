///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Marche à suivre pour Mepp (32 bits) sous Windows avec Visual Studio Express 2010 :
----------------------------------------------------------------------------------


1a) récupérer et installer DAEMON Tools Lite afin de pouvoir installer les images 'iso' des logiciels ci-dessous:
http://download.gforge.liris.cnrs.fr/meppbin/windows/DTLite4402-0131.exe

1b) Seulement si vous êtes sous Windows XP, récupérer et installer Windows Installer 4.5 Redistributable:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/WindowsXP-KB942288-v3-x86.exe (3 Mo)

1c) récupérer et installer Visual Studio Express 2010:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/VS2010ExpressFRA.iso (728 Mo)

2a) récupérer les dépendances du projet MEPP (headers & libs, CMake & CMake-gui):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/MEPP/mepp_prebuilt_binaries_vs2010_x86_v01.rar (584 Mo)
et décompresser l'archive dans le répertoire de votre choix (exemple: C:\dev)

2b) récupérer la mise à jour des dépendances du projet MEPP (support vidéo ffmpeg):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/MEPP/mepp_prebuilt_binaries_vs2010_x86_ffmpeg_addon.rar (9 Mo)
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

3b) installer Graphviz : http://download.gforge.liris.cnrs.fr/meppbin/windows/graphviz-2.26.3.msi

Note: Graphviz est utilisé par Doxygen pour la génération des images des graphes de dépendances. C'est l'outil 'dot.exe' qui est appelé.
Par exemple, si vous avez préalablement installé MATLAB, celui-ci utilise lui aussi un outil 'dot.exe' ce qui posera problème et aura pour conséquence d'avoir des images 'vides'.
Dans ce cas, il faut donc changer l'ordre de votre variable Path et positionner Graphviz avant MATLAB.

3c) redémarrer la machine pour la prise en compte des variables d'environnement ci-dessus

4) télécharger les sources de Mepp:
Par exemple, dans C:/mepp/SVN,
avec TortoiseSVN (SVN Checkout...): https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp/trunk

5) utiliser CMake-gui (dans C:\dev\_cmake-2.8.3.20110118_\bin)
 - renseigner le champ "Where is the source code" avec C:/mepp/SVN/trunk (attention, pas C:/mepp/SVN/trunk/src !)
 - renseigner le champ "Where to build the binaries" avec C:/mepp/SVN/trunk/build
 - cliquer sur Configure (en bas à gauche) et choisir comme Generator "Visual Studio 10"
 - activer/désactiver les composants que vous désirez ou non (premières lignes en haut toujours du type BUILD_component_nomducomposant, exemple: BUILD_component_Curvature)
 - cliquer sur Configure (en bas à gauche) à nouveau
 - cliquer sur Generate (en bas à gauche)
 - ouvrir avec Visual Studio* la solution mepp.sln générée dans C:/mepp/SVN/trunk/build puis compiler Mepp
 - se positionner sur le "sous-projet" mepp, faire un "bouton droit" puis cliquer sur "Définir comme projet de démarrage"
 
Note: attention, par défaut le projet se compile en Debug, à vous de basculer en Release si vous le souhaitez.

6) la documentation de Mepp (à venir...) et de votre composant au format Doxygen se génère également via Visual Studio



* Au premier lancement de Visual Studio 2010, il est conseillé d'activer cette option:
Menu Outils -> Paramètres -> Paramètres avancés

* Avec Visual Studio 2010 Ultimate, pour compiler Mepp en debug il peut être nécessaire de désactiver cette option:
Projet « mepp » -> bouton droit -> Propriétés -> Propriétés de configuration -> Outil Manifeste -> Sortie des commentaires -> Non
 