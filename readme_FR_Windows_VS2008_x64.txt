///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Marche à suivre pour Mepp (64 bits) sous Windows avec Visual Studio Express 2008 SP1 :
--------------------------------------------------------------------------------------


Note: vous devez disposer d'une machine 64 bits (Vista 64 bits ou Seven 64 bits)

1a) récupérer et installer DAEMON Tools Lite afin de pouvoir installer les images 'iso' des logiciels ci-dessous:
http://download.gforge.liris.cnrs.fr/meppbin/windows/DTLite4402-0131.exe

1b) récupérer et installer Visual Studio Express 2008 SP1:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2008/VS2008ExpressWithSP1FRAx1504731.iso (869 Mo)

1c) récupérer et installer MSDK for Windows 7 and .NET Framework 3.5 SP1:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2008/GRMSDKX_EN_DVD_7_and_3.5_SP1_(7.0).iso (1.4 Go)

1d) activer le nouveau SDK:
menu "Démarrer" -> "Tous les programmes" -> "Microsoft Windows SDK v7.0" -> "CMD Shell",
puis taper "WindowsSdkVer.exe -version:v7.0"

---

2a) récupérer les dépendances du projet MEPP (headers & libs, CMake & CMake-gui):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2008/MEPP/mepp_prebuilt_binaries_vs2008_x64_v02.rar (582 Mo)
et décompresser l'archive dans le répertoire de votre choix (exemple: C:\dev64 - Attention, ne pas utiliser un chemin 'exotique' du type "Documents and Settings\$USER$\My ...\My ...")

2b) récupérer la mise à jour des dépendances du projet MEPP (support vidéo ffmpeg):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2008/MEPP/mepp_prebuilt_binaries_vs2008_x64_ffmpeg_addon.rar (8 Mo)
et décompresser l'archive dans le même répertoire que ci-dessus (si besoin, répondre oui pour écraser des fichiers)

2c) patcher la base de registre car par défaut Visual Studio Express 2008 n'autorise pas la compilation 64 bits:
menu "Démarrer" -> "Tous les programmes" -> "Accessoires" -> bouton droit sur "Invite de commandes" puis "Exécuter en tant qu'administrateur" (important, sinon le script suivant échoue),
puis taper "cd C:\dev64\VCE64BIT_WIN7SDK" (si vous avez décompressé le fichier ci-dessus dans C:\dev64),
puis taper "setup_x64.bat"

---

3a) positionner 3 variables d’environnement (menu « Poste de travail » -> bouton droit -> Propriétés -> onglet « Avancé » -> bouton « Variables d’environnement » (en bas) -> puis dans la partie « Variables système » (en bas):
 - bouton « nouveau » : rajouter la variable QTDIR avec comme valeur :
C:\dev64\Qt_4.6.3_x64 si vous avez décompressé le fichier ci-dessus dans C:\dev64
 - bouton « nouveau » : rajouter la variable CGAL_DIR avec comme valeur :
C:\dev64\CGAL-3.6.1_x64 si vous avez décompressé le fichier ci-dessus dans C:\dev64
 - bouton « modifier » : rajouter au sein (à la fin par exemple) de la variable Path:
;C:\dev64\Qt_4.6.3_x64\bin (attention au ;)

Note: il se peut que vous disposiez d'un autre logiciel (exemple: MiKTeX) utilisant une autre version de Qt (et donc des dll incompatibles) ce qui provoquera une erreur au lancement de Mepp.
Dans ce cas, il faut donc changer l'ordre de votre variable Path et positionner C:\dev\qt-x.x.x\bin avant le logiciel en question.

3b) installer Graphviz : http://download.gforge.liris.cnrs.fr/meppbin/windows/graphviz-2.26.3.msi

Note: Graphviz est utilisé par Doxygen pour la génération des images des graphes de dépendances. C'est l'outil 'dot.exe' qui est appelé.
Par exemple, si vous avez préalablement installé MATLAB, celui-ci utilise lui aussi un outil 'dot.exe' ce qui posera problème et aura pour conséquence d'avoir des images 'vides'.
Dans ce cas, il faut donc changer l'ordre de votre variable Path et positionner Graphviz avant MATLAB.

3c) redémarrer la machine pour la prise en compte des variables d'environnement ci-dessus

4) télécharger les sources de Mepp:
Par exemple, dans C:/mepp/SVN,
avec TortoiseSVN (SVN Checkout...): https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp (en prenant bien soin de renseigner "votre_username_gforge" pour les membres du LIRIS) ou
http://scm.gforge.liris.cnrs.fr/public/mepp (pour les utilisateurs anonymes)

5) utiliser CMake-gui (dans C:\dev64\_cmake-2.8.3.20110118_\bin)
 - renseigner le champ "Where is the source code" avec C:/mepp/SVN/trunk (attention, pas C:/mepp/SVN/trunk/src !)
 - renseigner le champ "Where to build the binaries" avec C:/mepp/SVN/trunk/build
 - cliquer sur Configure (en bas à gauche) et choisir comme Generator "Visual Studio 9 2008 Win64" (attention, ne pas se tromper avec "Visual Studio 9 2008")
 - activer/désactiver les composants que vous désirez ou non (premières lignes en haut toujours du type BUILD_component_nomducomposant, exemple: BUILD_component_Curvature)
 - cliquer sur Configure (en bas à gauche) à nouveau
 - cliquer sur Generate (en bas à gauche)
 - ouvrir avec Visual Studio la solution mepp.sln générée dans C:/mepp/SVN/trunk/build puis compiler Mepp
 - se positionner sur le "sous-projet" mepp, faire un "bouton droit" puis cliquer sur "Définir comme projet de démarrage"
 
Note: attention, par défaut le projet se compile en Debug, à vous de basculer en Release si vous le souhaitez.
 
6) la documentation de Mepp (à venir...) et de votre composant au format Doxygen se génère également via Visual Studio
