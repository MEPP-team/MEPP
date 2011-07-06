///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

A) Créer une boîte de dialogue (avec Qt Designer):
--------------------------------------------------

1a) sous Windows, cet outil se trouve dans "C:\dev\qt-x.x.x\bin": designer.exe
1b) sous Linux, cet outil se trouve dans le menu "Applications > Programmation > Qt 4 Designer"
1c) sous Mac, cet outil se trouve via "Finder > Applications > MacPorts > Qt4 > Designer"

2) ouvrir le .ui de votre composant afin de placer vos widgets comme vous le désirez, pour se faire, vous pouvez lire une partie intéressante d'un tutoriel relatif à Qt Designer ici:
http://www.siteduzero.com/tutoriel-3-11360-modeliser-ses-fenetres-avec-qt-designer.html#ss_part_2

Ce dernier présente en effet notamment les notions importantes de "layouts" et "spacers" ainsi que la fenêtre d'édition des propriétés des widgets.

Vous pouvez également ouvrir les .ui des composants suivants afin de vous en inspirer:
 - dialSettings.ui du composants Canonical
 - dialSettings.ui du composants VSA
 - dialSettingsComp.ui* du composant Compression_Valence

Au cours de cette étape il sera surtout IMPORTANT de NOMMER CORRECTEMENT vos widgets à l'aide de l'inspecteur d'objet présent en haut de la colonne de droite de Qt Designer.

3) sauvegarder votre boîte de dialogue


-----


B) Invoquer votre boîte de dialogue et récupérer les valeurs saisies/renseignées au sein du code de votre composant:
--------------------------------------------------------------------------------------------------------------------

Dans trunk\src\components\'Votre_Catégorie'\'Votre_Composant'\src\mepp_component_Votre_Composant_plugin, utiliser un code similaire à celui-ci:

SettingsDialog_Votre_Composant dial;
if (dial.exec() == QDialog::Accepted)
{
	int iteration = dial.Iteration->value() // 'Iteration' étant un nom de widget (cf. étape ci-dessus)
	...


------------------------------------------------------------------------------------------
Note au sujet de l'utilisation de plusieurs boîtes de dialogue au sein de votre composant:
------------------------------------------------------------------------------------------

1) dupliquer les 3 fichiers suivants:
 - dialSettings_Votre_Composant.ui
 - dialSettings_Votre_Composant.cpp
 - dialSettings_Votre_Composant.hxx
 
 en
 
 - dialSettingsFonction_Votre_Composant.ui
 - dialSettingsFonction_Votre_Composant.cpp
 - dialSettingsFonction_Votre_Composant.hxx
 ('Fonction' représentant quelque chose en rapport avec cette nouvelle boîte de dialogue)
 
2) ouvrir "dialSettingsFonction_Votre_Composant.ui" avec Qt Designer afin de renommer la nouvelle boîte de dialogue à l'aide de l'inspecteur d'objet
(ici la propriété 'objectName' du QDialog devient "SettingsFonction" à la place de "Settings")

3) sauvegarder votre boîte de dialogue

4) rechercher et remplacer en respectant la case:
 - HEADER_MEPP_COMPONENT_VOTRE_COMPOSANT_PLUGIN_SETTINGS_H par HEADER_MEPP_COMPONENT_VOTRE_COMPOSANT_PLUGIN_SETTINGS_FONCTION_H au sein de "dialSettingsFonction_Votre_Composant.hxx"
 - ui_dialSettings_Votre_Composant.h par ui_dialSettingsFonction_Votre_Composant.h au sein de "dialSettingsFonction_Votre_Composant.hxx"
 - Ui_Settings par Ui_SettingsFonction au sein de "dialSettingsFonction_Votre_Composant.hxx"
 
 - dialSettings_Votre_Composant.hxx par dialSettingsFonction_Votre_Composant.hxx au sein de "dialSettingsFonction_Votre_Composant.cpp"
 
 - SettingsDialog_Votre_Composant par SettingsDialogFonction_Votre_Composant au sein des 2 fichiers "dialSettingsFonction_Votre_Composant.cpp(.hxx)"

5) dans trunk\src\components\'Votre_Catégorie'\'Votre_Composant'\src\mepp_component_Votre_Composant_plugin, rajouter sous la ligne 11 la ligne suivante:
#include "dialSettingsFonction_Votre_Composant.hxx"

puis utiliser un code similaire à celui-ci pour récupérer les valeurs saisies/renseignées au sein du code de votre composant:

SettingsDialogFonction_Votre_Composant dialFonction;
if (dialFonction.exec() == QDialog::Accepted)
{
	int iteration = dial.Iteration->value()
	...