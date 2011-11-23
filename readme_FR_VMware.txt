///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Vous trouverez ici http://download.gforge.liris.cnrs.fr/meppbin/vmware/ des machines virtuelles VMware prête à l'emploi avec VMware Player (gratuit) fournit au même endroit.

Il suffit de décompresser la machine de votre choix (format Winrar) puis de démarrer la machine avec VMware Player.


-------------------------------------------------------------------
Marche à suivre pour la machine virtuelle Ubuntu 11.10 / CGAL 3.8 :
-------------------------------------------------------------------

1) Une fois que la machine Linux est active, ouvrez un terminal (menu "Applications > Accessoires > Terminal") et procédez ainsi:
- cd Desktop
- svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp (en prenant bien soin de renseigner "votre_username_gforge" pour les membres du LIRIS) ou
svn checkout http://scm.gforge.liris.cnrs.fr/public/mepp (pour les utilisateurs anonymes)
- puis, pour les membres du LIRIS:
	- choisir "(t)emporaire" en appuyant sur la touche 't'
	- on vous propose ensuite d'entrer le mot de passe pour l'utilisateur 'mepp', appuyez alors une fois sur la touche "entrée" pour saisir à la place votre nom d'utilisateur GForge
	- saisissez ensuite votre mot de passe GForge associé
	- vous pouvez cliquer sur "annuler" sur la boite de dialogue relative au trousseau de clés
	- vous pouvez répondre "non" à la question du stockage en clair du mot de passe

2) Pour compiler, suivez les instructions 5a) ou 5b) du fichier "readme_FR_Linux.txt" présent à la racine du trunk.

---

Si besoin, voici les mots de passe de la machine virtuelle Linux:

Username: mepp
Password: mepp2011

Username: root
Password: mepp2011