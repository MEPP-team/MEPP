///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Vous trouverez ici http://download.gforge.liris.cnrs.fr/meppbin/vmware/ une machine virtuelle VMware prête à l'emploi avec VMware Player (gratuit) fournit au même endroit.

Il suffit de décompresser la machine (format 7z mais se décompresse aussi avec Winrar) puis de démarrer la machine avec VMware Player.

Une fois que la machine Linux est active, ouvrez un terminal (menu "Applications > Accessoires > Terminal") et procédez ainsi:
- cd Desktop
- svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp/trunk

- puis choisir "(t)emporaire" en appuyant sur la touche 't'
- on vous propose ensuite d'entrer le mot de passe pour l'utilisateur 'vadmin', appuyez alors une fois sur la touche "entrée" pour saisir à la place votre nom d'utilisateur GForge
- saisissez ensuite votre mot de passe GForge associé
- vous pouvez cliquer sur "annuler" sur la boite de dialogue relative au trousseau de clés
- vous pouvez répondre "non" à la question du stockage en clair du mot de passe

Pour compiler, suivez les instructions 6a) ou 6b) du fichier "readme_FR_Linux.txt" présent à la racine du trunk.

---

Si besoin, voici les mots de passe de la machine virtuelle Linux:

Username: vadmin
Password: mepp2010

Username: root
Password: mepp2010 