///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////


* Si vous êtes membres permanents ou temporaires (doctorants, stagiaires, ...) du LIRIS :
-----------------------------------------------------------------------------------------
(A) inscrivez vous sur la forge du laboratoire: https://gforge.liris.cnrs.fr/account/register.php
(B) demandez à rejoindre le projet Mepp: https://gforge.liris.cnrs.fr/project/request.php?group_id=39
(C) abonnez-vous aux listes de diffusion: https://gforge.liris.cnrs.fr/mail/?group_id=39
puis
(D) utilisez SVN pour récupérer MEPP


* Si vous êtes un utilisateur anonyme :
---------------------------------------
 utilisez SVN pour récupérer MEPP


---


SVN sous Windows :
------------------

1) installation et configuration de TortoiseSVN (Sub)Version Control:

- télécharger la dernière version en date pour Windows à l’adresse suivante : http://tortoisesvn.net/downloads.html
- installer le logiciel et redémarrer tout de suite (important)

- pour démarrer facilement avec TortoiseSVN, je vous conseille de lire ce petit tutoriel: http://liris.cnrs.fr/martial.tola/MEPP/doc_for_readme/Utiliser%20Subversion%20SVN%20avec%20Tortoise.pdf (source: Jonathan Petitcolas - http://www.jonathan-petitcolas.com/utiliser-subversion-svn-avec-tortoise/)
- pour une utilisation plus approfondie de TortoiseSVN, je vous recommande ce document: http://liris.cnrs.fr/martial.tola/MEPP/doc_for_readme/Tutoriel%20TortoiseSVN.pdf (source: Kevin Fardel - http://kevin.fardel.perso.esil.univmed.fr/documentation/TutorielTortoiseSVN.pdf)

2) récupération du projet MEPP:

- créer un dossier sur votre système qui hébergera votre version de MEPP (ex. : MEPP_SVN) : c’est dans ce dossier que nous allons procéder au "checkout" du projet
- faire un "bouton-droit" sur le dossier créé ci-dessus puis choisir "SVN Checkout" dans le menu contextuel
- dans la fenêtre qui s’ouvre, renseigner le champ "URL of repository" par l’url du projet MEPP sur la forge : https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp
(en prenant bien soin de renseigner "votre_username_gforge" pour les membres du LIRIS) ou http://scm.gforge.liris.cnrs.fr/public/mepp (pour les utilisateurs anonymes) puis valider par "Ok"
- pour les membres du LIRIS:
	- à la première utilisation vous pouvez obtenir une alerte de sécurité relatif à la clef SSH (et/ou à un certificat), répondre "Oui" pour stocker la clef dans le cache
	- une demande d'authentification vous est ensuite demandée
- le "checkout" commence alors


SVN sous Linux :
----------------

1) installation et configuration de SVN (Sub)Version Control:

- utiliser la commande "sudo apt-get install subversion" ou le "gestionnaire de paquets synaptic" pour installer le paquet "subversion"

2) récupération du projet MEPP:

- créer un dossier sur votre système qui hébergera votre version de MEPP (ex. : MEPP_SVN) : c’est dans ce dossier que nous allons procéder au "checkout" du projet
- utiliser la commande "svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp" (en prenant bien soin de renseigner "votre_username_gforge" pour les membres du LIRIS) ou "svn checkout http://scm.gforge.liris.cnrs.fr/public/mepp" (pour les utilisateurs anonymes)
- pour les membres du LIRIS:
	- choisir "(t)emporaire" en appuyant sur la touche 't'
	- on vous propose ensuite d'entrer le mot de passe pour l'utilisateur courant unix, appuyez alors une fois sur la touche "entrée" pour saisir à la place votre nom d'utilisateur GForge
	- saisissez ensuite votre mot de passe GForge associé
	- vous pouvez cliquer sur "annuler" sur la boite de dialogue relative au trousseau de clés
	- vous pouvez répondre "non" à la question du stockage en clair du mot de passe
- le "checkout" commence alors


SVN sous Mac OS X :
-------------------

1) installation et configuration de SVN (Sub)Version Control:

- installer Xcode: http://developer.apple.com/technologies/xcode.html
- installer MacPorts: http://www.macports.org/install.php
- mettre à jour MacPorts: sudo port -v selfupdate
- installer le paquet "subversion" avec MacPorts: sudo port install subversion

2) récupération du projet MEPP:

- créer un dossier sur votre système qui hébergera votre version de MEPP (ex. : MEPP_SVN) : c’est dans ce dossier que nous allons procéder au "checkout" du projet
- utiliser la commande "svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp" (en prenant bien soin de renseigner "votre_username_gforge" pour les membres du LIRIS) ou "svn checkout http://scm.gforge.liris.cnrs.fr/public/mepp" (pour les utilisateurs anonymes)
- pour les membres du LIRIS:
	- choisir "(t)emporaire" en appuyant sur la touche 't'
	- on vous propose ensuite d'entrer le mot de passe pour l'utilisateur courant unix, appuyez alors une fois sur la touche "entrée" pour saisir à la place votre nom d'utilisateur GForge
	- saisissez ensuite votre mot de passe GForge associé
	- vous pouvez répondre "non" à la question du stockage en clair du mot de passe
- le "checkout" commence alors